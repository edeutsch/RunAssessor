#!/usr/bin/env python3
import sys
import argparse
import os.path
import re
import itertools
import json
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
from numpy import exp
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

DEBUG = False

sys.path.append("C:\local\Repositories\GitHub\SpectralLibraries\lib")
from proforma_peptidoform import ProformaPeptidoform

from peptidoform import Peptidoform
from mass_reference import MassReference
from spectrum import Spectrum

# Define column offsets for peak_list. This dict-like behavior is a bit more efficient than using actual dicts
PL_I_PEAK = 0
PL_MZ = 1
PL_INTENSITY = 2
PL_INTERPRETATION_STRING = 3
PL_AGGREGATION_INFO = 4
PL_INTERPRETATIONS = 5
PL_ATTRIBUTES = 6
# Define column offsets for peak_list attributes
PLA_CHARGE = 0
PLA_N_ISOTOPES = 1
PLA_IS_ISOTOPE = 2
PLA_DELTA_PPM = 3
PLA_PARENT_PEAK = 4
PLA_N_NEUTRAL_LOSSES = 5
PLA_IS_NEUTRAL_LOSS = 6
PLA_IS_PRECURSOR = 7
PLA_IS_REPORTER = 8
PLA_DIAGNOSTIC_CATEGORY = 9
PLA_IS_DELETED = 10
# Define column offsets for interpretations
INT_MZ = 0
INT_REFERENCE_PEAK = 1
INT_INTERPRETATION_STRING = 2
INT_DELTA_PPM = 3
INT_SCORE = 4
INT_DELTA_SCORE = 5
INT_COMMONNESS_SCORE = 6
INT_DIAGNOSTIC_CATEGORY = 7



####################################################################################################
#### Define a crude set of expected isotope ratios for peptides of various masses
expected_isotope_ratios = {
    '100': [ 1.0, 0.038 ],
    '200': [ 1.0, 0.108 ],
    '300': [ 1.0, 0.167 ],
    '400': [ 1.0, 0.193 ],
    '500': [ 1.0, 0.230 ],
    '600': [ 1.0, 0.279 ],
    '700': [ 1.0, 0.353 ],
    '800': [ 1.0, 0.434 ],
    '900': [ 1.0, 0.483 ],
    '1000': [ 1.0, 0.483 ],
    '1100': [ 1.0, 0.542 ],
    '1200': [ 1.0, 0.595 ],
    '1300': [ 1.0, 0.658 ],
    '1400': [ 1.0, 0.716 ],
    '1500': [ 1.0, 0.716 ],
    '1600': [ 1.0, 0.819 ],
    '1700': [ 1.0, 0.857 ],
    '1800': [ 1.0, 0.927 ],
    '1900': [ 1.0, 0.986 ],
    '2000': [ 1.0, 1.049 ],
    '2100': [ 1.0, 1.098 ],
    '2200': [ 1.0, 1.172 ],
    '2300': [ 1.0, 1.172 ],
    '2400': [ 1.0, 1.253 ],
    '2500': [ 1.0, 1.301 ],
    '2600': [ 1.0, 1.361 ],
    '2700': [ 1.0, 1.413 ],
    '2800': [ 1.0, 1.413 ],
    '2900': [ 1.0, 1.476 ],
    '3000': [ 1.0, 1.535 ],
    '3100': [ 1.0, 1.638 ],
    '3200': [ 1.0, 1.675 ],
    '3300': [ 1.0, 1.746 ],
    '3400': [ 1.0, 1.805 ],
    '3500': [ 1.0, 1.831 ],
    '3600': [ 1.0, 1.868 ],
    '3700': [ 1.0, 1.916 ],
    '3800': [ 1.0, 1.990 ],
    '3900': [ 1.0, 2.035 ],
    '4000': [ 1.0, 2.112 ],
    '4100': [ 1.0, 2.182 ],
    '4200': [ 1.0, 2.182 ],
    '4300': [ 1.0, 2.249 ],
    '4400': [ 1.0, 2.352 ],
    '4500': [ 1.0, 2.352 ],
    '4600': [ 1.0, 2.480 ],
    '4700': [ 1.0, 2.518 ],
    '4800': [ 1.0, 2.588 ],
    '4900': [ 1.0, 2.647 ],
    '5000': [ 1.0, 2.673 ],
    '5100': [ 1.0, 2.710 ],
    '5200': [ 1.0, 2.758 ],
    '5300': [ 1.0, 2.832 ],
    '5400': [ 1.0, 2.914 ],
    '5500': [ 1.0, 2.914 ],
    '5600': [ 1.0, 2.962 ],
    '5700': [ 1.0, 3.022 ],
    '5800': [ 1.0, 3.074 ],
    '5900': [ 1.0, 3.137 ],
    '6000': [ 1.0, 3.196 ],
}



####################################################################################################
#### A crude lookup for the expected ratio between the monoisotope and the first isotopic peak
def get_mono_to_first_isotope_peak_ratio(mz):

    mz100 = round(mz/100)*100
    if mz100 < 100:
        mz100 = 100
    if mz100 > 6000:
        mz100 = 6000
    ratio = expected_isotope_ratios[str(mz100)][1]
    return ratio



####################################################################################################
#### SpectrumExaminer class
class SpectrumExaminer:
    '''
    - annotate()                Annotate a spectrum by calling a series of methods given a peptidoform
    - find_close_predicted_fragments()  For a given observed m/z, find all the potential matching predicted fragments within tolerance
    - index_peaks()             Put all the peaks into a dict keyed by integer mass to make lookups faster
    - find_close_ions()         Find the closest predicted fragment
    - add_interpretation()      Add an interpretation to a peak
    - identify_isotopes()       Loop through a spectrum and identify all the peaks that are isotopes of another
    - identify_precursors()     Identify all the peaks in a spectrum that might be precursor-related and annotate them
    - identify_neutral_losses() (NOT CURRENTLY USED?) Identify all peaks that might be a neutral loss of another
    - identify_low_mass_ions()  Identify all the low-mass ions in a spectrum that aren't specifically part of predicted fragments
    - identify_reporter_ions()  Identify known reporter ions, somewhat independently of whether they should be there or not
    - identify_complement_ions()  Find pairs of peaks that add up to the precursor (e.g. b2 and y6 ions for length 7 peptide)
    '''

    ####################################################################################################
    #### Constructor
    def __init__(self, mass_reference=None, verbose=0):

        # Set verbosity
        if verbose is None:
            verbose = 0
        self.verbose = verbose

        # If the mass_reference has not yet been set up or passed, then create it
        self.mass_reference = None
        if mass_reference is None:
            self.mass_reference = MassReference()
        else:
            self.mass_reference = mass_reference

        # Set up a dict for attributes related to the predicted spectrum
        self.spectrum_attributes = {}
        self.tolerance = 20.0

        # Set up a data structure for residuals
        self.residuals = {
            'absolute': {
                'ppm_deltas': [],
                'median_delta': 0.0,
                'siqr': 4.0,
            },
            'relative': {
                'ppm_deltas': [],
                'median_delta': 0.0,
                'siqr': 4.0,
            },
        }



    ####################################################################################################
    # Put all the peaks into a dict keyed by integer mass to make lookups faster
    def index_peaks(self, spectrum):

        # First clear a possibly existing index
        spectrum.peak_index = {}

        # Loop over all peaks, putting them in an integer bin
        for peak in spectrum.peak_list:
            int_peak_mz = int(peak[PL_MZ])
            if int_peak_mz not in spectrum.peak_index:
                spectrum.peak_index[int_peak_mz] = []
            spectrum.peak_index[int_peak_mz].append(peak)



    ####################################################################################################
    # Find the closest predicted fragment
    def find_close_ions(self, spectrum, search_mz, tolerance):

        # Override the input tolerance with the reference tolerance
        #tolerance = self.stats['mz_calibration_tolerance_ppm']

        debug = False
        if abs(search_mz - 997) < 0.01:
            debug = True
            print(f"***{search_mz}")

        # We will return a list of possible matches
        matches = []

        # Get the integer mass as a dict key
        int_search_mz = int(search_mz)

        look_in_lower_bin_too = False
        look_in_higher_bin_too = False

        if search_mz - int_search_mz < tolerance:
            look_in_lower_bin_too = True
        if 1.0 - ( search_mz - int_search_mz ) < tolerance:
            look_in_higher_bin_too = True

        #### If there's just one bin to look in and it's not even in the index, then we're done
        if not look_in_lower_bin_too and not look_in_higher_bin_too:
            if int_search_mz not in spectrum.peak_index:
                return matches
            else:
                peak_list = spectrum.peak_index[int_search_mz]

        #### But if we need to look in more than one bin, then create an enhanced peak_list
        else:
            peak_list = []
            if look_in_lower_bin_too:
                if int_search_mz - 1 in spectrum.peak_index:
                    peak_list.extend(spectrum.peak_index[int_search_mz - 1])
            if int_search_mz in spectrum.peak_index:
                peak_list.extend(spectrum.peak_index[int_search_mz])
            if look_in_higher_bin_too:
                if int_search_mz + 1 in spectrum.peak_index:
                    peak_list.extend(spectrum.peak_index[int_search_mz + 1])

        # Loop over all the peaks in this bin and add them to matches if they're within tolerance
        for peak in peak_list:
            i_peak = peak[PL_I_PEAK]
            mz = peak[PL_MZ]
            intensity = peak[PL_INTENSITY]
            interpretation_string = peak[PL_INTERPRETATION_STRING]
            delta = mz - search_mz
            delta_ppm = delta / search_mz * 1e6
            if delta_ppm < -1 * tolerance:
                continue
            if delta_ppm > tolerance:
                return matches

            # Compute a delta score based on distance from the search. FIXME
            delta_score = 1.0
            commonness_score = 1.0
            diagnostic_category = 'urk'

            score = commonness_score * delta_score * intensity / 70000.0
            match = [ mz, i_peak, interpretation_string, -1 * delta_ppm, score, delta_score, commonness_score, diagnostic_category ]
            matches.append(match)

        return matches



    ####################################################################################################
    # Add an interpretation to a peak
    def add_interpretation(self, peak, interpretation, diagnostic_category, residual_type=None):

        if peak[PL_INTERPRETATION_STRING] == '?':
            peak[PL_INTERPRETATION_STRING] = ''
        if len(peak[PL_INTERPRETATION_STRING]) > 0:
            peak[PL_INTERPRETATION_STRING] += ', '
        peak[PL_INTERPRETATION_STRING] += interpretation[INT_INTERPRETATION_STRING] + '/' + '{:.1f}'.format(interpretation[INT_DELTA_PPM]) + 'ppm'
        interpretation[INT_DIAGNOSTIC_CATEGORY] = diagnostic_category
        peak[PL_INTERPRETATIONS].append(interpretation)

        # If a residual_type was provided, store the residuals
        #if residual_type is not None:
        if residual_type is not None and peak[PL_INTENSITY] > 1000:
            self.residuals[residual_type]['ppm_deltas'].append(interpretation[INT_DELTA_PPM])



    ####################################################################################################
    # Loop through a spectrum and identify all the peaks that are isotopes of another
    def identify_isotopes(self, spectrum):

        #### Constants
        average_isotope_delta = 1.003355    # This is the official mass delta of carbon 13 over carbon 12 and seems to work best
        max_charge = 4
        debug = False

        #### Get some basic parameters
        n_peaks = spectrum.attributes['number of peaks']
        tolerance = self.tolerance

        # Loop through and identify isotopes
        for i_peak in range(n_peaks-1):

            # If this peak is already an isotope, no need to look further
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_ISOTOPE]:
                continue

            mz = spectrum.peak_list[i_peak][PL_MZ]
            if debug:
                print(f"Analyzing peak {i_peak} at {mz}")

            i_lookahead_peak = i_peak + 1
            i_isotope = 1
            charge = max_charge
            done = False
            pursue_more_isotopes = False
            while not done:
                lookahead_mz = spectrum.peak_list[i_lookahead_peak][PL_MZ]
                diff = lookahead_mz - mz
                delta = diff * charge - i_isotope * average_isotope_delta
                delta_ppm = delta / mz * 1e6
                if debug: print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} has diff={diff}, delta={delta}, delta_ppm={delta_ppm}")
                abs_delta_ppm = abs(delta_ppm)
                if abs_delta_ppm < tolerance:

                    #### Examine the intensities of the two peaks and compute how far off the expected intensity ratio they are
                    peak_intensity = spectrum.peak_list[i_peak][PL_INTENSITY]
                    lookahead_peak_intensity = spectrum.peak_list[i_lookahead_peak][PL_INTENSITY]
                    actual_ratio = lookahead_peak_intensity / peak_intensity
                    expected_ratio = get_mono_to_first_isotope_peak_ratio(mz*charge)
                    surprise_factor = actual_ratio / expected_ratio

                    if debug: print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} matches charge {charge} with abs_delta_ppm={abs_delta_ppm} with actual_ratio={actual_ratio} and expected_ratio={expected_ratio} and surprise_factor={surprise_factor:.3f}")

                    #### If the intensity ratio is not too far off the expected, then record it as an isotope
                    if surprise_factor < 4:
                        spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_CHARGE] = charge
                        spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_N_ISOTOPES] += 1
                        spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_CHARGE] = charge
                        spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_IS_ISOTOPE] = 1
                        spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_DELTA_PPM] = delta_ppm
                        spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_PARENT_PEAK] = i_peak

                        interpretation_string = f"isotope of peak {i_peak}"
                        commonness_score = 75
                        match = [ lookahead_mz, i_peak, interpretation_string, delta_ppm, 1.0, 1.0, commonness_score, 'isotope' ]
                        self.add_interpretation(spectrum.peak_list[i_lookahead_peak],match,diagnostic_category='isotope',residual_type='relative')

                        pursue_more_isotopes = True
                        done = True

                    #### Or if the ratio is too far off, then this can't be an isotope
                    else:
                        if debug: print(f"    ==> Lookahead peak intensity is too far off the expected, ignore this potential association")

                elif charge == max_charge and delta < 0:
                    if debug: print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is too close even for charge {charge} with delta_ppm={delta_ppm}. Next lookahead peak")
                    charge = max_charge
                    i_lookahead_peak += 1
                    if i_lookahead_peak >= n_peaks:
                        done = True
                elif charge == 1 and delta > 0:
                    if debug: print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is too far even for charge {charge} with delta_ppm={delta_ppm}. Move on to the next peak")
                    done = True
                elif charge == 1 and delta < 0:
                    if debug: print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is closer than the next charge {charge} isotope with delta_ppm={delta_ppm}. Next lookahead peak")
                    charge = max_charge
                    i_lookahead_peak += 1
                    if i_lookahead_peak >= n_peaks:
                        done = True
                else:
                    if debug: print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is not charge {charge} isotope with delta_ppm={delta_ppm}")
                    pass

                #### Try the next lower charge
                if not done:
                    charge -= 1
                    if charge == 0:
                        done = True

            #### If we found an isotope at a particular charge, then pursue more at that charge
            if pursue_more_isotopes:
                done = False
                i_lookahead_peak += 1
                if i_lookahead_peak >= n_peaks:
                    done = True
                i_isotope += 1
                while not done:
                    lookahead_mz = spectrum.peak_list[i_lookahead_peak][PL_MZ]
                    diff = lookahead_mz - mz
                    delta = diff * charge - i_isotope * average_isotope_delta
                    delta_ppm = delta / mz * 1e6
                    abs_delta_ppm = abs(delta_ppm)
                    #print(f"  Look ahead at peak {i_lookahead_peak} at {lookahead_mz} to look for match at charge {charge} isotope {i_isotope} with diff={diff}, delta={delta}, abs_delta_ppm={abs_delta_ppm}")
                    if abs_delta_ppm < tolerance:
                        #print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} matches charge {charge} isotope {i_isotope} with abs_delta_ppm={abs_delta_ppm}")
                        spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_N_ISOTOPES] += 1
                        spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_CHARGE] = charge
                        spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_IS_ISOTOPE] = i_isotope
                        spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_DELTA_PPM] = delta_ppm
                        spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_PARENT_PEAK] = i_peak

                        interpretation_string = f"isotope {i_isotope} of peak {i_peak}"
                        commonness_score = 60
                        match = [ lookahead_mz, i_peak, interpretation_string, delta_ppm, 1.0, 1.0, commonness_score, 'isotope' ]
                        self.add_interpretation(spectrum.peak_list[i_lookahead_peak], match, diagnostic_category='isotope', residual_type='relative')

                        i_isotope += 1
                    elif delta < 0:
                        #print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is too close for isotope {i_isotope} with delta_ppm={delta_ppm}. Next lookahead peak")
                        pass
                    elif delta > 0:
                        #print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is too far for isotope {i_isotope} with delta_ppm={delta_ppm}. Done looking.")
                        done = True

                    i_lookahead_peak += 1
                    if i_lookahead_peak >= n_peaks:
                        done = True



    ####################################################################################################
    def delete_isotopes(self, spectrum):

        n_peaks = spectrum.attributes['number of peaks']

        # Loop through and mark isotopes as deleted
        for i_peak in range(n_peaks):

            # If this peak is an isotope then remove it
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_ISOTOPE]:
                spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_DELETED] = True



    ####################################################################################################
    #### Identify all the peaks in a spectrum that might be precursor-related and annotate them
    def identify_precursors(self, spectrum):

        add_artificial_precursor = False


        # If there's no index yet, build one
        if len(spectrum.peak_index) == 0:
            self.index_peaks(spectrum)

        # If there is not a known analyte and precursor, then nothing we can do
        if '1' not in spectrum.analytes:
            return
        if 'precursor_mz' not in spectrum.analytes['1']:
            return
        precursor_mz = spectrum.analytes['1']['precursor_mz']

        # Define some basic parameters
        n_peaks = spectrum.attributes['number of peaks']
        #tolerance = self.stats['mz_calibration_tolerance_ppm']
        tolerance = self.tolerance

        charge = spectrum.analytes['1']['charge state']
        charge_string = ''
        if charge > 1:
            charge_string = f"^{charge}"

        matches = self.find_close_ions(spectrum,precursor_mz,tolerance)
        #print(f"*** {matches}")
        #exit()
        for match in matches:

            i_match_peak = match[INT_REFERENCE_PEAK]
            spectrum.peak_list[i_match_peak][PL_ATTRIBUTES][PLA_IS_PRECURSOR] += 1
            match[INT_INTERPRETATION_STRING] = f"p{charge_string}"
            match[INT_COMMONNESS_SCORE] = 30
            self.add_interpretation(spectrum.peak_list[i_match_peak],match,diagnostic_category='nondiagnostic',residual_type=None)
            #spectrum['attributes']['has_unfragmented_precursor'] += 1
            #spectrum['msrun_attributes']['n_unfragmented_precursors'] += 1
            #spectrum['keep'][match['peak_number']] = 0

        # Otherwise add an artificial one at the end to gather neutral losses more easily
        if add_artificial_precursor is True and len(matches) == 0:
            i_peak = n_peaks
            mz = precursor_mz
            intensity = 0.1                                   # FIXME spectrum['intensity_profile'][1]
            interpretation_string = 'artificial precursor'
            aggregation_info = ''
            interpretation = [ mz, i_peak, interpretation_string, 0.0, 1.0, 1.0, 1.0, 'nondiagnostic' ]
            interpretations = [ interpretation ]
            attributes = [ charge, 0, 0, 0, -1, 0, 0, 1, 0, 'nondiagnostic' ]
            spectrum.peak_list.append( [ i_peak, mz, intensity, interpretation_string, aggregation_info, interpretations, attributes ] )
            spectrum.attributes['number of peaks'] += 1

        # Loop over all peaks looking for peaks in the isolation window and exclude them
        target_mz = precursor_mz
        lower_offset = 1.5
        upper_offset = 1.5
        if 'isolation window target m/z' in spectrum.attributes:
            target_mz = spectrum.attributes['isolation window target m/z']
        if 'isolation window lower offset' in spectrum.attributes:
            lower_offset = spectrum.attributes['isolation window lower offset']
        if 'isolation window upper offset' in spectrum.attributes:
            upper_offset = spectrum.attributes['isolation window upper offset']
        
        lower_bound = target_mz - lower_offset
        upper_bound = target_mz + upper_offset
        other_precursor_count = 2
        for i_peak in range(n_peaks):
            mz = spectrum.peak_list[i_peak][PL_MZ]
            if mz < lower_bound:
                continue
            if mz > upper_bound:
                break
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_ISOTOPE] > 0:
                continue

            # If is it already known to be a precursor, then skip it
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_PRECURSOR]:
                continue

            spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_PRECURSOR] += 1
            interpretation_string = f"{other_precursor_count}@p"
            commonness_score = 20
            match = [ mz, i_peak, interpretation_string, 0.0, 1.0, 1.0, commonness_score, 'nondiagnostic' ]
            self.add_interpretation(spectrum.peak_list[i_peak],match,diagnostic_category='nondiagnostic',residual_type=None)
            #spectrum['keep'][i_peak] = 0
            other_precursor_count += 1



    ####################################################################################################
    def delete_precursors(self, spectrum):

        n_peaks = spectrum.attributes['number of peaks']

        # Loop through and mark precursors as deleted
        for i_peak in range(n_peaks):

            # If this peak is an isotope then remove it
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_PRECURSOR]:
                spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_DELETED] = True



    ####################################################################################################
    #### Identify all peaks that might be a neutral loss of another
    def identify_neutral_losses(self, spectrum):

        n_peaks = spectrum.attributes['number of peaks']
        mzs_list = []
        tolerance = self.tolerance

        if len(spectrum.peak_index) == 0:
            self.index_peaks(spectrum)

        # Loop over all peaks looking for neutral losses
        for i_peak in range(n_peaks):

            #### If this peak is already labeled as an isotope, then we can ignore it
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_ISOTOPE]:
                continue

            mz = spectrum.peak_list[i_peak][PL_MZ]
            charge = spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_CHARGE]

            #print(f"Looking for neutral losses for peak with mz={mz}")

            for neutral_loss,neutral_loss_attrs in self.mass_reference.neutral_losses.items():
                delta_mass = neutral_loss_attrs['delta_mass']
                formula = neutral_loss_attrs['formula']
                test_charge = charge or 1 # FIXME shoulnd't we test at least charge 1 and 2?
                test_mz = mz - delta_mass / test_charge
                matches = self.find_close_ions(spectrum,test_mz,tolerance)
                for match in matches:
                    #print(f"  foundmatch at delta_ppm={match['delta_ppm']}")
                    i_match_peak = match[INT_REFERENCE_PEAK]
                    #### If this peak has already been classified as an isotope, then don't overide what we already know. Isotopes take precedence
                    if spectrum.peak_list[i_match_peak][PL_ATTRIBUTES][PLA_IS_ISOTOPE] == 0:
                        spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_CHARGE] = test_charge
                        spectrum.peak_list[i_match_peak][PL_ATTRIBUTES][PLA_CHARGE] = test_charge
                        spectrum.peak_list[i_match_peak][PL_ATTRIBUTES][PLA_IS_NEUTRAL_LOSS] += 1
                        spectrum.peak_list[i_match_peak][PL_ATTRIBUTES][PLA_PARENT_PEAK] = i_peak
                        #spectrum['attributes']['n_neutral_loss_peaks'] += 1
                    spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_N_NEUTRAL_LOSSES] += 1

                    interpretation_string = f"? # NL {mz}-{formula}"
                    match[INT_INTERPRETATION_STRING] = interpretation_string
                    match[INT_COMMONNESS_SCORE ] = 4.0
                    self.add_interpretation(spectrum.peak_list[i_match_peak],match,diagnostic_category='unexplained',residual_type='relative')

                    #if neutral_loss == 'phosphoric acid':
                    #    spectrum['attributes']['n_phospho_loss_peaks'] += 1


        #if spectrum['attributes']['n_phospho_loss_peaks'] > 1:
        #    spectrum['msrun_attributes']['n_phospho_loss_spectra'] += 1



    ####################################################################################################
    #### Identify all the low-mass ions in a spectrum that aren't specifically part of predicted fragments
    def identify_low_mass_ions(self, spectrum, peptidoform=None):

        tolerance = self.tolerance

        #### Create a dict of the residues that the peptide has
        contained_residues = {}
        if peptidoform is not None:
            for residue in list(peptidoform.peptide_sequence):
                contained_residues[residue] = True

        #### If there is not yet a peak index, create one
        if len(spectrum.peak_index) == 0:
            self.index_peaks(spectrum)

        #### Fetch a set of known low-mass ions
        low_mass_ions = self.mass_reference.low_mass_ions

        #### For each possible known low-mass ion, see if there's a match                           #### FIXME arent' there likely to be so many more low-mass ions than peaks? go other way?
        for low_mass_ion_name, low_mass_ion_mz in low_mass_ions.items():

            matches = self.find_close_ions(spectrum, low_mass_ion_mz, tolerance)

            if len(matches) > 0:
                spectrum.attributes['n_identified_peptide_low_mass_ions'] += 1

            for match in matches:
                i_match_peak = match[INT_REFERENCE_PEAK]
                match[INT_DELTA_PPM] = -1 * match[INT_DELTA_PPM]                                                        #### Weird things because I reversed the sign of the delta?????

                #### Special handling for an immonium ions
                if low_mass_ion_name[0] == 'I':
                    # If they correspond to a residue in the peptidoform, they are considered nondiagnostic
                    if low_mass_ion_name[1] in contained_residues:
                        match[INT_INTERPRETATION_STRING] = low_mass_ion_name
                        match[INT_COMMONNESS_SCORE] = 85
                        self.add_interpretation(spectrum.peak_list[i_match_peak], match, diagnostic_category='nondiagnostic', residual_type='absolute')
                    #### But if they don't correspond to a residue in the peptidoform, they are considered contamination and annotated a little differently
                    else:
                        match[INT_INTERPRETATION_STRING] = f"0@{low_mass_ion_name}"
                        match[INT_COMMONNESS_SCORE] = 40
                        self.add_interpretation(spectrum.peak_list[i_match_peak], match, diagnostic_category='contamination', residual_type='absolute')

                #### Special handling for chemical formulas
                elif low_mass_ion_name[0] == 'f':
                    match[INT_INTERPRETATION_STRING] = low_mass_ion_name
                    match[INT_COMMONNESS_SCORE] = 40
                    self.add_interpretation(spectrum.peak_list[i_match_peak], match, diagnostic_category='contamination', residual_type='absolute')

                #### Otherwise record it as a contamination ion
                else:
                    # If the form is already like y1{K} or b2{LL}, then just prepend with 0@
                    if '{' in low_mass_ion_name:
                        match[INT_INTERPRETATION_STRING] = f"0@{low_mass_ion_name}"
                    # Else for more complicated things like names, enclude in curly braces
                    else:
                        match[INT_INTERPRETATION_STRING] = '0@_'+ '{' + f"{low_mass_ion_name}" + '}'
                    match[INT_COMMONNESS_SCORE] = 40
                    self.add_interpretation(spectrum.peak_list[i_match_peak],match,diagnostic_category='contamination',residual_type='absolute')



    ####################################################################################################
    #### Identify known reporter ions, somewhat independently of whether they should be there or not
    def identify_reporter_ions(self, spectrum):

        tolerance = self.tolerance
        if len(spectrum.peak_index) == 0:
            self.index_peaks(spectrum)

        reporter_ions = self.mass_reference.reporter_ions

        # Keep a list of reporter ions that we found to use later for looking at precursor losses
        found_reporter_ions = {}

        for reporter_ion_name,reporter_ion_attributes in reporter_ions.items():

            #print(f"Searching for {reporter_ion_name}")
            matches = self.find_close_ions(spectrum, reporter_ion_attributes['mz'], tolerance)

            if len(matches) > 0:
                spectrum.attributes['n_identified_reporter_ions'] += 1

            for match in matches:

                i_match_peak = match[INT_REFERENCE_PEAK]
                match[INT_DELTA_PPM] = -1 * match[INT_DELTA_PPM]

                match[INT_INTERPRETATION_STRING] = f"r[{reporter_ion_name}]"

                #### Special logic to move neutral losses and gains out of the brackets (e.g. r[TMT6plex+H2O] to r[TMT6plex]+H2O)
                regexp_match = re.match(r'r\[(.+)?([\-\+].+)\]$', match[INT_INTERPRETATION_STRING])
                if regexp_match:
                    match[INT_INTERPRETATION_STRING] = f"r[{regexp_match.group(1)}]{regexp_match.group(2)}"

                match[INT_COMMONNESS_SCORE] = 60
                spectrum.peak_list[i_match_peak][PL_ATTRIBUTES][PLA_IS_REPORTER] += 1
                self.add_interpretation(spectrum.peak_list[i_match_peak], match, diagnostic_category='nondiagnostic', residual_type='absolute')

                # Record that we found it for use later
                found_reporter_ions[reporter_ion_name] = reporter_ions[reporter_ion_name]

        # Now loop through again for the ones that we found looking for precursor losses
        precursor_mz = None
        precursor_charge = None
        try:
            precursor_mz = spectrum.analytes['1']['precursor_mz']
            precursor_charge = spectrum.analytes['1']['charge state']
        except:
            #### Can't go any farther without a precursor mz
            return
        if precursor_charge < 2:
            return

        new_precursor_charge = 1
        precursor_mz = precursor_mz * precursor_charge - self.mass_reference.atomic_masses['proton'] * ( precursor_charge - 1 )

        # Define the list of possble neutral loss combinations to search for
        possible_loss_set_list = [ [ 'carbon monoxide' ], [ 'carbon monoxide', 'ammonia' ] ]

        for reporter_ion_name,reporter_ion_attributes in found_reporter_ions.items():

            for possible_loss_set in possible_loss_set_list:

                loss_name = ''
                loss_mass = 0.0
                for loss in possible_loss_set:
                    loss_name += f"-{self.mass_reference.neutral_losses[loss]['formula']}"
                    loss_mass += self.mass_reference.neutral_losses[loss]['delta_mass']


                #print(f"Searching for p - {reporter_ion_name} {loss_name}")
                search_mz = precursor_mz - reporter_ion_attributes['mz'] - loss_mass
                matches = self.find_close_ions(spectrum, search_mz, tolerance)

                if len(matches) > 0:
                    spectrum.attributes['n_identified_reporter_ions'] += 1

                for match in matches:

                    i_match_peak = match[INT_REFERENCE_PEAK]
                    match[INT_INTERPRETATION_STRING] = f"p-{reporter_ion_name}{loss_name}"
                    match[INT_COMMONNESS_SCORE] = 38
                    spectrum.peak_list[i_match_peak][PL_ATTRIBUTES][PLA_IS_REPORTER] += 1
                    self.add_interpretation(spectrum.peak_list[i_match_peak], match, diagnostic_category='nondiagnostic', residual_type=None)



    ####################################################################################################
    def identify_complement_ions(self, spectrum):

        debug = False

        n_peaks = spectrum.attributes['number of peaks']
        if len(spectrum.peak_index) == 0:
            self.index_peaks(spectrum)

        tolerance = self.tolerance
        precursor_mz = spectrum.analytes['1']['precursor_mz']
        precursor_charge = spectrum.analytes['1']['charge state']
        precursor_mass = precursor_mz * precursor_charge

        # Loop over all peaks looking for complement ions based on the supplied precursor
        for i_peak in range(n_peaks):

            mz = spectrum.peak_list[i_peak][PL_MZ]

            complement_mz = precursor_mass - mz

            if debug:
                print(f"Looking for complement ions for peak with mz={mz}")

            matches = self.find_close_ions(spectrum, complement_mz, tolerance)
            for match in matches:
                i_match_peak = match[INT_REFERENCE_PEAK]
                delta_ppm=match[INT_DELTA_PPM]
                if debug:
                    print(f"  found match at delta_ppm={delta_ppm}")

                # The precursor will be its own self-complement
                if i_match_peak == i_peak:
                    interpretation_string = f"? # unfragmented precursor"
                else:
                    interpretation_string = f"? # complement ion of {i_match_peak}"

                match[INT_INTERPRETATION_STRING] = interpretation_string
                match[INT_COMMONNESS_SCORE ] = 1.0
                self.add_interpretation(spectrum.peak_list[i_peak], match, diagnostic_category='unknown', residual_type='relative')



    ####################################################################################################
    def analyze_mass_defects(self, spectrum, show_plot=False, delete_outliers=False):

        debug = True

        # The allowed deviation from the line for a singly charged fragment ion
        sigma = 0.06

        n_peaks = spectrum.attributes['number of peaks']
        if len(spectrum.peak_index) == 0:
            self.index_peaks(spectrum)

        precursor_mz = spectrum.analytes['1']['precursor_mz']
        precursor_charge = spectrum.analytes['1']['charge state']
        precursor_mass = precursor_mz * precursor_charge
        if precursor_charge == 1:
            precursor2plus = None
        elif precursor_charge == 2:
            precursor2plus = precursor_mz
            # If the integer part of the precursor_ms is even, make it odd to get the 2+ fragment ions line right
            if int(precursor2plus)%2 == 0:
                precursor2plus += 1.0
        else:
            precursor2plus = None

        mzs = [ 0.0 ]
        mass_defects = [ 0.0 ]
        labels = [ 'bottom']

        singley_charged_slope = ( precursor_mass - int(precursor_mass) ) / precursor_mass

        # Loop over all peaks analyzing the mass defects
        for i_peak in range(n_peaks):

            mz = spectrum.peak_list[i_peak][PL_MZ]
            mass_defect = mz - int(mz)
            label = str(i_peak)
            if abs(mass_defect - singley_charged_slope * mz) > sigma:
                label = 'X'
                interpretation_string = f"? # too far from 1+ fragment line"
                match = [ mz, i_peak, interpretation_string, 0.0, 0.0, 0.0, 1.0, 'unknown' ]
                self.add_interpretation(spectrum.peak_list[i_peak], match, diagnostic_category='unknown', residual_type='relative')
                if delete_outliers:
                    spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_DELETED] = True

            mzs.append(mz)
            mass_defects.append(mass_defect)
            labels.append(label)

        mzs.append(precursor_mass)
        mass_defects.append(precursor_mass - int(precursor_mass))
        labels.append('top')

        if show_plot:
            import matplotlib.pyplot as plt

            plt.plot( [mzs[0], mzs[-1] ], [ mass_defects[0], mass_defects[-1] ], color='gray')
            plt.plot( [mzs[0], mzs[-1] ], [ mass_defects[0] + sigma, mass_defects[-1] + sigma ], color='gray')
            plt.plot( [mzs[0], mzs[-1] ], [ mass_defects[0] - sigma, mass_defects[-1] - sigma ], color='gray')
            if precursor2plus is not None:
                plt.plot( [113.084/2, precursor2plus], [ 113.084/2 - int(113.084/2), precursor2plus/2 - int(precursor2plus/2) ], color='black')

            plt.scatter(mzs, mass_defects, color='green')
            plt.title(f"Spectrum mass defect plot")
            plt.xlabel('mz')
            plt.ylabel('mass defect')
            plt.grid(True)

            for i_peak in range(len(mzs)):
                plt.text(mzs[i_peak]+5, mass_defects[i_peak] - 0.02, labels[i_peak], size=8, color='black')

            plt.show()



####################################################################################################
#### Gaussian function used during curve fitting procedures
def gaussian_function(x, a, x0, sigma):
    return a * exp( -(x-x0)**2 / (2*sigma**2) )



####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Class representing a peptidoform')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('--usi', action='store', default=None, help='USI to process')
    argparser.add_argument('--tolerance', action='store', default=None, type=float, help='Tolerance in ppm for annotation')
    argparser.add_argument('--input_json_filename', action='store', default=None, type=float, help='Filename of an input json file')
    argparser.add_argument('--show_all_annotations', action='count', help='If set, show all the potential annotations, not just the final one' )
    argparser.add_argument('--examine', action='count', help='If set, examine the spectrum to see what can be learned' )
    argparser.add_argument('--plot', action='count', help='If set, make a nice figure' )
    params = argparser.parse_args()

    # Set verbose mode
    verbose = params.verbose
    if verbose is None:
        verbose = 0

    #### Flag for showing all annotations
    show_all_annotations = False
    if params.show_all_annotations is not None and params.show_all_annotations > 0:
        show_all_annotations = True

    #### Get the user-supplied USI
    if params.usi is None or params.usi == '':
        print("ERROR: A USI must be supplied with the --usi flag")
        return
    usi_string = params.usi

    #### Parse the USI to get USI metadata
    if usi_string.startswith('mzspec'):
        sys.path.append("C:\local\Repositories\GitHub\PSI\SpectralLibraryFormat\implementations\python\mzlib")
        from universal_spectrum_identifier import UniversalSpectrumIdentifier
        usi = UniversalSpectrumIdentifier(usi_string)
        peptidoform_string = usi.peptidoform_string
        charge = usi.charge
        if verbose:
            print("Parsed information from the USI:")
            print(json.dumps(usi.__dict__,sort_keys=True,indent=2))

    else:
        print("ERROR: USI is malformed: {usi_string}")
        return

    #### Fetch the spectrum
    spectrum = Spectrum()
    spectrum.fetch_spectrum(usi_string)

    #### Need to do this as apparently the peptidoform that comes back from usi is a dict, not an object?
    peptidoform = None
    if peptidoform_string is not None and peptidoform_string != '':
        peptidoform = ProformaPeptidoform(peptidoform_string)


    # Examine the spectrum
    if params.examine:
        examiner = SpectrumExaminer()
        examiner.identify_isotopes(spectrum)
        #examiner.delete_isotopes(spectrum)
        spectrum.compute_spectrum_metrics()
        examiner.index_peaks(spectrum)
        examiner.identify_complement_ions(spectrum)
        examiner.identify_low_mass_ions(spectrum)
        examiner.identify_reporter_ions(spectrum)
        examiner.identify_precursors(spectrum)
        #examiner.delete_precursors(spectrum)
        examiner.identify_neutral_losses(spectrum)
        #examiner.analyze_residuals(spectrum)
        #examiner.rescore_interpretations(spectrum)

        show_plot = False
        if params.plot is not None and params.plot > 0:
            show_plot = True
        examiner.analyze_mass_defects(spectrum, show_plot=show_plot, delete_outliers=True)
        print(spectrum.show(show_all_annotations=show_all_annotations))
        return



#### For command line usage
if __name__ == "__main__": main()
