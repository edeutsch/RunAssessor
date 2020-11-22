#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import re
import itertools
import requests
import json
import copy
import pandas as pd
import numpy

from peptidoform import Peptidoform
from mass_reference import MassReference
from spectrum import Spectrum

sys.path.append("G:\Repositories\GitHub\HUPO-PSI\SpectralLibraryFormat\implementations\python\mzlib")
from universal_spectrum_identifier import UniversalSpectrumIdentifier

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
#### SpectrumAnnotator class
class SpectrumAnnotator:

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

        # Set up a list for the predicted fragments
        self.predicted_fragments_list = []
        self.predicted_fragments_index = {}

        # Set up a dict for attributes related to the predicted spectrum
        self.spectrum_attributes = {}

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
    #### Annotate the spectrum using a variety of different strategies given a peptidoform
    def annotate(self, spectrum, peptidoform_string, charge):

        self.identify_isotopes(spectrum)
        self.identify_low_mass_ions(spectrum)
        self.identify_reporter_ions(spectrum)
        self.identify_neutral_losses(spectrum)
        self.annotate_peptidoform(spectrum, peptidoform_string=peptidoform_string, charge=charge)
        self.identify_precursors(spectrum)
        self.analyze_residuals(spectrum)
        self.rescore_interpretations(spectrum)


    ####################################################################################################
    #### Predict all the fragment ions for the provided peptidoform
    def predict_fragment_ions(self, peptidoform=None, charge=1, fragmentation_type='HCD', skip_internal_fragments=False):

        # Use the passed peptidoform or else the previously supplied one or return None
        if peptidoform is None:
            if self.peptidoform is None:
                return
            else:
                peptidoform = self.peptidoform
        else:
            self.peptidoform = peptidoform

        # Store some spectrum attributes
        self.spectrum_attributes['charge'] = charge
        self.spectrum_attributes['fragmentation_type'] = fragmentation_type

        # Clear out the predicted fragments so that this object can be reused without being recreated
        self.predicted_fragments_list = []
        self.predicted_fragments_index = {}

        # Ensure that there are at least some residues
        if len(peptidoform.residues) == 0:
            return

        # Define the series_list
        series_list = [ 'b', 'y' ]
        if fragmentation_type == 'HCD':
            series_list = [ 'a', 'b', 'y' ]
        else:
            eprint("ERROR: Unrecognized fragmentation type")
            return
        base_series_score = { 'y': 90, 'b': 80, 'a': 70, 'm': 60 }

        # Get handles for some needed reference masses
        masses = self.mass_reference.atomic_masses
        residue_masses = self.mass_reference.nr_aa_masses
        ion_series_attr = self.mass_reference.ion_series_attributes
        neutral_losses = self.mass_reference.neutral_losses
        neutral_losses_by_residue = self.mass_reference.neutral_losses_by_residue
        neutral_losses_by_formula = self.mass_reference.neutral_losses_by_formula
        terminal_modifications = self.mass_reference.terminal_modifications

        # Determine the terminal modification masses
        have_labile_nterm_mod = False
        terminal_mass_modifications = { 'nterm': 0.0, 'cterm': 0.0 }
        if peptidoform.nterm['name'] != '':
            if peptidoform.nterm['name'] in terminal_modifications:
                terminal_mass_modifications['nterm'] = terminal_modifications[peptidoform.nterm['name']]['mass']
                if terminal_modifications[peptidoform.nterm['name']]['is_labile']:
                    have_labile_nterm_mod = True
            else:
                eprint(f"ERROR: Unrecognized nterm mod {peptidoform.nterm['name']}")
        if peptidoform.cterm['name'] != '':
            if peptidoform.cterm['name'] in terminal_modifications:
                terminal_mass_modifications['cterm'] = terminal_modifications[peptidoform.cterm['name']]['mass']
            else:
                eprint(f"ERROR: Unrecognized cterm mod {peptidoform.cterm['name']}")

        # Initialize data for each series
        cumulative_mass = {}
        potential_losses = {}
        cumulative_residues = {}
        for series in series_list:
            cumulative_mass[series] = ion_series_attr[series]['mass'] + terminal_mass_modifications[ion_series_attr[series]['terminus_type']]
            #print(f"Add {terminal_mass_modifications[ion_series_attr[series]['terminus_type']]} to {series}")
            potential_losses[series] = {}
            cumulative_residues[series] = 0

        # Prepare to loop through all residues
        peptide_length = len(peptidoform.residues)
        all_annotations = {}

        # Main loop: iterate through each position, working from both ends simultaneously
        for i_residue in range(peptide_length):

            # Add additional series entries for internal ions
            if (i_residue > 0 and i_residue < peptide_length-2) or (have_labile_nterm_mod is True and i_residue > 0 and i_residue < peptide_length-2):
                series_name = f"m{i_residue+1}"
                series_list.append(series_name)
                series_type = 'm'
                #### Internal fragments don't get terminus
                #cumulative_mass[series_name] = ion_series_attr[series_type]['mass'] + terminal_mass_modifications[ion_series_attr[series_type]['terminus_type']]
                cumulative_mass[series_name] = ion_series_attr[series_type]['mass']
                potential_losses[series_name] = {}
                cumulative_residues[series_name] = 0

                # Also add the very common "a" type ion as a CO neutral loss for internal fragments
                loss_type = 'CO'
                potential_losses[series_name][loss_type] = 1

            #### On the last pass for the whole peptide, just compute the b bion, which will become a p down below
            charge_list = range(1,charge+1)
            if i_residue == peptide_length - 1 and have_labile_nterm_mod is False:
                series_list = [ 'b' ]
                charge_list = [ charge ]


            # Generate fragments for each ion series we expect to find
            for series in series_list:

                series_type = series[0]
                cumulative_residues[series] += 1

                # And for each expected charge
                for i_charge in charge_list:

                    # Get the next residue
                    if ion_series_attr[series_type]['terminus_type'] == 'nterm':
                        residue = peptidoform.residues[i_residue]['string']
                    else:
                        residue = peptidoform.residues[peptide_length - 1 - i_residue]['string']

                    # If this residue is not recognized, this is a serious error
                    if residue not in residue_masses:
                        eprint(f"ERROR: Unrecognized residue '{residue}'")
                        return

                    # Only compute certain things on the first pass
                    if i_charge == 1 or i_residue == peptide_length - 1:
                        # Update the cumulative mass
                        cumulative_mass[series] += residue_masses[residue]

                        #### If this is the precursor pass, also add in the other terminus
                        if i_residue == peptide_length - 1:
                            cumulative_mass[series] += terminal_mass_modifications['cterm']
                            cumulative_mass[series] += ion_series_attr['y']['mass']

                        # See if this residue can yield a neutral loss and store it if so
                        if residue in neutral_losses_by_residue:
                            for loss_type in neutral_losses_by_residue[residue]:
                                #print(loss_type)
                                loss_type_formula = loss_type['formula']
                                if loss_type_formula not in potential_losses[series]:
                                    potential_losses[series][loss_type_formula] = 0
                                potential_losses[series][loss_type_formula] += 1
                                #print(f"Adding an instance of {loss_type_formula}")

                    # Create a list of the possible neutral losses
                    # FIXME, probably should limit water and ammonia losses to 2 max each??
                    losses_list = []
                    for potential_neutral_loss_formula, potential_neutral_loss_number in potential_losses[series].items():
                        for i_loss in range(1,potential_neutral_loss_number+1):
                            losses_list.append(potential_neutral_loss_formula)

                    # Create a list of all the possible combinations of neutral losses (including no loss)
                    all_combinations = []
                    for r in range(len(losses_list) + 1):
                        combinations_object = itertools.combinations(losses_list, r)
                        combinations_list = list(combinations_object)
                        all_combinations += combinations_list

                    # Create the annotations for each combination of losses (including no loss)
                    for potential_neutral_loss_combination in all_combinations:
                        loss_string = ''
                        loss_mass = 0.0
                        potential_neutral_loss_combination_list = list(potential_neutral_loss_combination)
                        potential_neutral_loss_combination_list.sort()
                        for potential_neutral_loss_formula in potential_neutral_loss_combination_list:
                            if potential_neutral_loss_formula is not None:
                                loss_string += f"-{potential_neutral_loss_formula}"
                                loss_mass += neutral_losses_by_formula[potential_neutral_loss_formula]['delta_mass']

                        # Create the default interpretation
                        interpretation = f"{series}{i_residue+1}{loss_string}"

                        #### If this is the final pass for the precursor
                        if i_residue == peptide_length - 1 and series == 'b':
                            interpretation = f"p{loss_string}"

                        # But if this is an internal fragment, create that style
                        if series_type == 'm':
                            if cumulative_residues[series] > 1:
                                interpretation = f"{series}:{i_residue+1}{loss_string}"
                            # Unless it is only of length 1, then skip this entirely because immonium ions are handled staticly
                            else:
                                continue

                        # If the current fragment charge is > 1, add a component for that
                        if i_charge > 1:
                            interpretation += f"^{i_charge}"

                        # Avoid duplicate annotations when different permutations lead to the same thing
                        if interpretation not in all_annotations:
                            all_annotations[interpretation] = 1

                            # Compute the interpretation score
                            interpretation_score = base_series_score[series[0]] - ( i_charge - 1) * 6
                            if loss_string:
                                interpretation_score -= loss_string.count('-') * 8

                            # Compute the final mz and score everything
                            mz = ( cumulative_mass[series] - loss_mass + masses['proton'] * i_charge ) / i_charge
                            self.predicted_fragments_list.append( [ mz, [ [ interpretation, interpretation_score ] ] ] )


        # Print out the resulting fragment list for debugging
        #for frag in self.predicted_fragments_list:
        #    print(frag)
        #exit()


        # Now sort the fragments by mz
        self.predicted_fragments_list.sort(key = lambda x: x[0])

        # Loop through to remove redundancy
        new_fragments_list = []
        previous_four_digit_mz = -1.0
        previous_interpretation = [ [ '?' ] ]
        for fragment in self.predicted_fragments_list:
            interpretation = fragment[1][0]
            four_digit_peak_mz = float(int(fragment[0] * 10000 + 0.5)) / 10000.0
            if four_digit_peak_mz == previous_four_digit_mz:
                #print(f" {four_digit_peak_mz}  {interpretation}  previous_interpretation={previous_interpretation}")
                # If there are two internal fragment ions with the same mz, then only keep the first.
                # FIXME. Maybe we ought to keep them and rescore them based on similar parts later? advanced feature
                if interpretation[0][0] == 'm' and previous_interpretation[0][0] == 'm':
                    continue
                new_fragments_list[-1][1].append(fragment[1][0])
            else:
                new_fragments_list.append(fragment)
            previous_four_digit_mz = four_digit_peak_mz
            previous_interpretation = interpretation

        self.predicted_fragments_list = new_fragments_list

        # Create an index for the predicted fragments to make lookups faster
        for fragment in self.predicted_fragments_list:
            int_peak_mz = int(fragment[0])
            if int_peak_mz not in self.predicted_fragments_index:
                self.predicted_fragments_index[int_peak_mz] = []
            self.predicted_fragments_index[int_peak_mz].append(fragment)


    ####################################################################################################
    #### Annotate the spectrum with the predicted fragments from a supplied peptidoform
    def annotate_peptidoform(self, spectrum, peptidoform_string, charge, skip_internal_fragments=False):

        tolerance = 20.0

        peptidoform = Peptidoform(peptidoform_string)

        # Store the stripped sequence if present
        stripped_sequence = ''
        if peptidoform.stripped_sequence is not None:
            stripped_sequence = peptidoform.stripped_sequence

        self.predict_fragment_ions(peptidoform=peptidoform, charge=charge, fragmentation_type='HCD', skip_internal_fragments=skip_internal_fragments)

        for peak in spectrum.peak_list:
            mz = peak[PL_MZ]

            # Have a look at the previously-annotated immonium ions and if they are for residues that are present here, strip the 0@
            # FIXME This is not going to work for IC[Carbamidomethyl] or IS[Phosho]
            if mz < 300 and len(peak[PL_INTERPRETATIONS]) > 0:
                for interpretation in peak[PL_INTERPRETATIONS]:
                    if interpretation[INT_INTERPRETATION_STRING].startswith('0@I'):
                        residue = interpretation[INT_INTERPRETATION_STRING][3]
                        if residue in stripped_sequence:
                            interpretation[INT_INTERPRETATION_STRING] = interpretation[INT_INTERPRETATION_STRING][2:]

            #print(f"Processing peak at {mz}")
            matches = self.find_close_predicted_fragments(mz,tolerance)
            if matches:
                diagnostic_category = 'diagnostic'
                for match in matches:
                    #peak[PL_INTERPRETATION_STRING] = f"{match[INT_INTERPRETATION_STRING]}/" + '{:.1f}'.format(match[INT_DELTA_PPM]) + 'ppm'
                    #peak[PL_INTERPRETATIONS].append(match)
                    if match[INT_INTERPRETATION_STRING].startswith('p'):
                        peak[PL_ATTRIBUTES][PLA_IS_PRECURSOR] += 1
                        diagnostic_category = 'nondiagnostic'

                    self.add_interpretation(peak,match,diagnostic_category=diagnostic_category,residual_type='absolute')


    ####################################################################################################
    # Find the closest predicted fragment
    def find_close_predicted_fragments(self,search_mz,tolerance):

        # Override the input tolerance with the reference tolerance
        #tolerance = self.stats['mz_calibration_tolerance_ppm']

        # We will return a list of possible matches
        matches = []

        # Get the integer mass as a dict key
        int_search_mz = int(search_mz)
        if int_search_mz not in self.predicted_fragments_index:
            return matches

        # Loop over all the peaks in this bin and add them to matches if they're within tolerance
        for predicted_fragment in self.predicted_fragments_index[int_search_mz]:
            fragment_mz = predicted_fragment[0]
            delta = fragment_mz - search_mz
            delta_ppm = delta / search_mz * 1e6
            if delta_ppm < -1 * tolerance:
                continue
            if delta_ppm > tolerance:
                return matches

            # Compute a delta score based on distance from the search. FIXME
            delta_score = 1.0

            # Loop over all the interpretations and add them to the list
            for interpretation in predicted_fragment[1]:
                interpretation_string = interpretation[0]
                commonness_score = interpretation[1]
                score = commonness_score * delta_score
                match = [ fragment_mz, -1, interpretation_string, delta_ppm, score, delta_score, commonness_score, 'unknown' ]
                matches.append(match)

        return matches


    ####################################################################################################
    # Put all the peaks into a dict by integer mass to make lookups faster
    def index_peaks(self,spectrum):

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
    def find_close_ions(self,spectrum,search_mz,tolerance):

        # Override the input tolerance with the reference tolerance
        #tolerance = self.stats['mz_calibration_tolerance_ppm']

        # We will return a list of possible matches
        matches = []

        # Get the integer mass as a dict key
        int_search_mz = int(search_mz)
        #if int_search_mz == 602:
        #    print(f"Checking 602")
        if int_search_mz not in spectrum.peak_index:
            return matches
        #if int_search_mz == 602:
        #    print(f"Found 602")

        # Loop over all the peaks in this bin and add them to matches if they're within tolerance
        for peak in spectrum.peak_index[int_search_mz]:
            i_peak = peak[PL_I_PEAK]
            mz = peak[PL_MZ]
            intensity = peak[PL_INTENSITY]
            interpretation_string = peak[PL_INTERPRETATION_STRING]
            delta = mz - search_mz
            delta_ppm = delta / search_mz * 1e6
            #print(f"+++ i_peak={i_peak}, delta_ppm={delta_ppm}")
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
            #if int_search_mz == 126:
            #    print(f"   {match}")

        return matches


    ####################################################################################################
    # Put all the peaks into a dict by integer mass to make lookups faster
    def add_interpretation(self,peak,interpretation,diagnostic_category,residual_type=None):

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
    def identify_isotopes(self, spectrum):

        # Constants
        average_isotope_delta = 1.00291
        max_charge = 4
        debug = False

        # Define some basic parameters
        n_peaks = spectrum.attributes['number of peaks']
        #tolerance = self.stats['mz_calibration_tolerance_ppm']
        tolerance = 20

        # Loop through and identify isotopes
        #i_peak = 0
        #while i_peak < n_peaks - 1:
        for i_peak in range(n_peaks-1):

            # If this peak is already an isotope, no need to look further
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_ISOTOPE]:
                continue

            mz = spectrum.peak_list[i_peak][PL_MZ]
            if debug: print(f"Analyzing peak {i_peak} at {mz}")
            i_lookahead_peak = i_peak + 1
            i_isotope = 1
            i_charge = 0
            charge = max_charge
            done = False
            pursue_more_isotopes = False
            while not done:
                lookahead_mz = spectrum.peak_list[i_lookahead_peak][PL_MZ]
                diff = lookahead_mz - mz
                delta = diff * charge - i_isotope * average_isotope_delta
                #delta = diff * charge - i_isotope - ( 1 - average_isotope_delta ) * charge
                delta_ppm = delta / mz * 1e6
                if debug: print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} has diff={diff}, delta={delta}, delta_ppm={delta_ppm}")
                abs_delta_ppm = abs(delta_ppm)
                if abs_delta_ppm < tolerance:
                    if debug: print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} matches charge {charge} with abs_delta_ppm={abs_delta_ppm}")
                    spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_CHARGE] = charge
                    spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_N_ISOTOPES] += 1
                    spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_CHARGE] = charge
                    spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_IS_ISOTOPE] = 1
                    spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_DELTA_PPM] = delta_ppm
                    spectrum.peak_list[i_lookahead_peak][PL_ATTRIBUTES][PLA_PARENT_PEAK] = i_peak

                    interpretation_string = f"isotope of peak {i_peak}"
                    commonness_score = 50
                    match = [ lookahead_mz, i_peak, interpretation_string, delta_ppm, 1.0, 1.0, commonness_score, 'isotope' ]
                    self.add_interpretation(spectrum.peak_list[i_lookahead_peak],match,diagnostic_category='isotope',residual_type='relative')

                    pursue_more_isotopes = True
                    done = True

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

                # Try a lower charge
                if not done:
                    charge -= 1
                    if charge == 0:
                        done = True

            # If we found an isotope at a particular charge, then pursue more at that charge
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
                        commonness_score = 45
                        match = [ lookahead_mz, i_peak, interpretation_string, delta_ppm, 1.0, 1.0, commonness_score, 'isotope' ]
                        self.add_interpretation(spectrum.peak_list[i_lookahead_peak],match,diagnostic_category='isotope',residual_type='relative')

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

            #i_peak += 1


    ####################################################################################################
    def remove_isotopes(self, spectrum):

        n_peaks = spectrum.attributes['number of peaks']
        new_peak_list = []

        # Loop through and remove isotopes
        i_new_peak = 0
        for i_peak in range(n_peaks):

            # If this peak is an isotope then remove it
            if not spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_ISOTOPE]:
                new_peak_list.append(spectrum.peak_list[i_peak])
                new_peak_list[i_new_peak][PL_I_PEAK] = i_new_peak
                i_new_peak += 1

        spectrum.peak_list = new_peak_list
        spectrum.attributes['number of peaks'] = len(spectrum.peak_list)


    ####################################################################################################
    def identify_precursors(self, spectrum):

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
        tolerance = 10

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
        if len(matches) == 0:
            i_peak = n_peaks
            mz = precursor_mz
            intensity = 10000                                   # FIXME spectrum['intensity_profile'][1]
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
    def remove_precursors(self, spectrum):

        n_peaks = spectrum.attributes['number of peaks']
        new_peak_list = []

        # Loop through and remove isotopes
        i_new_peak = 0
        for i_peak in range(n_peaks):

            # If this peak is an isotope then remove it
            if not spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_PRECURSOR]:
                new_peak_list.append(spectrum.peak_list[i_peak])
                new_peak_list[i_new_peak][PL_I_PEAK] = i_new_peak
                i_new_peak += 1

        spectrum.peak_list = new_peak_list
        spectrum.attributes['number of peaks'] = len(spectrum.peak_list)


    ####################################################################################################
    def identify_neutral_losses(self, spectrum):

        n_peaks = spectrum.attributes['number of peaks']
        mzs_list = []
        tolerance = 20

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
    def identify_complement_ions(self, spectrum):

        n_peaks = spectrum.attributes['number of peaks']
        if len(spectrum.peak_index) == 0:
            self.index_peaks(spectrum)

        tolerance = 20
        precursor_mz = spectrum.analytes['1']['precursor_mz']
        precursor_charge = spectrum.analytes['1']['charge state']
        precursor_mass = precursor_mz * precursor_charge

        # Loop over all peaks looking for neutral losses
        for i_peak in range(n_peaks):

            mz = spectrum.peak_list[i_peak][PL_MZ]

            #### Why the fudge factor? I dunno. But it seems to help?
            complement_mz = precursor_mass - mz - 0.0035

            print(f"Looking for complement ions for peak with mz={mz}")

            matches = self.find_close_ions(spectrum,complement_mz,tolerance)
            for match in matches:
                i_match_peak = match[INT_REFERENCE_PEAK]
                delta_ppm=match[INT_DELTA_PPM]
                print(f"  found match at delta_ppm={delta_ppm}")

                # The precursor will be its own self-complement
                if i_match_peak == i_peak:
                    interpretation_string = f"? # unfragmented precursor"
                else:
                    interpretation_string = f"? # complement ion of {i_match_peak}"

                match[INT_INTERPRETATION_STRING] = interpretation_string
                match[INT_COMMONNESS_SCORE ] = 1.0
                self.add_interpretation(spectrum.peak_list[i_peak],match,diagnostic_category='unknown',residual_type='relative')



    ####################################################################################################
    def identify_low_mass_ions(self, spectrum):

        tolerance = 20
        if len(spectrum.peak_index) == 0:
            self.index_peaks(spectrum)

        low_mass_ions = self.mass_reference.low_mass_ions

        for low_mass_ion_name,low_mass_ion_mz in low_mass_ions.items():

            matches = self.find_close_ions(spectrum,low_mass_ion_mz,tolerance)

            if len(matches) > 0:
                spectrum.attributes['n_identified_peptide_low_mass_ions'] += 1

            for match in matches:

                i_match_peak = match[INT_REFERENCE_PEAK]
                match[INT_INTERPRETATION_STRING] = f"0@{low_mass_ion_name}"
                match[INT_COMMONNESS_SCORE] = 40
                self.add_interpretation(spectrum.peak_list[i_match_peak],match,diagnostic_category='contamination',residual_type='absolute')

                # Add this match to the residuals list
                # [ 'type', 'scan_number', 'parent_mz', 'parent_intensity', 'child_mz', 'child_intensity', 'charge', 'delta', 'residual', 'residual_ppm', 'tag' ]
                #row = [ low_mass_ion_name, spectrum['scan_number'], low_mass_ion_mz, spectrum['intensity_array'][match['peak_number']], match['mz'],
                #    spectrum['intensity_array'][match['peak_number']], 1, low_mass_ion_mz - match['mz'], -1, match['delta_ppm'], low_mass_ion_name[:2] ]
                #spectrum['residuals_list'].append(row)


    ####################################################################################################
    def remove_low_mass_ions(self, spectrum):

        n_peaks = spectrum.attributes['number of peaks']
        new_peak_list = []

        # Loop through and remove isotopes
        i_new_peak = 0
        for i_peak in range(n_peaks):

            # If this peak is an external ion that won't be useful for de novo, don't add to the new list
            if spectrum.peak_list[i_peak][PL_INTERPRETATION_STRING].startswith('0@_'):
                pass
            elif spectrum.peak_list[i_peak][PL_INTERPRETATION_STRING].startswith('0@I'):
                pass
            else:
                new_peak_list.append(spectrum.peak_list[i_peak])
                new_peak_list[i_new_peak][PL_I_PEAK] = i_new_peak
                i_new_peak += 1

        spectrum.peak_list = new_peak_list
        spectrum.attributes['number of peaks'] = len(spectrum.peak_list)


    ####################################################################################################
    def identify_reporter_ions(self, spectrum):

        tolerance = 20.0
        if len(spectrum.peak_index) == 0:
            self.index_peaks(spectrum)

        reporter_ions = self.mass_reference.reporter_ions

        # Keep a list of reporter ions that we found to use later for looking at precursor losses
        found_reporter_ions = {}

        for reporter_ion_name,reporter_ion_attributes in reporter_ions.items():

            #print(f"Searching for {reporter_ion_name}")
            matches = self.find_close_ions(spectrum,reporter_ion_attributes['mz'],tolerance)

            if len(matches) > 0:
                spectrum.attributes['n_identified_reporter_ions'] += 1

            for match in matches:

                i_match_peak = match[INT_REFERENCE_PEAK]
                match[INT_INTERPRETATION_STRING] = reporter_ion_name
                match[INT_COMMONNESS_SCORE] = 60
                spectrum.peak_list[i_match_peak][PL_ATTRIBUTES][PLA_IS_REPORTER] += 1
                self.add_interpretation(spectrum.peak_list[i_match_peak],match,diagnostic_category='nondiagnostic',residual_type='absolute')

                # Record that we found it for use later
                found_reporter_ions[reporter_ion_name] = reporter_ions[reporter_ion_name]

                # Add this match to the residuals list
                # [ 'type', 'scan_number', 'parent_mz', 'parent_intensity', 'child_mz', 'child_intensity', 'charge', 'delta', 'residual', 'residual_ppm', 'tag' ]
                #row = [ reporter_ion_name, spectrum['scan_number'], reporter_ion_attributes['mz'], spectrum['intensity_array'][match['peak_number']], match['mz'],
                #    spectrum['intensity_array'][match['peak_number']], 1, reporter_ion_attributes['mz'] - match['mz'], -1, match['delta_ppm'], reporter_ion_name[:3] ]
                #spectrum['residuals_list'].append(row)

        # Now loop through again for the ones that we found looking for precursor losses
        precursor_mz = spectrum.analytes['1']['precursor_mz']
        precursor_charge = spectrum.analytes['1']['charge state']
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
                matches = self.find_close_ions(spectrum,search_mz,tolerance)

                if len(matches) > 0:
                    spectrum.attributes['n_identified_reporter_ions'] += 1

                for match in matches:

                    i_match_peak = match[INT_REFERENCE_PEAK]
                    match[INT_INTERPRETATION_STRING] = f"p-{reporter_ion_name}{loss_name}"
                    match[INT_COMMONNESS_SCORE] = 38
                    spectrum.peak_list[i_match_peak][PL_ATTRIBUTES][PLA_IS_REPORTER] += 1
                    self.add_interpretation(spectrum.peak_list[i_match_peak],match,diagnostic_category='nondiagnostic',residual_type=None)


    ####################################################################################################
    def compute_spectrum_metrics(self, spectrum):

        n_peaks = spectrum.attributes['number of peaks']
        intensity_list = []
        mz_list = []

        # Loop through and remove isotopes
        for i_peak in range(n_peaks):
            mz_list.append(spectrum.peak_list[i_peak][PL_MZ])
            intensity_list.append(spectrum.peak_list[i_peak][PL_INTENSITY])

        # Convert mzs to numpy array and compute metrics
        mz_array = numpy.asarray(mz_list)
        sorted_array = numpy.sort(mz_array)
        spectrum.attributes['minimum mz'] = sorted_array[0]
        spectrum.attributes['maximum mz'] = sorted_array[-1]

        # Convert intensities to numpy array and compute metrics
        intensity_array = numpy.asarray(intensity_list)
        sorted_array = numpy.sort(intensity_array)
        if len(intensity_list) >= 10:
            noise = sorted_array[3]
            signal = sorted_array[-3]
        else:
            index = int(len(intensity_list)/5)
            noise = sorted_array[index]
            signal = sorted_array[-1*index]

        spectrum.attributes['minimum intensity'] = sorted_array[0]
        spectrum.attributes['maximum intensity'] = sorted_array[-1]
        if spectrum.attributes['maximum intensity'] == 0:
            spectrum.attributes['maximum intensity'] = 0.01

        # Prevent division by 0
        if noise == 0:
            noise = 0.05

        # Assume that the smallest peaks are signal to noise of 3.0. An arbitrary guess
        # At least in HCD spectra, the smallest peaks are usually not noise.
        # They're the smallest peaks above the noise that the peak-finding algorithm selected
        noise = noise / 3.0

        signal_to_noise = signal / noise

        spectrum.attributes['noise level'] = noise
        spectrum.attributes['signal level'] = signal
        spectrum.attributes['signal to noise level'] = signal_to_noise


    ####################################################################################################
    def select_denovo_start_peaks(self, spectrum):

        n_peaks = spectrum.attributes['number of peaks']
        temp_peak_list = []
        max_intensity = spectrum.attributes['maximum intensity']
        precursor_mz = spectrum.analytes['1']['precursor_mz']
        minimum_mz = spectrum.attributes['minimum mz']

        first_quartile_mz = ( precursor_mz - minimum_mz ) / 2 + minimum_mz

        # Loop through and remove isotopes
        for i_peak in range(n_peaks):
            mz = spectrum.peak_list[i_peak][PL_MZ]
            intensity = spectrum.peak_list[i_peak][PL_INTENSITY]
            offset = 0.0
            if mz > precursor_mz:
                offset = 10.0
            elif mz > first_quartile_mz:
                offset = 5.0
            #temp_peak_list.append( [ i_peak, mz, intensity / max_intensity + offset ] )
            temp_peak_list.append( [ i_peak, mz, intensity / max_intensity ] )

        return sorted(temp_peak_list,key=lambda x: x[2], reverse=True)














    ####################################################################################################
    def interpret_de_novo(self, spectrum, sequencing_parameters):

        # Get the number of peaks and ensure that there is an index
        n_peaks = spectrum.attributes['number of peaks']

        #### Extract some sequence parameters
        precursor_mz = sequencing_parameters['precursor_mz']
        precursor_charge = sequencing_parameters['precursor_charge']
        tolerance = sequencing_parameters['tolerance']
        fragmentation_type = sequencing_parameters['fragmentation_type']

        verbose = 1

        # Get handles for some needed reference masses
        masses = self.mass_reference.atomic_masses
        residue_masses = self.mass_reference.nr_aa_masses
        ion_series_attr = self.mass_reference.ion_series_attributes
        neutral_losses = self.mass_reference.neutral_losses
        neutral_losses_by_residue = self.mass_reference.neutral_losses_by_residue
        neutral_losses_by_formula = self.mass_reference.neutral_losses_by_formula
        #terminal_modifications = self.mass_reference.terminal_modifications

        # Compute the full precursor molecular weight and store
        precursor_mw = ( precursor_mz - masses['proton'] ) * precursor_charge
        sequencing_parameters['precursor_mw'] = precursor_mw


        # Set up a list of potential termini to consider
        potential_termini = []
        ion_series_list = [ 'b', 'y']
        terminus_notations = { 'nterm': 'n-', 'cterm': '-c' }
        for ion_series in ion_series_list:
            terminus_type = ion_series_attr[ion_series]['terminus_type']
            terminus_mass = ion_series_attr[ion_series]['mass']
            potential_termini.append( {
                'ion_series': ion_series,
                'terminus_type': terminus_type,
                'terminus_name': terminus_notations[terminus_type],
                'terminus_mass': terminus_mass } )
            for terminus_name,terminus_attributes in self.mass_reference.terminal_modifications.items():
                if terminus_attributes['terminus'] == terminus_type:
                    if terminus_type == 'nterm':
                        terminus_name = f"[{terminus_name}]-"
                    if terminus_type == 'cterm':
                        terminus_name = f"-[{terminus_name}]"
                    potential_termini.append( {
                        'ion_series': ion_series,
                        'terminus_type': terminus_type,
                        'terminus_name': terminus_name,
                        'terminus_mass': terminus_mass + terminus_attributes['mass'],
                    } )
        sequencing_parameters['potential_termini'] = potential_termini
        #for potential_terminus in potential_termini:
        #    print(potential_terminus)
        #exit()


        # Create the spectral network and add the nodes
        spectrum.network = { 'nodes': {}, 'edges': {}, 'paths': {} }
        for i_peak in range(n_peaks):
            id = str(i_peak)
            node = {
                'id': id,
                'sort_order': i_peak,
                'type': 'observed',
                'mz': spectrum.peak_list[i_peak][PL_MZ],
                'intensity': spectrum.peak_list[i_peak][PL_INTENSITY],
                'charge': 0.0,
                'mass_defect': spectrum.peak_list[i_peak][PL_MZ] - int(spectrum.peak_list[i_peak][PL_MZ]),
                'upward_edges': {},
                'downward_edges': {},
            }
            spectrum.network['nodes'][id] = node
        spectrum.network['nodes']['bottom'] = {
            'id': 'bottom',
            'sort_order': -1,
            'type': 'terminus',
            'mz': 0.0,
            'intensity': spectrum.attributes['maximum intensity'],
            'charge': 1.0,
            'mass_defect': 0.0,
            'upward_edges': {},
            'downward_edges': {},
        }
        spectrum.network['nodes']['top'] = {
            'id': 'top',
            'sort_order': 1000000,
            'type': 'terminus',
            'mz': precursor_mz * precursor_charge,
            'intensity': spectrum.attributes['maximum intensity'],
            'charge': 1.0,
            'mass_defect': precursor_mz * precursor_charge - int(precursor_mz * precursor_charge),
            'upward_edges': {},
            'downward_edges': {},
        }

        # Create an index for the network
        self.reindex_network(spectrum)


        #### Then add the edges
        for i_peak in range(n_peaks):
            self.link_nodes(spectrum, sequencing_parameters, start_i_peak=i_peak, direction='upward')
            self.link_nodes(spectrum, sequencing_parameters, start_i_peak=i_peak, direction='downward')
        self.link_nodes(spectrum, sequencing_parameters, start_i_peak='bottom', direction='upward')
        self.link_nodes(spectrum, sequencing_parameters, start_i_peak='top', direction='downward')


        # Show some stuff
        if True:
            for node_id,node in spectrum.network['nodes'].items():
                print(f"{node_id}\t{node['mz']}\t{node['intensity']}")
                print(json.dumps(node, indent=2, sort_keys=True))
            self.plot_network(spectrum)
            #exit()

        verbose = 2

        #### Loop over all peaks and show some info
        finished_sequences = []
        active_sequences = []
        followed_peaks = {}
        nodes = list(spectrum.network['nodes'].keys())
        sorted_nodes = sorted(nodes, key=lambda x: spectrum.network['nodes'][x]['sort_order'])
        for node_id in sorted_nodes:
            node = spectrum.network['nodes'][node_id]
            if verbose: print(f"Processing node {node_id}")

            done = False
            active_sequences = [ { 'score': 0.0, 'current_node_id': node_id, 'sequence': [] } ]

            if node_id not in followed_peaks:
                while not done:

                    new_active_sequences = []
                    if verbose == 2:
                        print("==================================================================")
                        print(f"n active_sequences = {len(active_sequences)}")
                        print(f"active_sequences = ", active_sequences)
                    for active_sequence in active_sequences:

                        current_node_id = active_sequence['current_node_id']
                        current_node = spectrum.network['nodes'][current_node_id]
                        current_node_mz = current_node['mz']
                        if verbose == 2:
                            print(f"###### Currently examining node {current_node_id} with active sequence {active_sequence}")
                            print(json.dumps(current_node, indent=2, sort_keys=True))

                        #### If there's nowhere to go, then terminate here
                        if len(current_node['upward_edges']) == 0:
                            active_sequence['sequence'].append( '(+' + '{:.4f}'.format( precursor_mw + masses['proton'] - current_node_mz ) + ')' )
                            finished_sequences.append(active_sequence)
                            if verbose == 2:
                                print(f"  xxxx There are no upward edges, so terminate here at {current_node_mz}")
                                print('        xxxx',active_sequence)

                        #### Loop over all the available upward edges
                        for edge_id in current_node['upward_edges']:
                            edge = spectrum.network['edges'][edge_id]
                            next_node_id = edge['node2_id']
                            if next_node_id == current_node_id:
                                next_node_id = edge['node1_id']
                            next_node = spectrum.network['nodes'][next_node_id]
                            if verbose == 2:
                                print(f"  ++++ Currently examining edge {edge_id}, next_node_id={next_node_id}")
                                print(json.dumps(edge, indent=2, sort_keys=True))

                            new_active_sequence = copy.deepcopy(active_sequence)

                            if len(new_active_sequence['sequence']) == 0:
                                if current_node_id == 'bottom':
                                    new_active_sequence['sequence'].append( '*' + edge['terminus_name'] )
                                    new_active_sequence['sequence'].append( edge['label'] )
                                    new_active_sequence['current_node_id'] = next_node_id
                                    new_active_sequence['score'] += edge['score']
                                    new_active_sequences.append(new_active_sequence)
                                    if verbose == 2: print('    ++',new_active_sequence)
                                else:
                                    new_active_sequence['sequence'].append( '*(+' + '{:.4f}'.format(current_node_mz) + ')' )
                                    new_active_sequence['sequence'].append( edge['label'] )
                                    new_active_sequence['score'] += edge['score']
                                    new_active_sequence['current_node_id'] = next_node_id
                                    new_active_sequences.append(new_active_sequence)
                                    if verbose == 2: print('    ++',new_active_sequence)
                            else:
                                if next_node_id == 'top':
                                    #print(json.dumps(current_node, indent=2, sort_keys=True))
                                    #print(json.dumps(edge, indent=2, sort_keys=True))
                                    #print(json.dumps(next_node, indent=2, sort_keys=True))
                                    #exit()
                                    new_active_sequence['sequence'].append( edge['label'] )
                                    new_active_sequence['sequence'].append( edge['terminus_name'] + 'x' )
                                    new_active_sequence['score'] += edge['score']
                                    finished_sequences.append(new_active_sequence)
                                    if verbose == 2: print('    ====',new_active_sequence)
                                elif len(next_node['upward_edges']) == 0:
                                    new_active_sequence['sequence'].append( '(+' + '{:.4f}'.format( precursor_mw + masses['proton'] - next_node['mz'] ) + ')x' )
                                    finished_sequences.append(new_active_sequence)
                                    if verbose == 2:
                                        print(f"  xxxx There are no upward edges, so terminate here at {current_node_mz}")
                                        print('        xxxx',active_sequence)
                                else:
                                    new_active_sequence['sequence'].append( edge['label'] )
                                    new_active_sequence['score'] += edge['score']
                                    new_active_sequence['current_node_id'] = next_node_id
                                    new_active_sequences.append(new_active_sequence)
                                    if verbose == 2: print('    ++',new_active_sequence)

                            followed_peaks[next_node_id] = 1

                    active_sequences = new_active_sequences
                    if verbose: print(f"Number of active sequences: {len(active_sequences)}, finished sequences: {len(finished_sequences)}")
                    #for active_sequence in active_sequences:
                    #    print(active_sequence)
                    #print('**',active_sequences)

                    if len(active_sequences) == 0:
                        done = True



        #### Sort all results
        sorted_finished_sequences = sorted(finished_sequences,key=lambda x: x['score'], reverse=True)



        #### Loop over and report cterm and nterm results separately
        cterm_sequences = []
        nterm_sequences = []
        complete_sequences = []
        for sequence in sorted_finished_sequences:
            if len(sequence['sequence']) <= 2:
                continue

            if sequence['sequence'][0].endswith('-') and sequence['sequence'][-1].endswith('-'):
                continue
            if sequence['sequence'][0].startswith('-') and sequence['sequence'][-1].startswith('-'):
                continue

            if sequence['sequence'][0].startswith('-'):
                sequence['sequence'].reverse()
            if sequence['sequence'][-1].endswith('-'):
                sequence['sequence'].reverse()

            if sequence['sequence'][0].endswith('-') and sequence['sequence'][-1].startswith('-'):
                complete_sequences.append(sequence)
            elif sequence['sequence'][-1].startswith('-'):
                cterm_sequences.append(sequence)
            elif sequence['sequence'][0].endswith('-'):
                nterm_sequences.append(sequence)

        print("#### Best complete sequences:")
        counter = 0
        for sequence in complete_sequences:
            score = ''
            if counter < 2:
                self.annotate(spectrum, peptidoform_string=''.join(sequence['sequence']), charge=spectrum.analytes['1']['charge state'])
                #score = spectrum.attributes['psm_score']['diagnostic']['percent'] * spectrum.attributes['psm_score']['diagnostic']['count']
                score = spectrum.attributes['psm_score']['diagnostic']['percent']
            print(f"{sequence['score']}\t{sequence['sequence']}\t{''.join(sequence['sequence'])}\t{score}")
            counter += 1
            if counter > 50:
                break

        print("#### Best cterm sequences:")
        counter = 0
        for sequence in cterm_sequences:
            print(f"{sequence['score']}\t{sequence['sequence']}\t{''.join(sequence['sequence'])}")
            counter += 1
            if counter > 50:
                break

        print("#### Best nterm sequences:")
        counter = 0
        for sequence in nterm_sequences:
            print(f"{sequence['score']}\t{sequence['sequence']}\t{''.join(sequence['sequence'])}")
            counter += 1
            if counter > 50:
                break














    ####################################################################################################
    # Put all the peaks into a dict by integer mass to make lookups faster
    def reindex_network(self,spectrum):

        # Create an empty index
        spectrum.network['mz_index'] = {}

        # Loop over all nodes, putting them into an integer bin
        for node_id,node in spectrum.network['nodes'].items():
            int_node_mz = int(node['mz'])
            if int_node_mz not in spectrum.network['mz_index']:
                spectrum.network['mz_index'][int_node_mz] = []
            spectrum.network['mz_index'][int_node_mz].append(node_id)


    ####################################################################################################
    # Plot the spectrum network
    def plot_network(self,spectrum):

        import matplotlib.pyplot as plt

        x = []
        y = []

        # Loop over all nodes, creating x and y
        for node_id,node in spectrum.network['nodes'].items():
            x.append(node['mz'])
            y.append(node['mass_defect'])

        plt.scatter(x, y, color='green')
        plt.title(f"Spectrum mass defect plot")
        plt.xlabel('mz')
        plt.ylabel('mass defect')
        plt.grid(True)

        for node_id,node in spectrum.network['nodes'].items():
            plt.text(node['mz']+5, node['mass_defect'] - 0.02, str(node_id), size=8, color='black')

        plt.show()


    ####################################################################################################
    # Find the nodes in the network within the tolerance
    def find_close_nodes(self,spectrum,search_mz,tolerance):

        #### FIXME. In order to speed this up, ensure that all the peaks in each bin are sorted in order
        #### and then stop looking once we pass the +tolerance
        #### FIXME. It would probably be best to sort the matches in order of bestness

        matches = []

        # Get the integer mass as a dict key
        int_search_mz = int(search_mz)
        if int_search_mz not in spectrum.network['mz_index']:
            return matches

        # Loop over all the peaks in this bin and add them to matches if they're within tolerance
        for node_id in spectrum.network['mz_index'][int_search_mz]:
            mz = spectrum.network['nodes'][node_id]['mz']
            delta = search_mz - mz
            delta_ppm = delta / search_mz * 1e6

            if delta_ppm < -1 * tolerance or delta_ppm > tolerance:
                continue

            # Compute a delta score based on the delta and the tolerance
            delta_score = 1.0
            commonness_score = 1.0

            intensity = spectrum.network['nodes'][node_id]['intensity'] / spectrum.attributes['maximum intensity']
            intensity_score = intensity
            if intensity >= 0.05:
                intensity_score = intensity * 15.0 - 0.70
            if intensity >= 0.10:
                intensity_score = intensity * 0.36 + 0.62

            match = {
                'mz': mz,
                'node_id': node_id,
                'delta': delta,
                'delta_ppm': delta_ppm,
                'delta_score': delta_score,
                'commonness_score': commonness_score,
                'score': commonness_score * delta_score * intensity_score,
                'interpretation_string': '???',
            }

            matches.append(match)

        return matches


    ####################################################################################################
    # Starting from the specified peak, de novo in the specified direction
    def link_nodes(self, spectrum, sequencing_parameters, start_i_peak, direction):

        direction_multiplier = 1
        opposite_direction = 'downward'
        if direction == 'downward':
            direction_multiplier = -1
            opposite_direction = 'upward'

        #### Extract some sequence parameters
        precursor_mz = sequencing_parameters['precursor_mz']
        precursor_charge = sequencing_parameters['precursor_charge']
        tolerance = sequencing_parameters['tolerance']
        precursor_mw = sequencing_parameters['precursor_mw']
        potential_termini = sequencing_parameters['potential_termini']

        verbose = 0
        if start_i_peak == 36:
            verbose = 2

        # Get handles for some needed reference masses
        masses = self.mass_reference.atomic_masses
        residue_masses = self.mass_reference.nr_aa_masses
        ion_series_attr = self.mass_reference.ion_series_attributes

        start_i_peak = str(start_i_peak)

        # Create a starting sequence record
        mz = spectrum.network['nodes'][start_i_peak]['mz']
        sequence = {
            'terminus_type': 'unknown',
            'start_i_peak': start_i_peak,
            'start_peak_mz': mz,
            'direction': direction_multiplier,
            'working_mass': mz,
            'complement_mass': precursor_mw - ( mz - masses['proton'] ),
            'residues': [],
            'preliminary_score': spectrum.network['nodes'][start_i_peak]['intensity'] / spectrum.attributes['maximum intensity'],
        }
        if direction == 'downward':
            sequence['complement_mass'] = mz - masses['proton']


        #### Set up real or stub terminus information to test
        if start_i_peak == 'top' or start_i_peak == 'bottom':
            potential_termini = sequencing_parameters['potential_termini']
        else:
            potential_termini = [ {
                'ion_series': 'any',
                'terminus_type': None,
                'terminus_name': '',
                'terminus_mass': 0.0 } ]

        #### Loop over all potential termini (or just a stub non-terminus if not a terminal edge)
        for potential_terminus in potential_termini:

            terminus_mass = potential_terminus['terminus_mass']
            if potential_terminus['terminus_type'] is not None:
                terminus_mass += masses['proton']
            #print(f"Considering potential_terminus name={potential_terminus['terminus_name']} mass={potential_terminus['terminus_mass']}")

            #### Loop over single hops first, then double hops
            hops = [ 'single_residues', 'double_residues' ]
            successful_single_hops = {}
            for hop in hops:

                #### Set the residue mass options depending on hop
                if hop == 'single_residues':
                    if verbose: print(f"Looking to extend {sequence['working_mass']} {direction} with single residues")
                    aa_masses = residue_masses
                else:
                    if verbose: print(f"Looking to extend {sequence['working_mass']} {direction} with residue pairs")
                    aa_masses = self.mass_reference.nr_double_aa_masses

                #### Loop over all the potential residue masses
                for residue,residue_mass in aa_masses.items():

                    #### Skip I because it is the same as L
                    if residue == 'I':
                        continue
                    #### Skip a pair if one of the residues in the pair matches a successful single hop
                    if hop == 'double_residues':
                        pair = self.mass_reference.nr_double_aa_masses[residue]['residues']
                        skip_residue = False
                        for item in pair:
                            if item in successful_single_hops:
                                #print(f"Skipping {residue} because there is an {item} already in successful_single_hops")
                                skip_residue = True
                        if skip_residue:
                            continue
                        residue_mass = self.mass_reference.nr_double_aa_masses[residue]['mass']

                    lookahead_mz = sequence['working_mass'] + ( residue_mass + terminus_mass ) * direction_multiplier
                    if verbose: print(f"  Look for {residue} at {lookahead_mz}...")
                    matches = self.find_close_nodes(spectrum,lookahead_mz,tolerance)
                    if verbose: print(matches)
                    if len(matches) > 1:
                        print("WARNING: There are multiple matches!")
                        print(matches)
                        matches = [ matches[0] ]
                    for match in matches:
                        i_match_node = match['node_id']
                        if verbose: print(f"  Found match for {residue} at {match['mz']} with delta_ppm {match['delta_ppm']} ppm and score {match['score']}")

                        #### Since we always find edges from the top and bottom nodes specially with potential termini, if we accidentally land on one, ignore it
                        if i_match_node == 'top' or i_match_node == 'bottom':
                            continue

                        node1_id = str(start_i_peak)
                        node2_id = str(i_match_node)
                        node1_direction = direction
                        node2_direction = opposite_direction
                        if node1_id > node2_id:
                            node2_id = str(start_i_peak)
                            node1_id = str(i_match_node)
                            node2_direction = direction
                            node1_direction = opposite_direction

                        edge_id = f"{node1_id}-{node2_id}({potential_terminus['terminus_name']}{residue})"
                        edge = {
                            'edge_id': edge_id,
                            'node1_id': node1_id,
                            'node2_id': node2_id,
                            'direction': direction,
                            'score': match['score'],
                            'delta_ppm': match['delta_ppm'],
                            'type': 'residue',
                            'label': residue,
                        }
                        #### If this is a terminal edge, add that information
                        if potential_terminus['terminus_type'] is not None:
                            edge['terminus_type'] = potential_terminus['terminus_type']
                            edge['terminus_name'] = potential_terminus['terminus_name']
                            edge['score'] += 1.0

                        if edge_id not in spectrum.network['edges']:
                            spectrum.network['edges'][edge_id] = edge

                        spectrum.network['nodes'][node1_id][node1_direction+'_edges'][edge_id] = match['score']
                        spectrum.network['nodes'][node2_id][node2_direction+'_edges'][edge_id] = match['score']

                        #### If this is a single hop, record that it was a successful so we don't try it with a double hop
                        if hop == 'single_residues':
                            successful_single_hops[residue] = 1










    ####################################################################################################
    def analyze_residuals(self, spectrum):

        show_interactive_plots = 0

        for residual_type in [ 'relative','absolute' ]:
            #print(f"Analyzing {residual_type} residuals")
            residuals = self.residuals[residual_type]['ppm_deltas']
            #print(residuals)
            n_residuals = len(residuals)
            sorted_residuals = sorted(residuals)
            median = sorted_residuals[int(n_residuals/2)]
            q1 = sorted_residuals[int(n_residuals * 0.25 )]
            q3 = sorted_residuals[int(n_residuals * 0.75 )]
            siqr = ( q3 - q1 ) / 2.0
            #print(f"  n={n_residuals}, median={median}, q1={q1}, q3={q3}, siqr={siqr}")
            if show_interactive_plots:
                import matplotlib.pyplot as plt
                x = range(n_residuals)
                plt.scatter(x,residuals)
                plt.plot([0,n_residuals],[median,median])
                plt.plot([0,n_residuals],[q1,q1])
                plt.plot([0,n_residuals],[q3,q3])
                plt.show()

        spectrum.attributes['mass_accuracy'] = {
            'offset': median,
            'siqr': siqr,
            'is_optimized': True,
            'best_tolerance': siqr,
            'middle_tolerance': 2 * siqr,
            'outer_tolerance': 5 * siqr,
            'max_tolerance': 10.0,
        }

        if spectrum.attributes['mass_accuracy']['outer_tolerance'] > spectrum.attributes['mass_accuracy']['max_tolerance']:
            spectrum.attributes['mass_accuracy']['max_tolerance'] = spectrum.attributes['mass_accuracy']['outer_tolerance'] + 5


        if show_interactive_plots:
            import matplotlib.pyplot as plt
            import numpy as np
            best_tolerance = spectrum.attributes['mass_accuracy']['best_tolerance']
            max_tolerance = spectrum.attributes['mass_accuracy']['max_tolerance']
            outer_tolerance = spectrum.attributes['mass_accuracy']['outer_tolerance']
            x_curve = np.arange(100) / 100 * ( 12 - 0) + 0
            y_curve = np.arange(100) * 0.0
            i = 0
            c = ( outer_tolerance - 0.1 * best_tolerance ) / ( outer_tolerance - best_tolerance )
            for x in x_curve:
                if x < best_tolerance:
                    y_curve[i] = 1
                else:
                    y_curve[i] = -0.9/(outer_tolerance-best_tolerance) * x + c
                if x > outer_tolerance:
                    y_curve[i] = -0.1/(max_tolerance-outer_tolerance) * x + 0.1 * max_tolerance / (max_tolerance-outer_tolerance)
                if y_curve[i] < 0:
                    y_curve[i] = 0
                i += 1
            plt.scatter(x_curve,y_curve)
            plt.plot([10,10],[0,1])
            plt.show()



    ####################################################################################################
    def rescore_interpretations(self, spectrum):

        #### If the spectrum mass accuracy information has not been optimized, then nothing to do
        mass_accuracy = spectrum.attributes['mass_accuracy']
        if mass_accuracy['is_optimized'] is False:
            return

        best_tolerance = mass_accuracy['best_tolerance']
        outer_tolerance = mass_accuracy['outer_tolerance']
        max_tolerance = mass_accuracy['max_tolerance']

        total_ion_current = 0.0
        categories = [ 'contamination', 'nondiagnostic', 'diagnostic', 'unexplained', 'unknown' ]
        metrics = [ 'intensity', 'count', 'percent' ]
        psm_score = {}
        for category in categories:
            psm_score[category] = {}
            for metric in metrics:
                psm_score[category][metric] = 0.0


        # Loop over all peaks shifting and rescoring the peak interpretations
        for i_peak in range(spectrum.attributes['number of peaks']):

            peak = spectrum.peak_list[i_peak]
            intensity = peak[PL_INTENSITY]
            best_score = 0

            #### Loop over the interpretations
            for interpretation in peak[PL_INTERPRETATIONS]:

                #### Unless the peak is a foreign precursor, correct delta for the previously-computed offset
                match = re.match(r'\d+\@p',interpretation[INT_INTERPRETATION_STRING])
                if match:
                    pass
                else:
                    interpretation[INT_DELTA_PPM] -= mass_accuracy['offset']

                #### Compute the absolute value of the delta ppm to use for the delta score
                abs_delta_ppm = abs(interpretation[INT_DELTA_PPM])

                # Compute a delta score
                x = abs_delta_ppm
                c = ( outer_tolerance - 0.1 * best_tolerance ) / ( outer_tolerance - best_tolerance )
                delta_score = 0.0
                if x < best_tolerance:
                    delta_score = 1.0
                else:
                    delta_score = -0.9 / (outer_tolerance - best_tolerance) * x + c
                    #print(f"**{delta_score}")
                if x > outer_tolerance:
                    delta_score = -0.1/(max_tolerance-outer_tolerance) * x + 0.1 * max_tolerance / (max_tolerance-outer_tolerance)
                    #print(f"##{delta_score}")
                if delta_score < 0.0:
                    delta_score = 0.0

                #if delta_score > 1:
                #    print(f"Yipe! {delta_score}, {x}")
                #    print(f"{best_tolerance}, {outer_tolerance}, {max_tolerance}")
                #    exit()


                interpretation[INT_DELTA_SCORE] = delta_score
                interpretation[INT_SCORE] = interpretation[INT_DELTA_SCORE] * interpretation[INT_COMMONNESS_SCORE]

                if interpretation[INT_SCORE] > best_score:
                    peak[PL_INTERPRETATION_STRING] = interpretation[INT_INTERPRETATION_STRING] + '/' + '{:.1f}'.format(interpretation[INT_DELTA_PPM]) + 'ppm'
                    best_score = interpretation[INT_SCORE]
                    peak[PL_ATTRIBUTES][PLA_DIAGNOSTIC_CATEGORY] = interpretation[INT_DIAGNOSTIC_CATEGORY]

            # Resolve isotopes
            if peak[PL_INTERPRETATION_STRING].startswith('isotope'):
                parent_peak = spectrum.peak_list[ peak[PL_ATTRIBUTES][PLA_PARENT_PEAK] ]

                # Inherit the category from the parent
                diagnostic_category = parent_peak[PL_ATTRIBUTES][PLA_DIAGNOSTIC_CATEGORY]
                peak[PL_ATTRIBUTES][PLA_DIAGNOSTIC_CATEGORY] = diagnostic_category

                parent_peak_interpretation = parent_peak[PL_INTERPRETATION_STRING]
                if peak[PL_ATTRIBUTES][PLA_IS_ISOTOPE] == 1:
                    isotope_string = '+i'
                else:
                    isotope_string = '+' + str(peak[PL_ATTRIBUTES][PLA_IS_ISOTOPE]) + 'i'
                if parent_peak_interpretation[0] == '?':
                    peak[PL_INTERPRETATION_STRING] = f"? # ISO {parent_peak[PL_MZ]}{isotope_string}/" + '{:.1f}'.format(interpretation[INT_DELTA_PPM]) + 'ppm'
                else:
                    # Strip off the delta after the slash
                    isotope_interpretation = re.sub(r'/.+','',parent_peak_interpretation)
                    # See if there's a charge string
                    match = re.search(r'(\^\d+)',isotope_interpretation)
                    charge_string = ''
                    if match:
                        charge_string = match.group(1)
                        isotope_interpretation = re.sub(r'\^\d+','',isotope_interpretation)
                    else:
                        if peak[PL_ATTRIBUTES][PLA_CHARGE] > 1:
                            charge_string = f"^{peak[PL_ATTRIBUTES][PLA_CHARGE]}"
                    isotope_interpretation += f"{isotope_string}{charge_string}/" + '{:.1f}'.format(interpretation[INT_DELTA_PPM]) + 'ppm'
                    peak[PL_INTERPRETATION_STRING] = isotope_interpretation

            else:
                diagnostic_category = peak[PL_ATTRIBUTES][PLA_DIAGNOSTIC_CATEGORY]

            # Record the intensity under the appropriate bin
            psm_score[diagnostic_category]['intensity'] += intensity
            psm_score[diagnostic_category]['count'] += 1
            total_ion_current += intensity

        for key in psm_score:
            psm_score[key]['percent'] = psm_score[key]['intensity'] / total_ion_current
        psm_score['total_ion_current'] = total_ion_current
        spectrum.attributes['psm_score'] = psm_score


    ####################################################################################################
    #### Return a printable buffer string of the details of the peptidoform
    def show(self):

        buf = ''
        buf += f"Peptidoform_string={self.peptidoform.peptidoform_string}\n"
        buf += f"Charge={self.spectrum_attributes['charge']}\n"
        for i_peak in range(len(self.predicted_fragments_list)):
            interpretations_string = ''
            for interpretation in self.predicted_fragments_list[i_peak][1]:
                if interpretations_string == '':
                    interpretations_string = interpretation[0]
                else:
                    interpretations_string += ', ' + interpretation[0]
            buf += '  ' + '{:10.4f}'.format(self.predicted_fragments_list[i_peak][0]) + '  ' + interpretations_string + "\n"
        return buf


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Class representing a peptidoform')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('--usi', action='store', default=None, help='USI to process')
    argparser.add_argument('--annotate', action='count', help='If set, annotate the USI spectrum' )
    argparser.add_argument('--denovo', action='count', help='If set, interpret the USI spectrum de novo' )
    argparser.add_argument('--examine', action='count', help='If set, examine the spectrum to see what can be learned' )
    argparser.add_argument('peptidoform_string', type=str, nargs='*', help='Optional peptidoform strings to parse')
    params = argparser.parse_args()

    # Set verbose mode
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    # If there are supplied peptidoform strings, use those, else set a default example
    peptidoform_strings_list = params.peptidoform_string
    if len(peptidoform_strings_list) == 0:
        peptidoform_strings_list = [ '[Acetyl]-GS[Phospho]PVLPHEPAK' ]
        peptidoform_strings_list = [ 'GQEY[Phospho]LLLEK' ]
        peptidoform_strings_list = [ ]

    #### Loop over all the peptidoform strings and decompose them
    for peptidoform_string in peptidoform_strings_list:
        print('**************************************************************************')
        peptidoform = Peptidoform(peptidoform_string, verbose=verbose)
        print(peptidoform.show())
        annotator = SpectrumAnnotator()
        annotator.predict_fragment_ions(peptidoform, charge=1)
        print(annotator.show())
    #return


    #### Get the user-supplied USI
    usi_string = ''
    if params.usi:
        usi_string = params.usi

    #### Or, use an example one
    if not usi_string:
        #usi_string = 'mzspec:PXD007058:SF_200217_pPeptideLibrary_pool1_HCDOT_rep1:scan:4336:GQEY[Phospho]LILEK/2'
        usi_string = 'mzspec:PXD005336:Dinaciclib_00952_B05_P008348_B00_A00_R1:scan:16419:LLSILSR/2'
        usi_string = 'mzspec:PXD000865::0223_F1_R1_P73514B20_TMT8:scan:19806:[TMT6plex]-LSEQEGSR/2' # nice annotation. Lots of immonium ions. several probably could be annotated with NIST system
        #usi_string = 'mzspec:PXD002255:ES_XP_Ubi_97H_HCD_349:scan:9617:LAEIYVNSSFYK/2'
        #usi_string = 'mzspec:PXD000612:20120224_EXQ5_KiSh_SA_LabelFree_HeLa_Phospho_EGF_rep4_FT3:scan:17080:VAPLS[Phospho]PGK/2' # Very nice annotation!
        #usi_string = 'mzspec:PXD000612:20120413_EXQ5_KiSh_SA_LabelFree_HeLa_pY_pervandate_rep2:scan:24932:RGIPNY[Phospho]EFK/2' # Nice. But the isotopes are peculiar on this one!
        usi_string = 'mzspec:PXD015910:TCGA_AR-A0TV_C8-A12Z_AO-A0JJ_117C_P_BI_20130831_H-PM_f04:scan:6346:[iTRAQ4plex]-AQETS[Phospho]GEEISK[iTRAQ4plex]/2'
        usi_string = 'mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2'
        usi_string = 'mzspec:PXD000966:CPTAC_CompRef_00_iTRAQ_12_5Feb12_Cougar_11-10-11.mzML:scan:11850:[UNIMOD:214]-YYWGGLYSWDMSK[UNIMOD:214]/2'
        usi_string = 'mzspec:PXD000966:CPTAC_CompRef_00_iTRAQ_12_5Feb12_Cougar_11-10-11.mzML:scan:11850:[iTRAQ4plex]-YYWGGLYSWDMSK[iTRAQ4plex]/2'
        usi_string = 'mzspec:PXD002286:081213-Wittie-Phos-1A:scan:14366:MVC[Carbamidomethyl]S[Phospho]PVTVR/2'
        #usi_string = 'mzspec:MSV000085202:210320_SARS_CoV_2_T:scan:131256:DLPQGFSALEPLVDLPIGIN[HexNac]ITR/2'
        usi_string = 'mzspec:PXD010154:01283_A02_P013187_S00_N09_R1:scan:30190:ELVISYLPPGM[Oxidation]ASK/2'
        usi_string = 'mzspec:PXD900000:201102-HeLa-TMT-01:scan:100033:[TMT]-VLPSIVNEVLK[TMT]/2'
        usi_string = 'mzspec:PXD900000:201102-HeLa-TMT-01:scan:100053:[TMT]-NLPFDFTWK[TMT]/2'
        usi_string = 'mzspec:PXD000561:Adult_NKcells_Gel_Elite_78_f01:scan:9715:IPEISIQDM[UNIMOD:35]TAQVTSPSGK/2'

    #### Fetch the spectrum
    spectrum = Spectrum()
    spectrum.fetch_spectrum(usi_string)

    usi = UniversalSpectrumIdentifier(usi_string)
    print(json.dumps(usi.__dict__,sort_keys=True,indent=2))
    peptidoform_string = usi.peptidoform
    charge = int(usi.charge)

    annotator = SpectrumAnnotator()

    # Prepare to annotate the spectrum de novo
    if params.denovo:
        # First, clean the spectrum of stuff that is not useful for de novo sequencing
        annotator.identify_isotopes(spectrum)
        annotator.remove_isotopes(spectrum)
        annotator.identify_precursors(spectrum)
        annotator.remove_precursors(spectrum)
        annotator.identify_low_mass_ions(spectrum)
        annotator.remove_low_mass_ions(spectrum)
        # After cleaning, recompute the metrics on what's left and re-index
        annotator.compute_spectrum_metrics(spectrum)
        annotator.index_peaks(spectrum)

        if False:
            annotator.identify_neutral_losses(spectrum)
            annotator.annotate_peptidoform(spectrum, peptidoform_string=peptidoform_string, charge=charge)
            annotator.identify_precursors(spectrum)
            annotator.analyze_residuals(spectrum)
            annotator.rescore_interpretations(spectrum)
            print(spectrum.show(show_all_annotations=False))
            exit()

        # Set some basic input parameters
        sequencing_parameters = {
            'fragmentation_type': 'HCD',
            'tolerance': 30.0,
            'precursor_mz': spectrum.analytes['1']['precursor_mz'],
            'precursor_charge': spectrum.analytes['1']['charge state'],
        }
        # And go!
        annotator.interpret_de_novo(spectrum, sequencing_parameters=sequencing_parameters)


    # Annotate the spectrum
    if params.annotate:
        annotator.annotate(spectrum, peptidoform_string=peptidoform_string, charge=charge)
        print(spectrum.show(show_all_annotations=False))


    # Examine the spectrum
    if params.examine:
        #annotator.identify_isotopes(spectrum)
        #annotator.remove_isotopes(spectrum)
        annotator.compute_spectrum_metrics(spectrum)
        annotator.index_peaks(spectrum)
        annotator.identify_complement_ions(spectrum)

        #annotator.identify_low_mass_ions(spectrum)
        #annotator.identify_reporter_ions(spectrum)
        #annotator.identify_neutral_losses(spectrum)
        #annotator.identify_precursors(spectrum)
        #annotator.analyze_residuals(spectrum)
        #annotator.rescore_interpretations(spectrum)
        print(spectrum.show(show_all_annotations=False))

    #### Display the final result
    #print(spectrum.show(show_all_annotations=False))


#### For command line usage
if __name__ == "__main__": main()
