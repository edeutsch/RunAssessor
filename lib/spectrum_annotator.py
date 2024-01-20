#!/usr/bin/env python3
import sys
import os
import argparse
import os.path
import re
import itertools
import requests
import json
import copy
import pandas as pd
import numpy as np
import numpy
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
#### Calculate all non-redundant permutations of a list of potential neutral losses
def get_nr_permutations(input_list, max_of_each=3):

    #### Preparation
    all_combinations = set()

    #### Reduce the list of all possible losses to a maximum of N of each type
    reduced_loss_list = []
    loss_counts = {}
    for item in input_list:
        if item not in loss_counts:
            loss_counts[item] = 0
        loss_counts[item] += 1
        if loss_counts[item] <= max_of_each:
            reduced_loss_list.append(item)

    sorted_losses_list = sorted(reduced_loss_list)

    #### Generate permutations
    for set_size in range(len(sorted_losses_list) + 1):
        combinations_tuple_list = sorted(list(itertools.combinations(sorted_losses_list, set_size)))
        combinations_set_list = [ list(i) for i in combinations_tuple_list ]
        combinations_sorted_set_list = [ sorted(item) for item in combinations_set_list ]

        for combination_set in combinations_sorted_set_list:
            combination_set_str = ';'.join(combination_set)
            all_combinations.add(combination_set_str)

    all_combinations = [ item.split(';') for item in all_combinations ]
    return(all_combinations)





####################################################################################################
#### SpectrumAnnotator class
class SpectrumAnnotator:
    '''
    - annotate()                Annotate a spectrum by calling a series of methods given a peptidoform
    - predict_fragment_ions()   Predict all potential fragment ions for the provided peptidoform
    - annotate_peptidoform()    Annotate the spectrum with the predicted fragments from a supplied peptidoform
    - compute_spectrum_score()  !!! FIXME Not quite clear what is going on here
    - find_close_predicted_fragments()  For a given observed m/z, find all the potential matching predicted fragments within tolerance
    - index_peaks()             Put all the peaks into a dict keyed by integer mass to make lookups faster
    - find_close_ions()         Find the closest predicted fragment
    - add_interpretation()      Add an interpretation to a peak
    - identify_isotopes()       Loop through a spectrum and identify all the peaks that are isotopes of another
    - identify_precursors()     Identify all the peaks in a spectrum that might be precursor-related and annotate them
    - identify_neutral_losses() (NOT CURRENTLY USED?) Identify all peaks that might be a neutral loss of another
    - identify_low_mass_ions()  Identify all the low-mass ions in a spectrum that aren't specifically part of predicted fragments
    - identify_reporter_ions()  Identify known reporter ions, somewhat independently of whether they should be there or not
    - compute_spectrum_metrics()    Compute a set of metrics from a spectrum such as SNR
    - analyze_residuals()       Analyze and potentially plot a set of residuals of a spectrum
    - rescore_interpretations() Rescore all the potential interpretations of a peak to select a winner
    - show()                    Return a printable buffer string of the details of the peptidoform and the annotations of all peaks
    - plot()                    Plot the spectrum and its annotations in a nice publishable way
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

        # Set up a list for the predicted fragments
        self.predicted_fragments_list = []
        self.predicted_fragments_index = {}

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
    #### Annotate a spectrum by calling a series of methods given a peptidoform
    def annotate(self, spectrum, peptidoform, charge, tolerance=None):

        if tolerance is None:
            tolerance = self.tolerance
        else:
            self.tolerance = tolerance

        spectrum.compute_spectrum_metrics()
        self.identify_isotopes(spectrum)
        self.identify_low_mass_ions(spectrum, peptidoform)
        self.identify_reporter_ions(spectrum)
        #self.identify_neutral_losses(spectrum)
        self.annotate_peptidoform(spectrum, peptidoform=peptidoform, charge=charge)
        self.identify_precursors(spectrum)
        self.analyze_residuals(spectrum)
        self.rescore_interpretations(spectrum)


    ####################################################################################################
    #### Predict all potential fragment ions for the provided peptidoform
    def predict_fragment_ions(self, peptidoform=None, charge=1, fragmentation_type='HCD', skip_internal_fragments=False):

        #skip_internal_fragments = True
        if DEBUG:
            eprint("DEBUG: Entering predict_fragment_ions")

        # Use the passed peptidoform object or else the previously supplied one or return None
        if peptidoform is None:
            if self.peptidoform is None:
                return
            else:
                peptidoform = self.peptidoform
        else:
            self.peptidoform = peptidoform

        debug = False

        special_annotation_rules = {}

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
        residue_masses = self.mass_reference.aa_masses
        ion_series_attr = self.mass_reference.ion_series_attributes
        neutral_losses = self.mass_reference.neutral_losses
        neutral_losses_by_residue = self.mass_reference.neutral_losses_by_residue
        neutral_losses_by_formula = self.mass_reference.neutral_losses_by_formula
        ref_terminal_modifications = self.mass_reference.terminal_modifications

        # Determine the terminal modification masses
        have_labile_nterm_mod = False                                       # FIXME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        terminal_mass_modifications = { 'nterm': 0.0, 'cterm': 0.0 }
        if peptidoform.terminal_modifications is not None and 'nterm' in peptidoform.terminal_modifications:
            nterm_attrs = peptidoform.terminal_modifications['nterm']
            terminal_mass_modifications['nterm'] = nterm_attrs['delta_mass']
            if nterm_attrs['modification_name'].startswith('TMT'):
                special_annotation_rules[nterm_attrs['modification_name']] = True
            if debug:
                eprint(f"INFO: Add n-terminal mass modification {nterm_attrs['modification_name']} as {terminal_mass_modifications['nterm']}")

        if peptidoform.terminal_modifications is not None and 'cterm' in peptidoform.terminal_modifications:
            cterm_attrs = peptidoform.terminal_modifications['cterm']
            terminal_mass_modifications['cterm'] = cterm_attrs['delta_mass']
            if debug:
                eprint(f"INFO: Add c-terminal mass modification {cterm_attrs['modification_name']} as {terminal_mass_modifications['cterm']}")

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
        peptide_length = len(peptidoform.peptide_sequence)
        all_annotations = {}

        # Main loop: iterate through each position, working from both ends simultaneously
        for i_residue in range(peptide_length):

            if DEBUG:
                eprint(f"DEBUG: Processing i_residue={i_residue} for peptide_length={peptide_length}")

            # Add additional series entries for internal ions
            if not skip_internal_fragments:
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

            #### On the last pass for the whole peptide, just compute b and y bions, but they will be relabeled a p (precursor)
            #### Also, just for charge 1 and the precursor charge. This is a little hokey. it should probably broken out to a
            #### separate system that works independent of the peptidoform
            charge_list = range(1, charge + 1)
            if i_residue == peptide_length - 1 and have_labile_nterm_mod is False:
                series_list = [ 'b', 'y' ]
                series_list = [ 'b' ]
                charge_list = [ 1, charge ]


            # Generate fragments for each ion series we expect to find
            for series in series_list:

                series_type = series[0]
                cumulative_residues[series] += 1

                # And for each expected charge
                for i_charge in charge_list:

                    if DEBUG:
                        eprint(f"DEBUG:   - Processing {series} in {series_list}, {i_charge} in {charge_list}")

                    # Get the next residue
                    residue_offset = i_residue + 1
                    if ion_series_attr[series_type]['terminus_type'] == 'cterm':
                        residue_offset = peptide_length - i_residue
                    residue = peptidoform.residues[residue_offset]['residue_string']
                    base_residue = residue
                    if len(residue) > 1:
                        base_residue = peptidoform.residues[residue_offset]['base_residue']
                    residue_mass = residue_masses[base_residue]
                    if len(residue) > 1:
                        residue_mass += peptidoform.residues[residue_offset]['delta_mass']

                    # Only compute certain things on the first pass
                    if i_charge == 1:
                        # Update the cumulative mass
                        cumulative_mass[series] += residue_mass

                        #### If this is the precursor pass, also add in the other terminus
                        if i_residue + 1 == peptide_length and series == 'b':
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
                    all_combinations = get_nr_permutations(losses_list, max_of_each=2)

                    # Create the annotations for each combination of losses (including no loss)
                    if DEBUG:
                        #eprint(f"DEBUG:       - Processing {len(all_combinations)} combinations")
                        pass
                    for potential_neutral_loss_combination_list in all_combinations:
                        loss_string = ''
                        loss_mass = 0.0
                        for potential_neutral_loss_formula in potential_neutral_loss_combination_list:
                            if potential_neutral_loss_formula != '':
                                loss_string += f"-{potential_neutral_loss_formula}"
                                loss_mass += neutral_losses_by_formula[potential_neutral_loss_formula]['delta_mass']

                        # Create the default interpretation
                        interpretation = f"{series}{i_residue + 1}{loss_string}"

                        #### If this is the final pass for the precursor
                        if i_residue + 1 == peptide_length:
                            interpretation = f"p{loss_string}"

                        #### If there is an n-terminal mod, then add a precursor with a loss of the n-terminal mod during the y series
                        #print(f"{i_residue}, {peptide_length - 1}. {series}, {list(peptidoform.terminal_modifications.keys())} -- {special_annotation_rules}")
                        #if i_residue + 1 == peptide_length and series == 'b' and i_charge == 1 and peptidoform.terminal_modifications is not None and 'nterm' in peptidoform.terminal_modifications:
                        if i_residue + 1 == peptide_length and i_charge == 1 and peptidoform.terminal_modifications is not None and 'nterm' in peptidoform.terminal_modifications:
                            #print(losses_list)
                            #print(all_combinations)
                            #exit()

                            #### If there are special annotation rules, apply those
                            apply_special_rules = None
                            possible_special_rules = [ 'TMT6plex', 'TMTpro' ]
                            for possible_special_rule in possible_special_rules:
                                if possible_special_rule in special_annotation_rules:
                                    apply_special_rules = possible_special_rule
                            if apply_special_rules and i_charge == 1:
                                for special_ion_name, special_ion_data in self.mass_reference.special_label_losses[apply_special_rules].items():
                                    special_interpretation = f"p-{special_ion_name}{loss_string}"
                                    special_interpretation_score = 55
                                    mz = cumulative_mass[series] - special_ion_data['mz'] - loss_mass + masses['proton'] * i_charge
                                    self.predicted_fragments_list.append( [ mz, [ [ special_interpretation, special_interpretation_score ] ] ] )

                            mz = cumulative_mass[series] - terminal_mass_modifications['nterm'] - loss_mass + masses['proton'] * i_charge
                            special_interpretation = f"p-[{peptidoform.terminal_modifications['nterm']['modification_name']}]{loss_string}"
                            special_interpretation_score = 56
                            self.predicted_fragments_list.append( [ mz, [ [ special_interpretation, interpretation_score ] ] ] )

                        # But if this is an internal fragment, create that style
                        if series_type == 'm':
                            if cumulative_residues[series] > 1:
                                interpretation = f"{series}:{i_residue + 1}{loss_string}"
                            # Unless it is only of length 1, then skip this entirely because that is immonium, and immonium ions are handled statically
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
                            #print(mz,interpretation)


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
    #### Annotate the spectrum with the predicted fragments from a supplied proforma peptidoform string
    def annotate_peptidoform(self, spectrum, peptidoform, charge, skip_internal_fragments=False, tolerance=None):

        #### If peptidoform is specified as None, then there's nothing we can do here
        if peptidoform is None:
            return

        if tolerance is None:
            tolerance = self.tolerance
        else:
            self.tolerance = tolerance

        if not isinstance(peptidoform, ProformaPeptidoform):
            eprint(f"ERROR: Passed peptidoform is not an object, but instead is type {type(peptidoform)}")
            return

        stripped_sequence = peptidoform.peptide_sequence
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
            matches = self.find_close_predicted_fragments(mz, tolerance)
            if matches:
                diagnostic_category = 'diagnostic'
                for match in matches:
                    #peak[PL_INTERPRETATION_STRING] = f"{match[INT_INTERPRETATION_STRING]}/" + '{:.1f}'.format(match[INT_DELTA_PPM]) + 'ppm'
                    #peak[PL_INTERPRETATIONS].append(match)
                    if match[INT_INTERPRETATION_STRING].startswith('p'):
                        peak[PL_ATTRIBUTES][PLA_IS_PRECURSOR] += 1
                        diagnostic_category = 'nondiagnostic'

                    self.add_interpretation(peak, match, diagnostic_category=diagnostic_category, residual_type='absolute')


    ####################################################################################################
    #### FIXME Not quite clear what is going on here
    def compute_spectrum_score(self, spectrum, peptidoform, charge):

        tolerance = self.tolerance

        # Store the stripped sequence if present
        stripped_sequence = ''
        if peptidoform.stripped_sequence is not None:
            stripped_sequence = peptidoform.stripped_sequence

        self.predict_fragment_ions(peptidoform=peptidoform, charge=charge, fragmentation_type='HCD', skip_internal_fragments=True)

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
            matches = self.find_close_predicted_fragments(mz, tolerance)
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
    # For a given observed m/z, find all the potential matching predicted fragments within tolerance
    def find_close_predicted_fragments(self,observed_mz, tolerance):

        # Override the input tolerance with the reference tolerance
        #tolerance = self.stats['mz_calibration_tolerance_ppm']

        # We will return a list of possible matches
        matches = []

        # Get the integer mass as a dict key
        int_observed_mz = int(observed_mz)
        if int_observed_mz not in self.predicted_fragments_index:
            return matches

        # Loop over all the peaks in this bin and add them to matches if they're within tolerance
        for predicted_fragment in self.predicted_fragments_index[int_observed_mz]:
            fragment_mz = predicted_fragment[0]
            #delta = fragment_mz - observed_mz                          #### REVERSEDHERE
            delta = observed_mz - fragment_mz
            delta_ppm = delta / observed_mz * 1e6
            #if delta_ppm < -1 * tolerance:                          #### REVERSEDHERE
            if delta_ppm > tolerance:
                continue
            #if delta_ppm > tolerance:                          #### REVERSEDHERE
            if delta_ppm < -1 * tolerance:
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

        # We will return a list of possible matches
        matches = []

        # Get the integer mass as a dict key
        int_search_mz = int(search_mz)
        if int_search_mz not in spectrum.peak_index:
            return matches

        # Loop over all the peaks in this bin and add them to matches if they're within tolerance
        for peak in spectrum.peak_index[int_search_mz]:
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

        # Constants
        average_isotope_delta = 1.003355    # This is the official mass delta of carbon 13 over caerbon 12 and seems to work best
        max_charge = 4
        debug = False

        # Define some basic parameters
        n_peaks = spectrum.attributes['number of peaks']
        #tolerance = self.stats['mz_calibration_tolerance_ppm']
        tolerance = self.tolerance

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

            #i_peak += 1


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
                        match[INT_COMMONNESS_SCORE] = 95
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
                match[INT_DELTA_PPM] = -1 * match[INT_DELTA_PPM]                                                        #### Weird things because I reversed the sign of the delta?????

                match[INT_INTERPRETATION_STRING] = f"r[{reporter_ion_name}]"
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
    #### Analyze and potentially plot a set of residuals of a spectrum
    def analyze_residuals(self, spectrum):

        show_interactive_plots = 0
        show_debugging_info = False

        for residual_type in [ 'relative','absolute' ]:
            residuals = self.residuals[residual_type]['ppm_deltas']
            if show_debugging_info:
                print(f"Analyzing {residual_type} residuals")
                print(residuals)

            #### If there are not enough residuals, then just return without further ado
            n_residuals = len(residuals)
            if n_residuals < 3:
                if show_debugging_info:
                    print(f"Not enough residuals with n_residuals={n_residuals}")
                return

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

            if show_interactive_plots:
                min = -20.0
                max = 20.0
                binsize = 1
                n_bins = int( (max - min) / binsize )
                y_values, x_floor, patches = plt.hist( residuals, n_bins, [min,max])
                plt.xlabel('ppm delta')
                plt.ylabel('N')
                plt.title(f"Residuals for {residual_type} annotations")
                plt.xlim(min, max)
                plt.grid(True)
                print('****', x_floor)
                popt,pcov = curve_fit(gaussian_function, x_floor[0:-1], y_values, p0=[100.0, 0.0 ,binsize])
                plt.plot(x_floor + 0.5, gaussian_function(x_floor,*popt), 'r:')
                print(popt)
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
    #### Rescore all the potential interpretations of a peak to select a winner
    def rescore_interpretations(self, spectrum):

        #### If the spectrum mass accuracy information has not been optimized, then nothing to do
        mass_accuracy = spectrum.attributes['mass_accuracy']
        #if mass_accuracy['is_optimized'] is False:
        #    return

        best_tolerance = mass_accuracy['best_tolerance']
        outer_tolerance = mass_accuracy['outer_tolerance']
        max_tolerance = mass_accuracy['max_tolerance']

        total_ion_current = 0.0
        categories = [ 'contamination', 'nondiagnostic', 'diagnostic', 'unexplained' ]
        metrics = [ 'intensity', 'count', 'fraction' ]
        psm_score = {}
        for category in categories:
            psm_score[category] = {}
            for metric in metrics:
                psm_score[category][metric] = 0.0


        # Loop over all peaks shifting and rescoring the peak interpretations
        for i_peak in range(spectrum.attributes['number of peaks']):

            peak = spectrum.peak_list[i_peak]
            intensity = peak[PL_INTENSITY]
            best_score = 0.0
            diagnostic_category = 'unexplained'

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

            #### If the best score is 0, then there's no annotation
            if best_score == 0.0:
                peak[PL_INTERPRETATION_STRING] = '?'
                peak[PL_ATTRIBUTES][PLA_DIAGNOSTIC_CATEGORY] = 'unexplained'


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
                    parent_peak[PL_INTERPRETATION_STRING] = '?' + str(parent_peak[PL_I_PEAK])
                    peak[PL_INTERPRETATION_STRING] = '?' + str(parent_peak[PL_I_PEAK]) + isotope_string + '/' + '{:.1f}'.format(interpretation[INT_DELTA_PPM]) + 'ppm'
                    #peak[PL_INTERPRETATION_STRING] = f"? # ISO {parent_peak[PL_MZ]}{isotope_string}/" + '{:.1f}'.format(interpretation[INT_DELTA_PPM]) + 'ppm'
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
            psm_score[key]['fraction'] = psm_score[key]['intensity'] / total_ion_current
        psm_score['total_ion_current'] = total_ion_current
        spectrum.attributes['psm_score'] = psm_score


    ####################################################################################################
    #### Return a printable buffer string of the details of the peptidoform and the annotations of all peaks
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
    #### Plot the spectrum and its annotations in a nice publishable way
    def plot(self, spectrum, peptidoform_string, charge):
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        from matplotlib.backends.backend_pdf import PdfPages

        #pdf = PdfPages(f"spectrum.pdf")
        fig = plt.figure(figsize=(8,8))
        gridspec1 = gridspec.GridSpec(1, 1)
        plot1 = fig.add_subplot(gridspec1[0])
        gridspec1.tight_layout(fig, rect=[0.02, 0.3, 1, 1])

        #### Compute the min and max for mz and intensities
        max_intensity = 0
        min_mz = 9999
        max_mz = 0
        for peak in spectrum.peak_list:
            mz = peak[PL_MZ]
            intensity = peak[PL_INTENSITY]
            if intensity > max_intensity:
                max_intensity = intensity
            if mz < min_mz:
                min_mz = mz
            if mz > max_mz:
                max_mz = mz
        xmin = int(min_mz / 30) * 30
        xmax = int(max_mz / 30) * 30 + 30
        xmin = int(min_mz) - 15
        xmax = int(max_mz) + 15

        #### Set up the main spectrum plot
        #plot1.plot( [0,1000], [0,1], color='tab:green' )
        plot1.set_xlabel('m/z', fontname='Arial')
        plot1.set_ylabel('Relative Intensity', fontname='Arial')
        plot1.set_xlim([xmin, xmax])
        plot1.set_ylim([0,1])
        #ax[0,1].plot( [0,limit], [0,limit], '--', linewidth=1, color='gray')

        #### Set up the residuals plot
        gridspec2 = gridspec.GridSpec(1, 1)
        plot2 = fig.add_subplot(gridspec2[0])
        gridspec2.tight_layout(fig, rect=[0.02, 0.12, 1, 0.3])

        #plot2.plot( [0,1000], [-10,10], color='tab:green' )
        #plot2.set_xlabel('m/z')
        plot2.set_ylabel('delta (PPM)', fontname='Arial')
        plot2.set_xlim([xmin, xmax])
        plot2.set_xticklabels([])
        plot2.set_ylim([-16,16])
        plot2.plot( [0,xmax], [0,0], '--', linewidth=0.6, color='gray')

        #### Set up colors for different types of ion and a grid to track where items have been layed out
        colors = { 'b': 'tab:blue', 'a': 'tab:green', 'y': 'tab:red', '0': 'tab:gray', 'I': 'tab:orange', '?': 'tab:gray', 'p': 'tab:pink', 'm': 'tab:brown', 'r': 'tab:purple' }
        blocked = np.zeros((xmax,100))

        #### Write the peptide sequence to the plot in order to annotate it later
        stripped_peptide = re.sub(r'\[.+?\]','',peptidoform_string)
        stripped_peptide = stripped_peptide.replace('-','')
        residues = list(stripped_peptide)
        sequence_offset = 300
        sequence_gap = (xmax-xmin)/45
        sequence_height = 0.8
        counter = 0
        for residue in residues:
            plot1.text(sequence_offset+counter*sequence_gap, sequence_height, residue, fontsize='x-large', ha='center', va='bottom', color='black', fontname='Arial')
            counter += 1

        #### Loop over all peaks and plot them in the color that they would be annotated with
        counter = 0
        annotations = []
        saved_residuals = []
        for peak in spectrum.peak_list:
            i_peak = peak[PL_I_PEAK]
            mz = peak[PL_MZ]
            intensity = peak[PL_INTENSITY] / max_intensity * 0.93
            interpretations_string = peak[PL_INTERPRETATION_STRING]
            mz_delta = 99
            match = re.search(r'/([\+\-\d\.]+)ppm', interpretations_string)
            if match:
                mz_delta = float(match.group(1))
                mz_delta += spectrum.attributes['mass_accuracy']['offset']

            ion_type = interpretations_string[0]
            if ion_type in [ '1','2','3','4','5','6','7','8','9' ]:
                color = 'tab:olive'
            else:
                color = colors[ion_type]

            print( '{:4d}'.format(i_peak) + '{:10.4f}'.format(mz) + '{:10.1f}'.format(intensity*10000) + '  ' + interpretations_string + '   ' + ion_type + f"   {mz_delta}" )

            plot1.plot( [mz,mz], [0,intensity], color=color, linewidth=0.8 )
            blocked[int(mz)-1:int(mz)+1,0:int(intensity*100)] = 1

            annotation_string = ''
            match = re.match(r'(.+?)/', interpretations_string)
            if match:
                annotation_string = match.group(1)

            should_label = True
            match = re.search(r'\+[\d]?i', interpretations_string)
            if match:
                should_label = False
                #print(f"skip {interpretations_string}")

            if should_label:
                annotations.append( { 'mz': mz, 'intensity': intensity, 'annotation_string': annotation_string, 'color': color } )

            if color not in [ 'tab:gray', 'tab:olive']:
                markersize = 0.5
                show = True
                if intensity > 0.05:
                    markersize = intensity * 2.0 + 2.0
                match = re.search(r'\+[\d]?i', interpretations_string)
                if match:
                    show = False
                match = re.search(r'\-H2O', interpretations_string)
                if match:
                    show = False
                match = re.search(r'\-NH3', interpretations_string)
                if match and color not in [ 'tab:orange' ]:
                    show = False
                if show: 
                    plot2.plot( [mz,mz], [mz_delta,mz_delta], marker='s', markersize=markersize, color=color )
                    saved_residuals.append( { 'mz': mz, 'mz_delta': mz_delta, 'markersize': markersize, 'color': color } )

            match = re.match(r'([aby])([\d]+)/', interpretations_string)
            if match:
                series = match.group(1)
                ordinal = int(match.group(2))
                if series == 'y':
                    x = sequence_offset + ( len(residues) - ordinal - 0.45 ) * sequence_gap
                    y = sequence_height + 0.01
                    plot1.plot( [x,x,x+sequence_gap*0.2], [y,y-(intensity/10.0+0.005),y-(intensity/10.0+0.005)], color='tab:red')
                if series == 'b':
                    x = sequence_offset + ( ordinal - .5 ) * sequence_gap
                    y = sequence_height + 0.035
                    plot1.plot( [x,x,x-sequence_gap*0.2], [y,y+(intensity/10.0+0.005),y+(intensity/10.0+0.005)], color='tab:blue')

            counter += 1

        #### Sort all the annotations by intensity so that we annotate the most intense ions first and then only lower ones if there room
        annotations.sort(key=lambda x: x['intensity'], reverse=True)
        counter = 0
        for annotation in annotations:
            mz = annotation['mz']
            intensity = annotation['intensity']
            annotation_string = annotation['annotation_string']
            color = annotation['color']
            if blocked[int(mz-5.5):int(mz+5),int(intensity*100)+1:int(intensity*100+len(annotation_string)*1.5)].sum() == 0:
                #plot1.plot([int(mz-5.5), int(mz-5.5), int(mz+5), int(mz+5), int(mz-5.5)],
                #           [(int(intensity*100)+1)/100.0, (int(intensity*100)+len(annotation_string)*1.5)/100.0, (int(intensity*100+len(annotation_string)*1.5))/100.0,
                #            (int(intensity*100)+1)/100.0, (int(intensity*100)+1)/100.0], color='gray')
                plot1.text(mz, intensity  + 0.01, annotation_string, fontsize='x-small', ha='center', va='bottom', color=color, rotation=90, fontname='Arial')
                blocked[int(mz-5.5):int(mz+5),int(intensity*100)+1:int(intensity*100+len(annotation_string)*1.5)] = 1
            else:
                #print(f"{mz}\t{intensity*100}\t{annotation_string} is blocked")
                pass
            counter += 1

        #### Plot a little P where the precursor m/z is
        precursor_mz = None
        try:
            precursor_mz = spectrum.analytes['1']['precursor_mz']
        except:
            pass

        if precursor_mz:
            plot1.text(precursor_mz, -0.003, 'P', fontsize='small', ha='center', va='top', color='red', fontname='Arial')


        #### Set up the third plot, nominally for the precursor window
        gridspec3 = gridspec.GridSpec(1, 1)
        plot3 = fig.add_subplot(gridspec3[0])
        gridspec3.tight_layout(fig, rect=[0.02, -0.02, 1, 0.15])

        plot3.set_xlim([xmin, xmax])
        plot3.set_xticklabels([])
        plot3.set_ylim([-16,16])
        plot3.plot( [0,xmax], [0,0], '--', linewidth=0.6, color='gray')

        plt.savefig('AnnotatedSpectrum.pdf',format='pdf')
        plt.savefig('AnnotatedSpectrum.svg',format='svg')
        plt.close()

        #plt.show()



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
    argparser.add_argument('--annotate', action='count', help='If set, annotate the USI spectrum' )
    argparser.add_argument('--show_all_annotations', action='count', help='If set, show all the potential annotations, not just the final one' )
    argparser.add_argument('--examine', action='count', help='If set, examine the spectrum to see what can be learned' )
    argparser.add_argument('--score', action='count', help='If set, score the spectrum with the supplied peptidoform' )
    argparser.add_argument('--plot', action='count', help='If set, make a nice figure' )
    params = argparser.parse_args()

    # Set verbose mode
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    #### Flag for showing all annotations
    show_all_annotations = False
    if params.show_all_annotations is not None and params.show_all_annotations > 0:
        show_all_annotations = True

    #### If there is a JSON to load from
    #if params.input_json_filename:
    #    with open(params.input_json_filename) as infile:


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

    #### Create an annotator object
    annotator = SpectrumAnnotator()

    # Annotate the spectrum
    if params.annotate:
        annotator.annotate(spectrum, peptidoform=peptidoform, charge=charge, tolerance=params.tolerance)
        print(spectrum.show(show_all_annotations=show_all_annotations))
        if params.plot:
            annotator.plot(spectrum, peptidoform_string=peptidoform_string, charge=charge)

    # Examine the spectrum
    if params.score:
        annotator.compute_spectrum_score(spectrum, peptidoform=peptidoform, charge=charge)

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
        print(spectrum.show(show_all_annotations=show_all_annotations))


#### For command line usage
if __name__ == "__main__": main()
