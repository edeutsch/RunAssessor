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
fontname = 'FreeSans'

sys.path.append("C:\local\Repositories\GitHub\SpectralLibraries\lib")
from proforma_peptidoform import ProformaPeptidoform

from peptidoform import Peptidoform
from mass_reference import MassReference
from spectrum import Spectrum
from spectrum_examiner import SpectrumExaminer


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
    def annotate(self, spectrum, peptidoforms, charges, tolerance=None):

        if tolerance is None:
            tolerance = self.tolerance
        else:
            self.tolerance = tolerance

        examiner = SpectrumExaminer()

        spectrum.compute_spectrum_metrics()
        examiner.identify_isotopes(spectrum)
        examiner.identify_low_mass_ions(spectrum, peptidoforms[0])
        examiner.identify_reporter_ions(spectrum)
        #examiner.identify_neutral_losses(spectrum)
        examiner.identify_precursors(spectrum)
        self.annotate_peptidoform(spectrum, peptidoforms=peptidoforms, charges=charges)
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
        elif fragmentation_type == 'CID':
            series_list = [ 'a', 'b', 'y' ]
        elif fragmentation_type == 'ETD':
            series_list = [ 'c', 'z' ]
        elif fragmentation_type == 'EThcD':
            series_list = [ 'a', 'b', 'c', 'y', 'z' ]
        else:
            eprint("ERROR: Unrecognized fragmentation type")
            return
        base_series_score = { 'y': 90, 'z': 85, 'b': 80, 'c': 75, 'a': 70, 'm': 60 }

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
        if False:
            for frag in self.predicted_fragments_list:
                print(frag)
            exit()

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
    def annotate_peptidoform(self, spectrum, peptidoforms, charges, skip_internal_fragments=False, tolerance=None):

        #### If peptidoform is specified as None, then there's nothing we can do here
        if peptidoforms is None:
            return

        if tolerance is None:
            tolerance = self.tolerance
        else:
            self.tolerance = tolerance

        stripped_sequence = ''
        for peptidoform in peptidoforms:
            if not isinstance(peptidoform, ProformaPeptidoform):
                eprint(f"ERROR: Passed peptidoform is not an object, but instead is type {type(peptidoform)}")
                return
            stripped_sequence += peptidoform.peptide_sequence

        #### Extract user_parameters settings
        fragmentation_type = 'HCD'
        try:
            user_parameters = spectrum.extended_data['user_parameters']
        except:
            user_parameters = {}
        if 'dissociation_type' in user_parameters and user_parameters['dissociation_type'] is not None:
            fragmentation_type = user_parameters['dissociation_type']
            if fragmentation_type not in [ 'HCD', 'EThcD' ]:
                fragmentation_type = 'HCD'

        i_peptidoform = 0
        n_peptidoforms = len(peptidoforms)
        for peptidoform in peptidoforms:
            self.predict_fragment_ions(peptidoform=peptidoform, charge=charges[i_peptidoform], fragmentation_type=fragmentation_type, skip_internal_fragments=skip_internal_fragments)

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

                matches = self.find_close_predicted_fragments(mz, tolerance)
                if matches:
                    diagnostic_category = 'diagnostic'
                    for match in matches:
                        #peak[PL_INTERPRETATION_STRING] = f"{match[INT_INTERPRETATION_STRING]}/" + '{:.1f}'.format(match[INT_DELTA_PPM]) + 'ppm'
                        #peak[PL_INTERPRETATIONS].append(match)
                        if match[INT_INTERPRETATION_STRING].startswith('p'):
                            peak[PL_ATTRIBUTES][PLA_IS_PRECURSOR] += 1
                            diagnostic_category = 'nondiagnostic'

                        if n_peptidoforms > 1:
                            match[INT_INTERPRETATION_STRING] = f"{i_peptidoform+1}@{match[INT_INTERPRETATION_STRING]}"

                        self.add_interpretation(peak, match, diagnostic_category=diagnostic_category, residual_type='absolute')

            i_peptidoform += 1



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

        # We will return a list of possible matches
        matches = []

        # Get the integer mass as a dict key
        int_observed_mz = int(observed_mz)

        look_in_lower_bin_too = False
        look_in_higher_bin_too = False

        if observed_mz - int_observed_mz < tolerance:
            look_in_lower_bin_too = True
        if 1.0 - ( observed_mz - int_observed_mz ) < tolerance:
            look_in_higher_bin_too = True

        #### If there's just one bin to look in and it's not even in the index, then we're done
        if not look_in_lower_bin_too and not look_in_higher_bin_too:
            if int_observed_mz not in self.predicted_fragments_index:
                return matches
            else:
                predicted_fragment_list = self.predicted_fragments_index[int_observed_mz]

        #### But if we need to look in more than one bin, then create an enhanced predicted_fragment_list
        else:
            predicted_fragment_list = []
            if look_in_lower_bin_too:
                if int_observed_mz - 1 in self.predicted_fragments_index:
                    predicted_fragment_list.extend(self.predicted_fragments_index[int_observed_mz - 1])
            if int_observed_mz in self.predicted_fragments_index:
                predicted_fragment_list.extend(self.predicted_fragments_index[int_observed_mz])
            if look_in_higher_bin_too:
                if int_observed_mz + 1 in self.predicted_fragments_index:
                    predicted_fragment_list.extend(self.predicted_fragments_index[int_observed_mz + 1])

        # Loop over all the peaks in this bin and add them to matches if they're within tolerance
        for predicted_fragment in predicted_fragment_list:
            fragment_mz = predicted_fragment[0]
            delta = observed_mz - fragment_mz
            delta_ppm = delta / observed_mz * 1e6
            if delta_ppm > tolerance:
                continue
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

        if self.tolerance > 20:
            best_tolerance = self.tolerance * 0.5
            outer_tolerance = self.tolerance * 0.9
            max_tolerance = self.tolerance

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

                #### If this annotation is an isotope of another peak, but that peak is the TMT-related ion, then severely penalize this isotope possibility
                #### because they don't have isotopes and it's much more likely to be another reporter ion or report ion echo
                if interpretation[INT_INTERPRETATION_STRING].startswith('isotope'):
                    parent_peak = spectrum.peak_list[ peak[PL_ATTRIBUTES][PLA_PARENT_PEAK] ]
                    parent_annotation = parent_peak[PL_INTERPRETATION_STRING]
                    if 'TMT' in parent_annotation:
                        interpretation[INT_COMMONNESS_SCORE] *= 0.1

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
    def plot(self, spectrum, peptidoform, charge, xmin=None, xmax=None, mask_isolation_width=None, ymax=None, write_files=None):
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        import io

        #eprint(f"Font location: {matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')}")

        #### Extract user_parameters settings
        try:
            user_parameters = spectrum.extended_data['user_parameters']
        except:
            user_parameters = {}
        if 'xmin' in user_parameters and user_parameters['xmin'] is not None:
            try:
                xmin = float(user_parameters['xmin'])
            except:
                pass
        if 'xmax' in user_parameters and user_parameters['xmax'] is not None:
            try:
                xmax = float(user_parameters['xmax'])
            except:
                pass
        if 'ymax' in user_parameters and user_parameters['ymax'] is not None:
            try:
                ymax = float(user_parameters['ymax'])
            except:
                pass

        show_sequence = True
        if 'show_sequence' in user_parameters and user_parameters['show_sequence'] is not None:
            if not user_parameters['show_sequence']:
                show_sequence = False
        show_b_and_y_flags = True
        if 'show_b_and_y_flags' in user_parameters and user_parameters['show_b_and_y_flags'] is not None:
            if not user_parameters['show_b_and_y_flags']:
                show_b_and_y_flags = False

        include_third_plot = False
        figure_height = 6
        figure_width = 8
        spectrum_viewport = [0.02, 0.2, 1.01, 1.01]
        residuals_viewport = [0.02, -0.04, 1.01, 0.21]
        if include_third_plot:
            figure_height = 7.5
            spectrum_viewport = [0.02, 0.28, 1.01, 1.01]
            residuals_viewport = [0.02, 0.10, 1.01, 0.28]
            third_plot_viewport = [0.02, -0.03, 1.01, 0.15]

        fig = plt.figure(figsize=(figure_width, figure_height))
        gridspec1 = gridspec.GridSpec(1, 1)
        plot1 = fig.add_subplot(gridspec1[0])
        gridspec1.tight_layout(fig, rect=spectrum_viewport)

        #### Try to get the precursor_mz
        precursor_mz = None
        try:
            precursor_mz = spectrum.analytes['1']['precursor_mz']
        except:
            pass

        #### Compute the min and max for mz and intensities
        max_intensity = 0
        max_non_precursor_intensity = 0
        min_mz = 9999
        max_mz = 0
        for peak in spectrum.peak_list:
            mz = peak[PL_MZ]
            intensity = peak[PL_INTENSITY]
            interpretations_string = peak[PL_INTERPRETATION_STRING]
            if precursor_mz is not None and mask_isolation_width is not None and mz > precursor_mz - mask_isolation_width/2 and mz < precursor_mz + mask_isolation_width/2:
                continue
            if intensity > max_intensity:
                max_intensity = intensity
            match = re.match(r'\d+\@p', interpretations_string)
            if len(interpretations_string) > 0 and interpretations_string[0] != 'p' and interpretations_string[0] != 'r' and not match and intensity > max_non_precursor_intensity:
                max_non_precursor_intensity = intensity
            if mz < min_mz:
                min_mz = mz
            if mz > max_mz:
                max_mz = mz
        if xmin is None:
            xmin = int(min_mz) - 15
        else:
            xmin = int(xmin + 0.5)
        if xmax is None:
            xmax = int(max_mz) + 15
        else:
            xmax = int(xmax + 0.5)
        xscale = (xmax - xmin) / 1200.0

        if ymax is None or ymax < 0.01:
            try:
                ymax = float(ymax)
            except:
                ymax = 1.075

        #### Set up the main spectrum plot
        #plot1.plot( [0,1000], [0,1], color='tab:green' )
        plot1.set_xlabel('m/z', fontname=fontname)
        plot1.set_ylabel('Relative Intensity', fontname=fontname)
        plot1.set_xlim([xmin, xmax])
        plot1.set_ylim([0,ymax])
        #ax[0,1].plot( [0,limit], [0,limit], '--', linewidth=1, color='gray')

        #### Set up the residuals plot
        gridspec2 = gridspec.GridSpec(1, 1)
        plot2 = fig.add_subplot(gridspec2[0])
        gridspec2.tight_layout(fig, rect=residuals_viewport)

        #plot2.plot( [0,1000], [-10,10], color='tab:green' )
        #plot2.set_xlabel('m/z')
        plot2.set_ylabel('delta (PPM)', fontname=fontname)
        plot2.set_xlim([xmin, xmax])
        plot2.set_xticklabels([])
        plot2.set_ylim([-16,16])
        plot2.plot( [0,xmax], [0,0], '--', linewidth=0.6, color='gray')

        #### Set up colors for different types of ion and a grid to track where items have been layed out
        colors = { 'b': 'tab:blue', 'a': 'tab:green', 'y': 'tab:red', '0': 'violet', '_': 'tab:gray', 'I': 'gold', '?': 'tab:gray', 'p': 'tab:pink', 'm': 'tab:brown', 'r': 'tab:purple', 'f': 'tab:purple', 'c': 'tab:orange', 'z': 'c' }
        blocked = np.zeros((xmax,100))

        #### Prepare to write the peptide sequence to the plot, although only write it later once we figure out where there's room
        stripped_peptide = peptidoform.peptide_sequence
        residues = list(stripped_peptide)
        modified_residues = peptidoform.residue_modifications
        #sequence_gap = (xmax-xmin)/45 # for x-large font
        sequence_gap = (xmax-xmin)/55
        sequence_height = 0.85 * ymax
        sequence_offset = xmin + 2 * sequence_gap
        original_sequence_offset = sequence_offset

        #### Loop over all peaks and plot them in the color that they would be annotated with
        counter = 0
        annotations = []
        saved_residuals = []
        all_peaks = []
        all_flags = []
        for peak in spectrum.peak_list:
            i_peak = peak[PL_I_PEAK]
            mz = peak[PL_MZ]
            intensity = peak[PL_INTENSITY] / max_intensity
            interpretations_string = peak[PL_INTERPRETATION_STRING]
            if precursor_mz is not None and mask_isolation_width is not None and mz > precursor_mz - mask_isolation_width/2 and mz < precursor_mz + mask_isolation_width/2:
                continue

            mz_delta = 99
            match = re.search(r'/([\+\-\d\.]+)ppm', interpretations_string)
            if match:
                mz_delta = float(match.group(1))
                mz_delta += spectrum.attributes['mass_accuracy']['offset']

            try:
                ion_type = interpretations_string[0]
            except:
                ion_type = '?'

            if len(interpretations_string) > 1 and interpretations_string[1] == '@':
                ion_type = interpretations_string[2]
            if ion_type in [ '1','2','3','4','5','6','7','8','9' ]:
                color = 'tab:olive'
            elif ion_type in colors:
                color = colors[ion_type]
            else:
                color = 'black'

            #print( '{:4d}'.format(i_peak) + '{:10.4f}'.format(mz) + '{:10.1f}'.format(intensity*10000) + '  ' + interpretations_string + '   ' + ion_type + f"   {mz_delta}" )

            #### Store the peak for later plotting
            all_peaks.append({ 'mz': mz, 'intensity': intensity, 'color': color } )

            annotation_string = interpretations_string
            match = re.match(r'(.+?)/', interpretations_string)
            if match:
                annotation_string = match.group(1)

            should_label = True
            match = re.search(r'\+[\d]?i', annotation_string)
            if match:
                should_label = False
                #print(f"skip {interpretations_string}")

            if annotation_string.startswith('?'):
                if 'show_unknown' in user_parameters and user_parameters['show_unknown']:
                    pass
                else:
                    should_label = False

            if should_label:
                annotation_priority = 1
                if ( annotation_string.startswith('y') or annotation_string.startswith('b') ) and annotation_string.count('-') == 0:
                    annotation_priority = 2
                annotations.append( { 'mz': mz, 'intensity': intensity, 'annotation_string': annotation_string, 'color': color, 'annotation_priority': annotation_priority } )

            #### For certain types of ions, record a massdiff residual to plot later
            if color not in [ 'tab:gray', 'tab:olive']:
                markersize = 0.5
                show = True
                if intensity > 0.05:
                    markersize = intensity * 2.0 + 2.0
                match = re.search(r'\+[\d]?i', interpretations_string)
                if match:
                    show = False
                #match = re.search(r'\-H2O', interpretations_string)
                #if match:
                #    show = False
                #match = re.search(r'\-NH3', interpretations_string)
                #if match and color not in [ 'tab:orange' ]:
                #    show = False
                if show: 
                    plot2.plot( [mz,mz], [mz_delta,mz_delta], marker='s', markersize=markersize, color=color )
                    saved_residuals.append( { 'mz': mz, 'mz_delta': mz_delta, 'markersize': markersize, 'color': color } )

            #### Decorate the sequence with little flags to indicate ion strength
            match = re.match(r'([aby])([\d]+)(-.+)?(\^\d)?', interpretations_string)
            if match:
                #print(f"+++{interpretations_string}")
                series = match.group(1)
                ordinal = int(match.group(2))
                flag_direction = 1.0
                flag_thickness = 0.9
                loss = match.group(3)
                flag_intensity = peak[PL_INTENSITY] / max_non_precursor_intensity
                if loss is not None:
                    flag_direction = -1.0
                    flag_thickness = 0.5
                if series == 'y':
                    x = sequence_offset + ( len(residues) - ordinal - 0.45 ) * sequence_gap - sequence_gap*0.02
                    y = sequence_height + 0.007 * ymax
                    all_flags.append( [ 'y', x, y, flag_intensity, 'tab:red', flag_direction, flag_thickness ] )
                    #plot1.plot( [x,x,x+sequence_gap*0.2], [y,y-(flag_intensity/10.0+0.005),y-(flag_intensity/10.0+0.005)], color='tab:red', linewidth=0.9)
                if series == 'b':
                    x = sequence_offset + ( ordinal - .5 ) * sequence_gap
                    y = sequence_height + 0.037 * ymax
                    all_flags.append( [ 'b', x, y, flag_intensity, 'tab:blue', flag_direction, flag_thickness ] )
                    #plot1.plot( [x,x,x-sequence_gap*0.2], [y,y+(flag_intensity/10.0+0.005),y+(flag_intensity/10.0+0.005)], color='tab:blue', linewidth=0.9)

            counter += 1

        #### First sort all the peaks in smallest to tallest and plot in that order
        all_peaks.sort(key=lambda x: x['intensity'], reverse=False)
        for peak in all_peaks:
            mz = peak['mz']
            intensity = peak['intensity']
            color = peak['color']
            plot1.plot( [mz,mz], [0,intensity], color=color, linewidth=0.6 )
            blocked[int(mz)-1:int(mz)+1,0:int(intensity*100)] = 1

        #### Sort all the annotations by intensity so that we annotate the most intense ions first and then only lower ones if there room
        annotations.sort(key=lambda x: (x['annotation_priority'], x['intensity']), reverse=True)
        xf = 1.8
        counter = 0
        for annotation in annotations:
            mz = annotation['mz']
            intensity = annotation['intensity']
            annotation_string = annotation['annotation_string']
            if annotation_string == '':
                continue
            if mz < xmin or mz > xmax:
                continue
            color = annotation['color']
            blocklen = len(annotation_string)
            if blocklen < 3:
                blocklen = 3
            #print(f"-- Try to annotate {mz}\t{intensity}\t{annotation_string}")
            if blocked[int(mz-7.5*xscale):int(mz+7*xscale),int(intensity*100)+1:int(intensity*100+blocklen*xf)].sum() == 0:
                #plot1.plot([int(mz-7.5*xscale), int(mz-7.5*xscale), int(mz+7*xscale), int(mz+7*xscale), int(mz-7.5*xscale)],
                #           [(int(intensity*100)+1)/100.0, (int(intensity*100+blocklen*xf))/100.0, (int(intensity*100+blocklen*xf))/100.0,
                #            (int(intensity*100)+1)/100.0, (int(intensity*100)+1)/100.0], color='gray', linewidth=0.2)
                plot1.text(mz, intensity  + 0.01*ymax, annotation_string, fontsize='x-small', ha='center', va='bottom', color=color, rotation=90, fontname=fontname)
                #blocked[int(mz-5.5):int(mz+5),int(intensity*100)+1:int(intensity*100+blocklen*1.5)] = 1
                blocked[int(mz-7.5*xscale):int(mz+7*xscale),int(intensity*100)+1:int(intensity*100+blocklen*xf)] = 1
                #print(f"   - easy")
            else:
                #### If there are more than 2 losses, don't work hard to squeeze it in
                if annotation_string.count('-') > 2:
                    continue
                #print(f"{mz}\t{intensity*100}\t{annotation_string} is blocked")
                reposition_attempts = [ [0.0, 3.0], [0.0, 6.0], [0.0, 10.0], [0.0, 15.0], [5.0, 3.0], [-5.0, 3.0], [5.0, 15.0], [-5.0, 15.0] ]
                found_a_spot = False
                for reposition_attempt in reposition_attempts:
                    x_offset = reposition_attempt[0]
                    y_offset = reposition_attempt[1]
                    if blocked[int(mz-7.5*xscale+x_offset):int(mz+7*xscale+x_offset),int(intensity*100+y_offset)+1:int(intensity*100+blocklen*xf+y_offset)].sum() == 0:
                        plot1.text(mz+x_offset, intensity + 0.01*ymax + y_offset/100.0*ymax, annotation_string, fontsize='x-small', ha='center', va='bottom', color=color, rotation=90, fontname=fontname)
                        blocked[int(mz-7.5*xscale+x_offset):int(mz+7*xscale+x_offset),int(intensity*100+y_offset)+1:int(intensity*100+blocklen*xf+y_offset)] = 1
                        plot1.plot( [mz,mz+x_offset], [intensity + 0.004*ymax, intensity + 0.01*ymax + (y_offset-0.5)/100.0*ymax], color='black', linewidth=0.2 )
                        #print(f"   - managed to find a spot at {reposition_attempt}")
                        found_a_spot = True
                        break
                if not found_a_spot:
                    #print(f"   - Failed to find a place")
                    pass

            counter += 1

        #### Plot a little P where the precursor m/z is
        if precursor_mz:
            plot1.text(precursor_mz, -0.003 * ymax, 'P', fontsize='small', ha='center', va='top', color='red', fontname=fontname)
        plot1.spines[['right', 'top']].set_visible(False)
        #plot2.spines[['right', 'top']].set_visible(False)

        #### Look for a clear spot to put the sequence
        if show_sequence is True:
            done = False
            x_extent = sequence_offset
            leftest_offset = None
            rightest_offset = None
            while not done:
    
                running_sum = 0.0
                for counter in range(len(residues)+2):
                    x_offset = sequence_offset + counter * sequence_gap
                    x_width = 0.5 * sequence_gap
                    y_offset = sequence_height
                    y_width = 0.07 * ymax
                    #plot1.plot([x_offset-x_width, x_offset-x_width, x_offset+x_width, x_offset+x_width, x_offset-x_width],
                    #        [y_offset-y_width, y_offset+y_width*1.5, y_offset+y_width*1.5, y_offset-y_width, y_offset-y_width], color='gray', linewidth=0.2)
                    running_sum += blocked[int(x_offset-x_width):int(x_offset+x_width),int((y_offset-y_width)*100):int((y_offset+y_width*1.5)*100)].sum()
                    x_extent = x_offset + x_width

                if running_sum == 0:
                    #print(f"*** Found space at sequence_offset={sequence_offset}, x_extent={x_extent}")
                    if leftest_offset is None:
                        leftest_offset = sequence_offset
                else:
                    #print(f"*** Blocked at sequence_offset={sequence_offset}, x_extent={x_extent}")
                    if leftest_offset is not None:
                        rightest_offset = sequence_offset
                        sequence_offset = ( leftest_offset + rightest_offset ) / 2.0
                        break
                sequence_offset += (xmax-xmin)/50.0

                if x_extent > xmax*.97:
                    #print(f"*** Reached the limit at sequence_offset={sequence_offset}, x_extent={x_extent}")
                    if leftest_offset is None:
                        sequence_offset -= (xmax-xmin)/50.0
                        sequence_offset = original_sequence_offset
                    else:
                        rightest_offset = sequence_offset
                        sequence_offset = ( leftest_offset + rightest_offset ) / 2.0
                    break

            #### Finally write out the sequence
            counter = 0
            if peptidoform.terminal_modifications is not None and 'nterm' in peptidoform.terminal_modifications:
                plot1.text(sequence_offset + 0*sequence_gap, sequence_height, '=', fontsize='large', ha='center', va='bottom', color='tab:orange', fontname=fontname)
                sequence_offset += sequence_gap * 0.8
            for residue in residues:
                #plot1.text(sequence_offset+counter*sequence_gap, sequence_height, residue, fontsize='x-large', ha='center', va='bottom', color='black', fontname=fontname)
                color = 'black'
                if modified_residues is not None and counter + 1 in modified_residues:
                    color= 'tab:orange'
                plot1.text(sequence_offset+counter*sequence_gap, sequence_height, residue, fontsize='large', ha='center', va='bottom', color=color, fontname=fontname)
                counter += 1
            plot1.text(sequence_offset + (counter+.2)*sequence_gap, sequence_height + 0.02*ymax, f"{charge}+", fontsize='medium', ha='center', va='bottom', color='black', fontname=fontname)

            #### Finally paint the flags
            if show_b_and_y_flags is True:
                for flag in all_flags:
                    series, x, y, intensity, color, flag_direction, flag_thickness = flag
                    x += ( sequence_offset - original_sequence_offset)
                    if series == 'y':
                        plot1.plot( [x,x,x+sequence_gap*0.2*flag_direction], [y,y-(intensity/10.0+0.005)*ymax,y-(intensity/10.0+0.005)*ymax], color='tab:red', linewidth=flag_thickness)
                    if series == 'b':
                        plot1.plot( [x,x,x-sequence_gap*0.2*flag_direction], [y,y+(intensity/10.0+0.005)*ymax,y+(intensity/10.0+0.005)*ymax], color='tab:blue', linewidth=flag_thickness)

        #with open('saved_residuals_calibrated.json', 'w') as outfile:
        #    outfile.write(json.dumps(saved_residuals))

        #### Set up the third plot, nominally for the precursor window
        if include_third_plot:
            gridspec3 = gridspec.GridSpec(1, 1)
            plot3 = fig.add_subplot(gridspec3[0])
            gridspec3.tight_layout(fig, rect=third_plot_viewport)

            plot3.set_ylabel('delta (PPM)', fontname=fontname)
            plot3.set_xlim([xmin, xmax])
            plot3.set_xticklabels([])
            plot3.set_ylim([-16,16])
            plot3.plot( [0,xmax], [0,0], '--', linewidth=0.6, color='gray')

            with open('saved_residuals_calibrated.json') as infile:
                saved_residuals = json.load(infile)
            for residual in saved_residuals:
                mz = residual['mz']
                mz_delta = residual['mz_delta']
                markersize = residual['markersize']
                color = residual['color']
                plot3.plot( [mz,mz], [mz_delta,mz_delta], marker='s', markersize=markersize, color=color )

            plot1.text(xmax-10, 0.98, 'A', fontname=fontname, fontsize=30, ha='right', va='top')
            plot2.text(xmax-6, 14, 'B', fontname=fontname, fontsize=20, ha='right', va='top')
            plot3.text(xmax-6, 14, 'C', fontname=fontname, fontsize=20, ha='right', va='top')


        #### Write out the figure to PDF and SVG
        if 'create_svg' in user_parameters and user_parameters['create_svg']:
            buffer = io.BytesIO()
            plt.savefig(buffer,format='svg')
            buffer.seek(0)
            content = buffer.read()
            spectrum.extended_data['svg'] = content.decode('utf-8')
        if 'create_pdf' in user_parameters and user_parameters['create_pdf']:
            buffer = io.BytesIO()
            plt.savefig(buffer,format='pdf')
            buffer.seek(0)
            content = buffer.read()
            spectrum.extended_data['pdf'] = content.decode('iso-8859-1')

        if write_files is not None:
            plt.savefig('AnnotatedSpectrum.pdf',format='pdf')
            plt.savefig('AnnotatedSpectrum.svg',format='svg')

        plt.close()



####################################################################################################
#### Gaussian function used during curve fitting procedures
def gaussian_function(x, a, x0, sigma):
    return a * exp( -(x-x0)**2 / (2*sigma**2) )


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Class representing a peptidoform')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('--usi', action='store', default=None, help='USI to process')
    argparser.add_argument('--input_json_filename', action='store', default=None, type=float, help='Filename of an input json file')
    argparser.add_argument('--tolerance', action='store', default=None, type=float, help='Tolerance in ppm for annotation')
    argparser.add_argument('--annotate', action='count', help='If set, annotate the USI spectrum' )
    argparser.add_argument('--examine', action='count', help='If set, examine the spectrum to see what can be learned' )
    argparser.add_argument('--simplify', action='count', help='If set, simplify the spectrum as much as possible' )
    argparser.add_argument('--score', action='count', help='If set, score the spectrum with the supplied peptidoform' )
    argparser.add_argument('--show_all_annotations', action='count', help='If set, show all the potential annotations, not just the final one' )
    argparser.add_argument('--plot', action='count', help='If set, make a nice figure' )
    argparser.add_argument('--write_files', action='count', help='If set, write the figures to files' )
    argparser.add_argument('--xmin', action='store', default=None, type=float, help='Set a manual x-axis (m/z) minimum' )
    argparser.add_argument('--xmax', action='store', default=None, type=float, help='Set a manual x-axis (m/z) maximum' )
    argparser.add_argument('--mask_isolation_width', action='store', default=None, type=float, help='When plotting, drop peaks within an isolation window with this full width' )
    argparser.add_argument('--ymax', action='store', default=None, type=float, help='Set a new ymax in order to compensate for very tall peaks (e.g, 0.5 or 1.2)' )
    argparser.add_argument('--show_sequence', action='store', default=True, type=str2bool, help='Set to false to suppress the peptide and its flags' )
    argparser.add_argument('--show_b_and_y_flags', action='store', default=True, type=str2bool, help='Set to false to suppress the peptide sequence flags' )
    params = argparser.parse_args()

    # Set verbose mode
    verbose = params.verbose
    if verbose is None:
        verbose = 0

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
        #print(json.dumps(usi.peptidoforms, indent=2))
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
    peptidoforms = []
    for usi_peptidoform in usi.peptidoforms:
        if usi_peptidoform['peptidoform_string'] is not None and usi_peptidoform['peptidoform_string'] != '':
            peptidoform = ProformaPeptidoform(usi_peptidoform['peptidoform_string'])
            peptidoforms.append(peptidoform)
        print(json.dumps(peptidoform.to_dict(),indent=2,sort_keys=True))

    # Annotate the spectrum
    if params.annotate:
        annotator = SpectrumAnnotator()
        spectrum.extended_data['user_parameters'] = {}
        spectrum.extended_data['user_parameters']['show_sequence'] = params.show_sequence
        spectrum.extended_data['user_parameters']['show_b_and_y_flags'] = params.show_b_and_y_flags
        print(json.dumps(spectrum.extended_data['user_parameters'], indent=2))
        annotator.annotate(spectrum, peptidoforms=peptidoforms, charges=usi.charges, tolerance=params.tolerance)
        print(spectrum.show(show_all_annotations=show_all_annotations, verbose=verbose))
        if params.plot:
            annotator.plot(spectrum, peptidoform=peptidoforms[0], charge=usi.charges[0], xmin=params.xmin, xmax=params.xmax,
                           mask_isolation_width=params.mask_isolation_width, ymax=params.ymax, write_files=params.write_files)

    # Score the spectrum
    if params.score:
        annotator.compute_spectrum_score(spectrum, peptidoform=peptidoforms[0], charge=usi.charges[0])


#### For command line usage
if __name__ == "__main__": main()
