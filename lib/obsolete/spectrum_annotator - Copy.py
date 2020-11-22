#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import re
import itertools

from peptidoform import Peptidoform
from mass_reference import MassReference


####################################################################################################
#### SpectrumAnnotator class
class SpectrumAnnotator:

    ####################################################################################################
    #### Constructor
    def __init__(self, peptidoform=None, verbose=0):

        # Set verbosity
        if verbose is None:
            verbose = 0
        self.verbose = verbose

        self.mzs_list = []
        self.interpretations_list = []
        self.spectrum_attributes = {}

        self.mass_reference = MassReference()

        self.peptidoform = peptidoform

        if peptidoform is not None:
            self.predict_fragment_ions()


    ####################################################################################################
    #### Predict all the fragment ions for the provided peptidoform
    def predict_fragment_ions(self, peptidoform=None, charge=1, fragmentation_type='HCD'):

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

        # Get handles for some needed reference masses
        masses = self.mass_reference.atomic_masses
        residue_masses = self.mass_reference.aa_masses
        ion_series_attr = self.mass_reference.ion_series_attributes
        neutral_losses = self.mass_reference.neutral_losses
        neutral_losses_by_residue = self.mass_reference.neutral_losses_by_residue
        neutral_losses_by_formula = self.mass_reference.neutral_losses_by_formula
        terminal_modifications = self.mass_reference.terminal_modifications

        # Determine the terminal modification masses
        terminal_mass_modifications = { 'nterm': 0.0, 'cterm': 0.0 }
        if peptidoform.nterm['name'] != '':
            if peptidoform.nterm['name'] in terminal_modifications:
                terminal_mass_modifications['nterm'] = terminal_modifications[peptidoform.nterm['name']]['mass']
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
            potential_losses[series] = {}
            cumulative_residues[series] = 0

        # Initialize a theoretical spectrum
        self.mzs_list = []
        self.interpretations_list = []
        peptide_length = len(peptidoform.residues)
        all_annotations = {}

        # Main loop: iterate through each position, working from both ends simultaneously
        for i_residue in range(peptide_length-1):

            # Add additional series for internal ions
            if i_residue > 0 and i_residue < peptide_length-2:
                series_name = f"m{i_residue+1}"
                series_list.append(series_name)
                series_type = 'm'
                cumulative_mass[series_name] = ion_series_attr[series_type]['mass'] + terminal_mass_modifications[ion_series_attr[series_type]['terminus_type']]
                potential_losses[series_name] = {}
                cumulative_residues[series_name] = 0

                # Also add the very common "a" type ion as a CO neutral loss
                loss_type = 'CO'
                potential_losses[series_name][loss_type] = 1

            # Generate fragments for each ion series we expect to find
            for series in series_list:

                series_type = series[0]
                cumulative_residues[series] += 1
                #print(series)

                # And for each expected charge
                for i_charge in range(1,charge+1):

                    # Get the next residue
                    if ion_series_attr[series_type]['terminus_type'] == 'nterm':
                        residue = peptidoform.residues[i_residue]['string']
                    else:
                        residue = peptidoform.residues[peptide_length - 1 - i_residue]['string']

                    # If this residue is not recognized, this is a serious error
                    if residue not in residue_masses:
                        eprint(f"ERROR: Unrecognized residue '{residue}'")
                        return

                    # Only compute certain things on the singly charged pass
                    if i_charge == 1:

                        # Update the cumulative mass
                        cumulative_mass[series] += residue_masses[residue]

                        # See if this residue can yield a neutral loss and store it if so
                        if residue in neutral_losses_by_residue:
                            loss_type = neutral_losses_by_residue[residue]['formula']
                            if loss_type not in potential_losses[series]:
                                potential_losses[series][loss_type] = 0
                            potential_losses[series][loss_type] += 1
                            #print(f"Adding an instance of {loss_type}")

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
                            #print(f"    {i_loss}")
                            if potential_neutral_loss_formula is not None:
                                loss_string += f"-{potential_neutral_loss_formula}"
                                loss_mass += neutral_losses_by_formula[potential_neutral_loss_formula]['delta_mass']

                        # Final m/z computation and store the interpretation
                        interpretation = f"{series}{i_residue+1}{loss_string}"
                        if series_type == 'm':
                            if cumulative_residues[series] > 1:
                                interpretation = f"{series}:{i_residue+1}{loss_string}"
                            else:
                                continue
                        if i_charge > 1:
                            interpretation += f"^{i_charge}"

                        # Avoid duplicate annotations when different permutations lead to the same thing
                        if interpretation not in all_annotations:
                            self.mzs_list.append( ( cumulative_mass[series] - loss_mass + masses['proton'] * i_charge ) / i_charge )
                            self.interpretations_list.append(interpretation)
                            all_annotations[interpretation] = 1




    ####################################################################################################
    #### Return a printable buffer string of the details of the peptidoform
    def show(self):

        buf = ''
        buf += f"Peptidoform_string={self.peptidoform.peptidoform_string}\n"
        buf += f"Charge={self.spectrum_attributes['charge']}\n"
        for i_peak in range(len(self.mzs_list)):
            buf += '  ' + '{:10.4f}'.format(self.mzs_list[i_peak]) + f": {self.interpretations_list[i_peak]}\n"
        return buf


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Class representing a peptidoform')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
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

    #### Loop over all the peptidoform strings and decompose them
    for peptidoform_string in peptidoform_strings_list:
        print('**************************************************************************')
        peptidoform = Peptidoform(peptidoform_string, verbose=verbose)
        print(peptidoform.show())
        annotator = SpectrumAnnotator()
        annotator.predict_fragment_ions(peptidoform)
        print(annotator.show())

#### For command line usage
if __name__ == "__main__": main()
