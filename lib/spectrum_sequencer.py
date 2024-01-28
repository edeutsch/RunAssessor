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
#### SpectrumSequencer class
class SpectrumSequencer:
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
    def remove_low_intensity_peaks(self, spectrum):

        n_peaks_to_keep_per_window = 7
        window_width = 100

        n_peaks = spectrum.attributes['number of peaks']
        new_peak_list = []
        tmp_peak_list = []
        prev_cutoff = None
        next_cutoff = None
        n_peaks_in_bin = 0
        intensity_threshold = 0
        debug = True

        # Loop through and keep the top N peaks
        for i_peak in range(n_peaks):

            mz = spectrum.peak_list[i_peak][PL_MZ]
            intensity = spectrum.peak_list[i_peak][PL_INTENSITY]

            if i_peak == 0:
                next_cutoff = int(mz/window_width) * window_width + window_width
                prev_cutoff = next_cutoff - window_width

            if mz > next_cutoff:
                if debug: print(f"Bin {prev_cutoff} thru {next_cutoff} has {n_peaks_in_bin} peaks")
                if n_peaks_in_bin == 0:
                    pass
                elif n_peaks_in_bin < n_peaks_to_keep_per_window:
                    new_peak_list.extend(tmp_peak_list)
                    if debug: print(f"  - keeping all {n_peaks_in_bin} peaks")
                else:
                    sorted_peaks = sorted(tmp_peak_list,key=lambda x: x[PL_INTENSITY], reverse=True)
                    top_peaks = sorted_peaks[0:n_peaks_to_keep_per_window]
                    resorted_peaks = sorted(top_peaks,key=lambda x: x[PL_MZ], reverse=False)
                    new_peak_list.extend(resorted_peaks)
                    if debug:
                        print(f"  - keeping top {len(resorted_peaks)} peaks")
                        #for ii in range(len(top_peaks)):
                        #    print(top_peaks[ii])
                tmp_peak_list = []
                prev_cutoff = next_cutoff
                next_cutoff += window_width
                n_peaks_in_bin = 0
                intensity_threshold = 0

            tmp_peak_list.append(spectrum.peak_list[i_peak])
            n_peaks_in_bin += 1

        if debug: print(f"Final bin thru {next_cutoff} has {n_peaks_in_bin} peaks")
        if n_peaks_in_bin == 0:
            pass
        elif n_peaks_in_bin < n_peaks_to_keep_per_window:
            new_peak_list.extend(tmp_peak_list)
            if debug: print(f"  - keeping all {n_peaks_in_bin} peaks")
        else:
            sorted_peaks = sorted(tmp_peak_list,key=lambda x: x[PL_INTENSITY], reverse=True)
            top_peaks = sorted_peaks[0:n_peaks_to_keep_per_window]
            resorted_peaks = sorted(top_peaks,key=lambda x: x[PL_MZ], reverse=False)
            new_peak_list.extend(resorted_peaks)
            if debug: print(f"  - keeping top {len(resorted_peaks)} peaks")

        spectrum.peak_list = new_peak_list
        spectrum.attributes['number of peaks'] = len(spectrum.peak_list)


    
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
    def create_peak_network(self, spectrum, sequencing_parameters):

        # Get the number of peaks and ensure that there is an index
        n_peaks = spectrum.attributes['number of peaks']

        #### Extract some sequence parameters
        precursor_mz = sequencing_parameters['precursor_mz']
        precursor_charge = sequencing_parameters['precursor_charge']
        tolerance = sequencing_parameters['tolerance']
        fragmentation_type = sequencing_parameters['fragmentation_type']

        verbose = 0

        #### If some labels were provided in the sequencing_parameters, then add them to the mass reference
        if 'labels' in sequencing_parameters and sequencing_parameters['labels'] is not None and len(sequencing_parameters['labels']) > 0:
            self.mass_reference.prepare_mass_tables(labels=sequencing_parameters['labels'])

        # Get handles for some needed reference masses
        masses = self.mass_reference.atomic_masses
        residue_masses = self.mass_reference.nr_aa_masses
        ion_series_attr = self.mass_reference.ion_series_attributes
        neutral_losses = self.mass_reference.neutral_losses
        neutral_losses_by_residue = self.mass_reference.neutral_losses_by_residue
        neutral_losses_by_formula = self.mass_reference.neutral_losses_by_formula

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

            #### Ignore deleted peaks
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_DELETED]:
                continue

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

        #### Add a node for bottom (i.e. 0)
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

        #### Add a node for the top (i.e. the precursor)
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
            #### Ignore deleted peaks
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_DELETED]:
                continue
            self.link_nodes(spectrum, sequencing_parameters, start_i_peak=i_peak, direction='upward')
            self.link_nodes(spectrum, sequencing_parameters, start_i_peak=i_peak, direction='downward')
        self.link_nodes(spectrum, sequencing_parameters, start_i_peak='bottom', direction='upward')
        self.link_nodes(spectrum, sequencing_parameters, start_i_peak='top', direction='downward')


        # Show some stuff
        if False:
            for node_id,node in spectrum.network['nodes'].items():
                print(f"{node_id}\t{node['mz']}\t{node['intensity']}")
                print(json.dumps(node, indent=2, sort_keys=True))
        if False:
            print(json.dumps(spectrum.network, indent=2, sort_keys=True))



    ####################################################################################################
    def interpret_de_novo(self, spectrum, sequencing_parameters):

        # Get the number of peaks and ensure that there is an index
        n_peaks = spectrum.attributes['number of peaks']

        #### Extract some sequence parameters
        precursor_mz = sequencing_parameters['precursor_mz']
        precursor_charge = sequencing_parameters['precursor_charge']
        tolerance = sequencing_parameters['tolerance']
        fragmentation_type = sequencing_parameters['fragmentation_type']

        verbose = 0

        #### If some labels were provided in the sequencing_parameters, then add them to the mass reference
        if 'labels' in sequencing_parameters and sequencing_parameters['labels'] is not None and len(sequencing_parameters['labels']) > 0:
            self.mass_reference.prepare_mass_tables(labels=sequencing_parameters['labels'])

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

            #### Ignore deleted peaks
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_DELETED]:
                continue

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

        #### Add a node for bottom (i.e. 0)
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

        #### Add a node for the top (i.e. the precursor)
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
            #### Ignore deleted peaks
            if spectrum.peak_list[i_peak][PL_ATTRIBUTES][PLA_IS_DELETED]:
                continue
            self.link_nodes(spectrum, sequencing_parameters, start_i_peak=i_peak, direction='upward')
            self.link_nodes(spectrum, sequencing_parameters, start_i_peak=i_peak, direction='downward')
        self.link_nodes(spectrum, sequencing_parameters, start_i_peak='bottom', direction='upward')
        self.link_nodes(spectrum, sequencing_parameters, start_i_peak='top', direction='downward')


        # Show some stuff
        if True:
            for node_id,node in spectrum.network['nodes'].items():
                print(f"{node_id}\t{node['mz']}\t{node['intensity']}")
                print(json.dumps(node, indent=2, sort_keys=True))
            #self.plot_network(spectrum)
            #exit()

        verbose = 0

        #### Set up some data containers
        active_sequences = []       # list of sequences that are still actively being extended
        finished_sequences = []     # list of sequences that cannot be extended any farther
        followed_peaks = {}         # list of peaks that have already been processed

        #### Set up containers for the tracking of progress that has already been made from each terminus
        sequencing_checkpoints = {
            'nterm': { 'list': [], 'dict': {} },
            'cterm': { 'list': [], 'dict': {} },
            'internal': { 'list': [], 'dict': {} }
        }

        #### Loop over all nodes (peaks) and start trying to sequence from there
        nodes = list(spectrum.network['nodes'].keys())
        sorted_nodes = sorted(nodes, key=lambda x: spectrum.network['nodes'][x]['sort_order'])
        for node_id in sorted_nodes:
            node = spectrum.network['nodes'][node_id]
            if verbose: print(f"Processing node {node_id}")

            done = False
            active_sequences = [ { 'score': 0.0, 'current_node_id': node_id, 'sequence': [], 'residues': [], 'score_per_100mz': 0.0, 'terminus_type': '', 'mass': 0.0 } ]

            #### Only if we haven't processed this node already do we continue
            if node_id not in followed_peaks:

                #### Extend as far as possible starting from here
                while not done:

                    new_active_sequences = []
                    if verbose == 2:
                        print("==================================================================")
                        print(f"n active_sequences = {len(active_sequences)}")
                        #print(f"active_sequences = ", active_sequences)
                        sorted_active_sequences = sorted(active_sequences,key=lambda x: x['score_per_100mz'], reverse=True)
                        for active_sequence in sorted_active_sequences:
                            print(active_sequence)

                    #### Loop over all active sequences and try to extend them
                    for active_sequence in active_sequences:

                        current_node_id = active_sequence['current_node_id']
                        current_node = spectrum.network['nodes'][current_node_id]
                        current_node_mz = current_node['mz']
                        if verbose == 2:
                            print(f"#######################################################################################################")
                            print(f"######################## Currently examining node {current_node_id} with active sequence {active_sequence}")
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

                            #### Get handles to the edge and the next node that we're extending to
                            edge = spectrum.network['edges'][edge_id]
                            next_node_id = edge['node2_id']
                            if next_node_id == current_node_id:
                                next_node_id = edge['node1_id']
                            next_node = spectrum.network['nodes'][next_node_id]

                            #### Debugging
                            if verbose == 2:
                                print(f"############# Currently examining edge {edge_id}, next_node_id={next_node_id}")
                                print("edge = " + json.dumps(edge, indent=2, sort_keys=True))

                            #### Clone the current active_sequence so that we can build on it
                            new_active_sequence = copy.deepcopy(active_sequence)

                            #### Special logic for when we're just starting out
                            if len(new_active_sequence['sequence']) == 0:

                                #### If we're starting from the bottom virtual node (i.e. at 0 m/z)
                                if current_node_id == 'bottom':
                                    new_active_sequence['terminus_type'] = edge['terminus_type']
                                    new_active_sequence['sequence'].append( edge['terminus_name'] )
                                    new_active_sequence['sequence'].append( edge['label'] )
                                    new_active_sequence['current_node_id'] = next_node_id
                                    new_active_sequence['score'] += edge['score'] + 0.4
                                    new_active_sequence['score_per_100mz'] += new_active_sequence['score'] / next_node['mz'] * 100
                                    new_active_sequence['mass'] += edge['mass']
                                    new_active_sequences.append(new_active_sequence)
                                    if verbose == 2: print('++ Created new_active_sequence = ',new_active_sequence)


                                    #### Create a checkpoint data component
                                    checkpoint_terminus_type = new_active_sequence['terminus_type']
                                    checkpoint_sequence_string = ''.join(new_active_sequence['sequence'])
                                    checkpoint = {
                                        'sequence_string': checkpoint_sequence_string,
                                        'mz': next_node['mz'],
                                        'mass': new_active_sequence['mass'],
                                        'residues': new_active_sequence['sequence'],
                                        'score': new_active_sequence['score'],
                                    }

                                    if checkpoint_sequence_string in sequencing_checkpoints[checkpoint_terminus_type]['dict']:
                                        if verbose == 2: print(f"== checkpoint {checkpoint_sequence_string} already exists")
                                    else:
                                        sequencing_checkpoints[checkpoint_terminus_type]['dict'][checkpoint_sequence_string] = checkpoint
                                        sequencing_checkpoints[checkpoint_terminus_type]['list'].append(checkpoint)
                                        if verbose == 2: print('++ Created new checkpoint = ',checkpoint)

                                    #### If this is a double, then we need to split it as well
                                    if len(edge['residues']) > 1:
                                        for split_residue in edge['residues']:
                                            split_active_sequence = copy.deepcopy(active_sequence)
                                            
                                            split_active_sequence['sequence'].append( edge['terminus_name'] )
                                            split_active_sequence['sequence'].append( split_residue )
                                            split_active_sequence['score'] += ( edge['score'] + 0.4 ) / 2
                                            split_active_sequence['mass'] += self.mass_reference.aa_masses[split_residue]

                                            checkpoint_sequence_string = ''.join(split_active_sequence['sequence'])
                                            checkpoint = {
                                                'sequence_string': checkpoint_sequence_string,
                                                'mz': '{:.4f}'.format(split_active_sequence['mass']),
                                                'mass': split_active_sequence['mass'],
                                                'residues': split_active_sequence['sequence'],
                                                'score': split_active_sequence['score'],
                                            }

                                            if checkpoint_sequence_string in sequencing_checkpoints[checkpoint_terminus_type]['dict']:
                                                if verbose == 2: print(f"## checkpoint {checkpoint_sequence_string} already exists")
                                            else:
                                                sequencing_checkpoints[checkpoint_terminus_type]['dict'][checkpoint_sequence_string] = checkpoint
                                                sequencing_checkpoints[checkpoint_terminus_type]['list'].append(checkpoint)
                                                if verbose == 2: print('@@ Created new checkpoint = ',checkpoint)


                                #### Otherwise we are starting somewhere in the middle
                                else:
                                    new_active_sequence['terminus_type'] = 'internal'
                                    new_active_sequence['sequence'].append( '(+' + '{:.4f}'.format(current_node_mz) + ')' )
                                    new_active_sequence['sequence'].append( edge['label'] )
                                    new_active_sequence['score'] += edge['score'] + 0.4
                                    new_active_sequence['score_per_100mz'] += new_active_sequence['score'] / next_node['mz'] * 100
                                    new_active_sequence['mass'] += edge['mass']
                                    new_active_sequence['current_node_id'] = next_node_id
                                    new_active_sequences.append(new_active_sequence)
                                    if verbose == 2: print('    ++',new_active_sequence)

                                    #### Need some sort of checkpoint here?? Maybe not because we're not going to get a full sequence here anyway?? FIXME


                            #### Continuing the sequence
                            else:

                                #### If we have reached the top, then we are done
                                if next_node_id == 'top':
                                    #print(json.dumps(current_node, indent=2, sort_keys=True))
                                    #print(json.dumps(edge, indent=2, sort_keys=True))
                                    #print(json.dumps(next_node, indent=2, sort_keys=True))
                                    #exit()
                                    new_active_sequence['sequence'].append( edge['label'] )
                                    new_active_sequence['sequence'].append( edge['terminus_name'] )
                                    new_active_sequence['score'] += edge['score'] + 0.4
                                    new_active_sequence['mass'] += edge['mass']
                                    finished_sequences.append(new_active_sequence)
                                    if verbose == 2: print('    ====',new_active_sequence)

                                #### If we're arrived at a place where we can go no further, finish off here
                                elif len(next_node['upward_edges']) == 0:
                                    new_active_sequence['sequence'].append( '(+' + '{:.4f}'.format( precursor_mw + masses['proton'] - next_node['mz'] ) + ')' )
                                    new_active_sequence['score_per_100mz'] += new_active_sequence['score'] / next_node['mz'] * 100

                                    finished_sequences.append(new_active_sequence)
                                    if verbose == 2:
                                        print(f"  xxxx There are no upward edges, so terminate here at {current_node_mz}")
                                        print('        xxxx',active_sequence)

                                #### We can continue so just continue the chain
                                else:
                                    new_active_sequence['sequence'].append( edge['label'] )
                                    new_active_sequence['score'] += edge['score'] + 0.4
                                    new_active_sequence['score_per_100mz'] += new_active_sequence['score'] / next_node['mz'] * 100
                                    new_active_sequence['mass'] += edge['mass']
                                    new_active_sequence['current_node_id'] = next_node_id
                                    new_active_sequences.append(new_active_sequence)
                                    if verbose == 2: print('    ++',new_active_sequence)


                                #### Create a checkpoint data component
                                if next_node_id == 'top' or len(next_node['upward_edges']) != 0:
                                    checkpoint_terminus_type = new_active_sequence['terminus_type']
                                    checkpoint_sequence_string = ''.join(new_active_sequence['sequence'])

                                    checkpoint = {
                                        'sequence_string': checkpoint_sequence_string,
                                        'mz': next_node['mz'],
                                        'mass': new_active_sequence['mass'],
                                        'residues': new_active_sequence['sequence'],
                                        'score': new_active_sequence['score'],
                                    }
                                    if checkpoint_sequence_string in sequencing_checkpoints[checkpoint_terminus_type]['dict']:
                                        if verbose == 2: print(f"== checkpoint {checkpoint_sequence_string} already exists")
                                        pass
                                    else:
                                        sequencing_checkpoints[checkpoint_terminus_type]['dict'][checkpoint_sequence_string] = checkpoint
                                        sequencing_checkpoints[checkpoint_terminus_type]['list'].append(checkpoint)
                                        if verbose == 2: print('** Created new checkpoint = ',checkpoint)


                                    #### If this is a double, then we need to split it as well
                                    if len(edge['residues']) > 1:
                                        for split_residue in edge['residues']:
                                            split_active_sequence = copy.deepcopy(active_sequence)
                                            
                                            split_active_sequence['sequence'].append( split_residue )
                                            split_active_sequence['score'] += ( edge['score'] + 0.4 ) / 2
                                            split_active_sequence['mass'] += self.mass_reference.aa_masses[split_residue]

                                            checkpoint_sequence_string = ''.join(split_active_sequence['sequence'])
                                            checkpoint = {
                                                'sequence_string': checkpoint_sequence_string,
                                                'mz': '{:.4f}'.format(split_active_sequence['mass']),
                                                'mass': split_active_sequence['mass'],
                                                'residues': split_active_sequence['sequence'],
                                                'score': split_active_sequence['score'],
                                            }

                                            if checkpoint_sequence_string in sequencing_checkpoints[checkpoint_terminus_type]['dict']:
                                                if verbose == 2: print(f"!! checkpoint {checkpoint_sequence_string} already exists")
                                            else:
                                                sequencing_checkpoints[checkpoint_terminus_type]['dict'][checkpoint_sequence_string] = checkpoint
                                                sequencing_checkpoints[checkpoint_terminus_type]['list'].append(checkpoint)
                                                if verbose == 2: print('$$ Created new checkpoint = ',checkpoint)


                            followed_peaks[next_node_id] = 1

                    active_sequences = new_active_sequences
                    if verbose: print(f"Number of active sequences: {len(active_sequences)}, finished sequences: {len(finished_sequences)}")
                    #for active_sequence in active_sequences:
                    #    print(active_sequence)
                    #print('**',active_sequences)

                    if len(active_sequences) == 0:
                        done = True

            #### For debugging, stop after the first node
            #break


        #### Sort all the checkpoints
        sequencing_checkpoints['nterm']['list'] = sorted(sequencing_checkpoints['nterm']['list'],key=lambda x: x['score'], reverse=True)
        sequencing_checkpoints['cterm']['list'] = sorted(sequencing_checkpoints['cterm']['list'],key=lambda x: x['score'], reverse=True)

        print("***************************************************************************")
        precursor_neutral_mass = precursor_mz * precursor_charge - masses['proton'] * precursor_charge
        print(f"precursor_mz={precursor_mz}, precursor_charge={precursor_charge}, precursor_neutral_mass={precursor_neutral_mass}")

        nterm_mass_index = {}
        nterm_sequence_index = {}
        merged_finished_sequences = []
        merged_finished_sequence_scores = {}

        print(f"#### Storing {len(sequencing_checkpoints['nterm']['list'])} nterm sequencing checkpoints in an index:")
        for checkpoint in sequencing_checkpoints['nterm']['list']:
            mass = checkpoint['mass']
            int_mass = int(mass)
            if int_mass not in nterm_mass_index:
                nterm_mass_index[int_mass] = []
            nterm_mass_index[int_mass].append(checkpoint)
            sequence_string = checkpoint['sequence_string']
            nterm_sequence_index[sequence_string] = checkpoint

        #### Try storing 0
        sequence_string = ''
        checkpoint = {'sequence_string': sequence_string, 'mz': 0.0, 'mass': 0.0, 'residues': [sequence_string], 'score': 0}
        nterm_mass_index[0] = []
        nterm_mass_index[0].append(checkpoint)
        nterm_sequence_index[sequence_string] = checkpoint

        noisy = False
        print(f"#### Top 50 of {len(sequencing_checkpoints['cterm']['list'])} cterm sequencing checkpoints:")
        counter = 0
        for checkpoint in sequencing_checkpoints['cterm']['list']:
            if False:
                print("-------------")
                print(f"{checkpoint['score']}\t{checkpoint['mass']}\t{checkpoint['sequence_string']}")
                print(checkpoint)
            neutral_mass = checkpoint['mass'] - masses['proton']
            #if checkpoint['residues'][-1].startswith('-') or checkpoint['residues'][-1].endswith('-'):
            if False:
                neutral_mass -= masses['proton']
                theoretical_mz = ( neutral_mass + masses['proton'] * precursor_charge ) / precursor_charge
                delta_mz = theoretical_mz - precursor_mz
                delta_mz_ppm = delta_mz / precursor_mz * 1e6
                #print(f"Complete sequencing theoretical_mz={theoretical_mz}, delta_mz={delta_mz}, delta_mz_ppm={delta_mz_ppm}")
            else:
                complement_mass = precursor_neutral_mass - checkpoint['mass'] + masses['proton'] * precursor_charge
                if noisy:
                    print("-------------")
                    print(f"{checkpoint['score']}\t{checkpoint['mass']}\t{checkpoint['sequence_string']}")
                    print(checkpoint)
                    print(f"Searching for complement mass {complement_mass}")
                int_complement_mass = int(complement_mass)
                if int_complement_mass in nterm_mass_index:
                    for matching_nterm in nterm_mass_index[int_complement_mass]:
                        if noisy:
                            print(matching_nterm)
                        delta_mz = complement_mass - matching_nterm['mass']
                        delta_mz_ppm = delta_mz / precursor_mz * 1e6
                        if abs(delta_mz_ppm) < 50:
                            if noisy:
                                print(f"Found a promising lead with delta_mz_ppm={delta_mz_ppm}. Try to extend")
                            first_matching_nterm = matching_nterm
                            deepest_matching_nterm = matching_nterm
                            residues_to_extend = reversed(checkpoint['residues'])
                            search_sequence = matching_nterm['sequence_string']

                            #### Search to see how far we can match the residues from the c terminus toward the n terminus
                            for search_residue in residues_to_extend:
                                search_sequence = search_sequence + search_residue
                                if noisy:
                                    print(f"Search for {search_sequence}")
                                #### If we find a match, then note this as the deepest yet              # FIXME This does not take into account AA pairs! It needs to!
                                if search_sequence in nterm_sequence_index:
                                    if noisy:
                                        print(f"  Found {nterm_sequence_index[search_sequence]}")
                                    deepest_matching_nterm = nterm_sequence_index[search_sequence]
                                #### Else, if we can't find it, then there's no need to continue any further
                                else:
                                    break

                            #### With what we've found, create a final sequence
                            missing_sequence = list(reversed(first_matching_nterm['residues']))
                            final_sequence = copy.deepcopy(checkpoint)
                            final_sequence['residues'].extend(missing_sequence)
                            final_sequence['sequence_string'] += ''.join(missing_sequence)
                            final_sequence['score'] += deepest_matching_nterm['score']
                            final_sequence['delta_mz_ppm'] = delta_mz_ppm
                            if final_sequence['sequence_string'] not in merged_finished_sequence_scores or final_sequence['score'] > merged_finished_sequence_scores[final_sequence['sequence_string']]:
                                merged_finished_sequences.append(final_sequence)
                                merged_finished_sequence_scores[final_sequence['sequence_string']] = final_sequence['score']
                #counter += 1
                #if counter > 50:
                #    break


        print("*****************************************")
        print(f"Displaying {len(merged_finished_sequences)} merged_finished_sequences")
        sorted_merged_finished_sequences = sorted(merged_finished_sequences,key=lambda x: x['score'], reverse=True)
        counter = 0
        first_score = None
        for merged_finished_sequence in sorted_merged_finished_sequences:
            merged_finished_sequence['residues'] = list(reversed(merged_finished_sequence['residues']))
            merged_finished_sequence['sequence_string'] = ''.join(merged_finished_sequence['residues'])
            if first_score is None:
                first_score = merged_finished_sequence['score']
            if ( first_score - merged_finished_sequence['score'] ) / first_score > 0.2:
                break
            print(f"{merged_finished_sequence['score']}\t{merged_finished_sequence['delta_mz_ppm']}\t{merged_finished_sequence['sequence_string']}")
            counter += 1
            if counter > 20:
                break


        return


        #### Sort all results
        sorted_finished_sequences = sorted(finished_sequences,key=lambda x: x['score_per_100mz'], reverse=True)

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
            if counter < 0:
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
            print(f"{sequence['score_per_100mz']}\t{sequence['score']}\t{sequence['sequence']}\t{''.join(sequence['sequence'])}")
            counter += 1
            if counter > 50:
                break

        print("#### Best nterm sequences:")
        counter = 0
        for sequence in nterm_sequences:
            print(f"{sequence['score_per_100mz']}\t{sequence['score']}\t{sequence['sequence']}\t{''.join(sequence['sequence'])}")
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

                #### Set a hop penalty to penalize double hops by a factor of 2
                hop_penalty = 1.0

                #### Set the residue mass options depending on hop
                if hop == 'single_residues':
                    if verbose: print(f"Looking to extend {sequence['working_mass']} {direction} with single residues")
                    aa_masses = residue_masses
                else:
                    if verbose: print(f"Looking to extend {sequence['working_mass']} {direction} with residue pairs")
                    #aa_masses = self.mass_reference.nr_double_aa_masses
                    aa_masses = self.mass_reference.double_aa_masses
                    hop_penalty = 0.5

                #### Loop over all the potential residue masses
                for residue,residue_mass in aa_masses.items():

                    #### Skip I because it is the same as L
                    if residue == 'I':
                        continue
                    #### Skip a pair if one of the residues in the pair matches a successful single hop
                    residues = [ residue ]
                    if hop == 'double_residues':
                        residues = aa_masses[residue]['residues']
                        skip_residue = False
                        for item in residues:
                            if item in successful_single_hops:
                                #print(f"Skipping {residue} because there is an {item} already in successful_single_hops")
                                skip_residue = True
                        if skip_residue:
                            continue
                        residue_mass = aa_masses[residue]['mass']

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
                            'score': match['score'] * hop_penalty,
                            'mass': residue_mass + terminus_mass,
                            'delta_ppm': match['delta_ppm'],
                            'type': 'residue',
                            'label': residue,
                            'residues': residues,
                        }
                        #### If this is a terminal edge, add that information
                        if potential_terminus['terminus_type'] is not None:
                            edge['terminus_type'] = potential_terminus['terminus_type']
                            edge['terminus_name'] = potential_terminus['terminus_name']
                            edge['score'] += 1.0

                        if edge_id not in spectrum.network['edges']:
                            spectrum.network['edges'][edge_id] = edge

                        spectrum.network['nodes'][node1_id][node1_direction+'_edges'][edge_id] = match['score'] * hop_penalty
                        spectrum.network['nodes'][node2_id][node2_direction+'_edges'][edge_id] = match['score'] * hop_penalty

                        #### If this is a single hop, record that it was a successful so we don't try it with a double hop
                        if hop == 'single_residues':
                            successful_single_hops[residue] = 1



####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Class representing a peptidoform')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('--usi', action='store', default=None, help='USI to process')
    argparser.add_argument('--precursor_mz', action='store', help='The precursor m/z of the ion' )
    argparser.add_argument('--labels', action='store', default=None, help='Quantitative labels to consider (e.g. TMT6plex)')
    argparser.add_argument('--input_json_filename', action='store', default=None, type=float, help='Filename of an input json file')
    argparser.add_argument('--tolerance', action='store', default=None, type=float, help='Tolerance in ppm for annotation')
    argparser.add_argument('--examine', action='count', help='If set, examine the spectrum to see what can be learned' )
    argparser.add_argument('--annotate', action='count', help='If set, annotate the USI spectrum' )
    argparser.add_argument('--denovo', action='count', help='If set, interpret the USI spectrum de novo' )
    argparser.add_argument('--create_network', action='count', help='If set, clean the spectrum and create a peak network' )
    argparser.add_argument('--score', action='count', help='If set, score the spectrum with the supplied peptidoform' )
    argparser.add_argument('--show_all_annotations', action='count', help='If set, show all the potential annotations, not just the final one' )
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

    #### If the user specified labels
    labels = []
    if params.labels is not None and params.labels != '':
        labels = params.labels.split(',')

    #### Need to do this as apparently the peptidoform that comes back from usi is a dict, not an object?
    peptidoform = None
    if peptidoform_string is not None and peptidoform_string != '':
        peptidoform = ProformaPeptidoform(peptidoform_string)

    #### If the USI specifies some quanititative labels, then add those to the list of things to consider
    if peptidoform_string is not None and peptidoform_string != '':
        if len(labels) == 0:
            possible_labels = [ 'TMT6plex', 'TMTpro', 'iTRAQ4plex', 'iTRAQ8plex' ]
            for possible_label in possible_labels:
                if f"[{possible_label}]" in peptidoform_string:
                    labels.append(possible_label)


    if params.create_network:
        examiner = SpectrumExaminer()
        sequencer = SpectrumSequencer()
        # First, clean the spectrum of stuff that is not useful for de novo sequencing
        examiner.identify_isotopes(spectrum)
        examiner.delete_isotopes(spectrum)
        examiner.identify_precursors(spectrum)
        examiner.delete_precursors(spectrum)
        examiner.identify_low_mass_ions(spectrum)
        examiner.analyze_mass_defects(spectrum, delete_outliers=True)
        #annotator.remove_low_mass_ions(spectrum)
        #annotator.remove_low_intensity_peaks(spectrum)
        # After cleaning, recompute the metrics on what's left and re-index
        spectrum.compute_spectrum_metrics()

        # Set some basic input parameters
        precursor_mz = spectrum.analytes['1']['precursor_mz']
        if params.precursor_mz is not None:
            precursor_mz = float(params.precursor_mz)
        sequencing_parameters = {
            'fragmentation_type': 'HCD',
            'tolerance': 20.0,
            'precursor_mz': precursor_mz,
            'precursor_charge': spectrum.analytes['1']['charge state'],
            'labels': labels
        }
        sequencer.create_peak_network(spectrum, sequencing_parameters=sequencing_parameters)
        print(json.dumps(spectrum.network, indent=2, sort_keys=True))
        return


    if params.denovo:
        examiner = SpectrumExaminer()
        sequencer = SpectrumSequencer()
        # First, clean the spectrum of stuff that is not useful for de novo sequencing
        examiner.identify_isotopes(spectrum)
        examiner.delete_isotopes(spectrum)
        examiner.identify_precursors(spectrum)
        examiner.delete_precursors(spectrum)
        examiner.identify_low_mass_ions(spectrum)
        examiner.analyze_mass_defects(spectrum, delete_outliers=True)
        #annotator.remove_low_mass_ions(spectrum)
        #annotator.remove_low_intensity_peaks(spectrum)
        # After cleaning, recompute the metrics on what's left and re-index
        spectrum.compute_spectrum_metrics()

        # Set some basic input parameters
        precursor_mz = spectrum.analytes['1']['precursor_mz']
        if params.precursor_mz is not None:
            precursor_mz = float(params.precursor_mz)
        sequencing_parameters = {
            'fragmentation_type': 'HCD',
            'tolerance': 20.0,
            'precursor_mz': precursor_mz,
            'precursor_charge': spectrum.analytes['1']['charge state'],
            'labels': labels
        }
        sequencer.interpret_de_novo(spectrum, sequencing_parameters=sequencing_parameters)


#### For command line usage
if __name__ == "__main__": main()
