#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import timeit
import re
import json
import numpy
import pickle
import gzip
from lxml import etree
import multiprocessing
import copy

from metadata_handler import MetadataHandler

#### Import technical modules and pyteomics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from pyteomics import mzml, auxiliary

masses = {
    'proton': 1.0072765,
    'hydrogen': 1.007825035,
    'oxygen': 15.99491463,
    'carbon': 12.0000000,
    'nitrogen': 14.0030740,
    'phosphorus': 30.973762,
    'sulfur': 31.9720707,
}

ion_series_attr = {
    'a': { 'terminus_type': 'nterm', 'mass': -1 * masses['carbon'] - masses['oxygen'], 'complement': 'a' },
    'b': { 'terminus_type': 'nterm', 'mass': 0.0, 'complement': 'y' },
    'c': { 'terminus_type': 'nterm', 'mass': 3 *  masses['hydrogen'] + masses['nitrogen'], 'complement': 'z' },
    'x': { 'terminus_type': 'cterm', 'mass': masses['carbon'] + 2 * masses['oxygen'], 'complement': 'a' },
    'y': { 'terminus_type': 'cterm', 'mass': 2 * masses['hydrogen'] + masses['oxygen'], 'complement': 'b'},
    'z': { 'terminus_type': 'cterm', 'mass': masses['oxygen'] - masses['nitrogen'], 'complement': 'c' }
}


####################################################################################################
#### SpectrumProcessor class
class SpectrumProcessor:

    ####################################################################################################
    #### Constructor
    def __init__(self, mzml_file, verbose=0):

        # Store input file
        self.mzml_file = mzml_file

        # Set verbosity
        if verbose is None: verbose = 0
        self.verbose = verbose

        # Store the provided metadata or create a template to work in
        metadata_handler = MetadataHandler()
        self.metadata = metadata_handler.create_template()
        self.metadata['files'][mzml_file] = {}

        # Create a place to store our composite spectra for analysis
        self.composite = {}

        self.aa_masses = None
        self.double_aa_masses = None
        self.nterm_modifications = None
        self.nterm_aa_masses = None
        self.nterm_double_aa_masses = None
        self.reporter_ions = None

        self.prepare_mass_tables()
        self.reference_peak_attrs = {}

        self.stats = {}
        self.run_mode = 'unspecified'

        self.user_parameters = {
            'first_pass_tolerance': 20,
        }


    ####################################################################################################
    def prepare_mass_tables(self):

        self.aa_masses = {
            'G': 57.0214636,
            'A': 71.0371136,
            'S': 87.0320282,
            'P': 97.0527636,
            'V': 99.0684136,
            'T': 101.0476782,
            'C': 103.0091854,
            'L': 113.0840636,
            #'I': 113.0840636,
            #'X': 113.0840636,
            'N': 114.0429272,
            #'O': 114.0793126,
            #'B': 114.534935,
            'D': 115.0269428,
            'Q': 128.0585772,
            'K': 128.0949626,
            #'Z': 128.550585,
            'E': 129.0425928,
            'M': 131.0404854,
            'H': 137.0589116,
            'F': 147.0684136,
            'R': 156.1011106,
            'Y': 163.0633282,
            'W': 186.0793126
        }

        # Have a version with no modifications
        nomod_aa_masses = copy.deepcopy(self.aa_masses)

        # Should make a catalog of equivalencies and not search them at the beginning, but add them in later during final evaluation
        # Would be more efficient to compute these rather than curate by hand
        # Q[Deamidated] = E
        # E[Glu->pyro-Glu] = Q[Gln->pyro-Glu]

        # Add some modifications
        modifications = {
            'Phospho': { 'delta_mass': 79.966331, 'residues': [ 'S', 'T', 'Y' ], 'neutral_losses': [ 97.976896 ] },
            #'Phospho': { 'delta_mass': 79.966331, 'residues': [ 'Y' ], 'neutral_losses': [ ] },
            'Pyrophospho': { 'delta_mass': 159.932662, 'residues': [ 'S', 'T', 'Y' ], 'neutral_losses': [ 177.943227, 79.966331 ] },
            'Carbamidomethyl': { 'delta_mass': 57.021464, 'residues': [ 'C' ] },
            #'Beta-methythiolation': { 'delta_mass': 45.987721, 'residues': [ 'C' ] },
            'Oxidation': { 'delta_mass': 15.994915, 'residues': [ 'M', 'W', 'P' ] },
            'Deamidated': { 'delta_mass': 0.984016, 'residues': [ 'N' ] }, # Q[Deamidated] is the same as E, so don't bother here
            #'Aceylation': { 'delta_mass': 42.010565, 'residues': [ 'K' ] },
            #'Methylthio': { 'delta_mass': 45.987721, 'residues': [ 'K', 'D', 'N', 'C' ] },
            #'Methylation': { 'delta_mass': 14.015650, 'residues': [ 'K', 'R' ] },
            #'Ubiquitination': { 'delta_mass': 114.042927, 'residues': [ 'K' ] },
            #'Carbamylation': { 'delta_mass': 43.005814, 'residues': [ 'K' ] },
            #'TMT6': { 'delta_mass': 229.162932, 'residues': [ 'K' ] },
            #'iTRAQ4': { 'delta_mass': 144.102063, 'residues': [ 'K' ] },
        }
        for modification in modifications:
            for residue in modifications[modification]['residues']:
                self.aa_masses[f"{residue}[{modification}]"] = self.aa_masses[residue] + modifications[modification]['delta_mass']

        # Add n-terminal modifications
        self.nterm_aa_modifications = {
            'Glu->pyro-Glu': { 'delta_mass': -18.010565, 'residues': [ 'E' ] },
            #'Gln->pyro-Glu': { 'delta_mass': -17.026549, 'residues': [ 'Q' ] } # only need to do one, because they end up the same mass
        }
        self.nterm_aa_masses = copy.deepcopy(self.aa_masses)
        nterm_residues = {}
        for modification in self.nterm_aa_modifications:
            for residue in self.nterm_aa_modifications[modification]['residues']:
                self.nterm_aa_masses[f"{residue}[{modification}]"] = self.aa_masses[residue] + self.nterm_aa_modifications[modification]['delta_mass']
                nterm_residues[f"{residue}[{modification}]"] = 1

        # Create all the amino acid pairs
        self.double_aa_masses = {}
        nomod_double_aa_masses = {}
        for aa1 in self.aa_masses:
            for aa2 in self.aa_masses:
                if aa2 + aa1 in self.double_aa_masses:
                    continue
                self.double_aa_masses[aa1 + aa2] = self.aa_masses[aa1] + self.aa_masses[aa2]
                if '[' not in aa1 + aa2:
                    nomod_double_aa_masses[aa1 + aa2] = self.aa_masses[aa1] + self.aa_masses[aa2]

        # Create all the n-terminal amino acid pairs
        self.nterm_double_aa_masses = {}
        for aa1 in self.nterm_aa_masses:
            for aa2 in self.nterm_aa_masses:
                if aa2 in nterm_residues:
                    continue
                if aa2 + aa1 in self.nterm_double_aa_masses:
                    continue
                self.nterm_double_aa_masses[aa1 + aa2] = self.nterm_aa_masses[aa1] + self.nterm_aa_masses[aa2]


        # Define some potential neutral losses
        # Note: it seems that Y[Phospho] does not lose -98! Could maybe use that information somewhere useful!
        self.neutral_losses = {
            'water': { 'formula': 'H2O', 'delta_mass': masses['hydrogen'] * 2 + masses['oxygen'], 'residues': [ 'S', 'T', 'E', 'D' ] },
            'ammonia': { 'formula': 'NH3', 'delta_mass': masses['nitrogen'] + masses['hydrogen'] * 3, 'residues': [ 'R', 'K', 'N', 'Q' ] },
            'carbon monoxide': { 'formula': 'CO', 'delta_mass': masses['carbon'] + masses['oxygen'], 'residues': [ 'b-ion' ] },
            'phosphoric acid': { 'formula': 'H3PO4', 'delta_mass': masses['hydrogen'] * 3 + masses['phosphorus'] + masses['oxygen'] * 4, 
                'residues': [ 'S[Phospho]', 'T[Phospho]' ] },
            'metaphosphoric acid': { 'formula': 'HPO3', 'delta_mass': masses['hydrogen'] * 1 + masses['phosphorus'] + masses['oxygen'] * 3, 'residues': [ 'Y[Phospho]' ] },
            'methanesulfenic acid': { 'formula': 'CH4OS', 'delta_mass': masses['carbon'] * 1 + masses['hydrogen'] * 4 + masses['oxygen'] * 1 + masses['sulfur'] * 1,
                'residues': [ 'M[Oxidation]' ] },
        }

        self.immonium_ions = {}
        for residue,residue_mass in self.aa_masses.items():
            immonium_name = f"Im{residue}A"
            immonium_mass = residue_mass - ( masses['carbon'] + masses['oxygen'] ) + masses['proton']
            self.immonium_ions[immonium_name] = immonium_mass

        self.low_mass_ions = copy.deepcopy(self.immonium_ions)
        for residue in [ 'R', 'K' ]:
            ion_name = f"y1({residue})"
            ion_mass = self.aa_masses[residue] + 2 * masses['hydrogen'] + masses['oxygen'] + masses['proton']
            self.low_mass_ions[ion_name] = ion_mass

            ion_name = f"y1({residue})-H20"
            ion_mass = self.aa_masses[residue] + masses['proton']
            self.low_mass_ions[ion_name] = ion_mass

        # Add to the low mass ions, all the b2 and a2 ions, which are very common and make for a nice calibrator
        # But exclude a few that seem to generate some side noise
        double_aas_to_skip = { 'QH', 'YY', 'LW' }
        for residue_pair,mass in nomod_double_aa_masses.items():
            if residue_pair in double_aas_to_skip or 'M' in residue_pair:
                continue
            ion_name = f"b2({residue_pair})"
            ion_mass = mass + masses['proton']
            self.low_mass_ions[ion_name] = ion_mass

            ion_name = f"a2({residue_pair})"
            ion_mass = mass - masses['carbon'] - masses['oxygen'] + masses['proton']
            self.low_mass_ions[ion_name] = ion_mass

        # Define reporter ions to look for
        self.reporter_ions = {
            'TMT6/10_126': { 'type': 'TMT', 'mz': 126.127726 },     # clean in TMT6
            'TMT6/10_127N': { 'type': 'TMT', 'mz': 127.124761 },    # clean in TMT6
            'TMT10_127C': { 'type': 'TMT', 'mz': 127.131081 },      # Need to first find and calibrate on TMT6 and only then do a narrow look here
            'TMT10_128N': { 'type': 'TMT', 'mz': 128.128116 },
            'TMT6/10_128C': { 'type': 'TMT', 'mz': 128.134436 },   # clean in TMT6
            'TMT6/10_129N': { 'type': 'TMT', 'mz': 129.131471 },   # clean in TMT6
            'TMT10_129C': { 'type': 'TMT', 'mz': 129.137790 },
            'TMT10_130N': { 'type': 'TMT', 'mz': 130.134825 },
            'TMT6/10_130C': { 'type': 'TMT', 'mz': 130.141145 },
            'TMT6/10_131': { 'type': 'TMT', 'mz': 131.138180 },
            'TMT11_131C': { 'type': 'TMT', 'mz': 131.138180 },
            'TMT_nterm': { 'type': 'TMT', 'mz': 229.162932 + masses['proton'] },

            # My old numbers from somewhere
            'iTRAQ4_114': { 'type': 'iTRAQ', 'mz': 114.11068 },
            'iTRAQ4_115': { 'type': 'iTRAQ', 'mz': 115.107715 },
            'iTRAQ4_116': { 'type': 'iTRAQ', 'mz': 116.111069 },
            'iTRAQ4_117': { 'type': 'iTRAQ', 'mz': 117.114424 },
            'iTRAQ4_nterm': { 'type': 'iTRAQ4', 'mz': 145.109 },

            # Jimmy's numbers from https://proteomicsresource.washington.edu/protocols03/isotopic_labeling.php
            #'iTRAQ4_114': { 'type': 'iTRAQ4', 'mz': 114.1112 },
            #'iTRAQ4_115': { 'type': 'iTRAQ4', 'mz': 115.1083 },
            #'iTRAQ4_116': { 'type': 'iTRAQ4', 'mz': 116.1116 },
            #'iTRAQ4_117': { 'type': 'iTRAQ4', 'mz': 117.1150 },
            #'iTRAQ4_nterm': { 'type': 'iTRAQ4', 'mz': (144.105918 + 144.099599 + 144.102063 + 144.102063)/4 + masses['proton'] },

            # My old numbers from somewhere
            'iTRAQ8_118': { 'type': 'iTRAQ8', 'mz': 118.111459 },   #confounder?
            'iTRAQ8_119': { 'type': 'iTRAQ8', 'mz': 119.114814 },
            'iTRAQ8_121': { 'type': 'iTRAQ8', 'mz': 121.121524 },
            'iTRAQ8_113': { 'type': 'iTRAQ8', 'mz': 113.107325 },   #confounder??

            # Numbers from Jimmy: 113.1078	114.1112	115.1082	116.1116	117.1149	118.1120	119.1153	121.1220
        }



    ####################################################################################################
    #### Read spectra
    def read_spectra(self, write_fragmentation_type_file=None):

        #### Set up information
        t0 = timeit.default_timer()
        stats = { 
            'n_spectra': 0, 
            'n_ms1_spectra': 0, 
            'n_ms2_spectra': 0,
            'n_HR_HCD_spectra': 0,
            'n_LR_HCD_spectra': 0,
            'n_HR_IT_CID_spectra': 0,
            'n_LR_IT_CID_spectra': 0,
            'n_HR_IT_ETD_spectra': 0,
            'n_LR_IT_ETD_spectra': 0,
            'n_HR_EThcD_spectra': 0,
            'n_LR_EThcD_spectra': 0,
            'n_HR_ETciD_spectra': 0,
            'n_LR_ETciD_spectra': 0,
            'high_accuracy_precursors': 'unknown', 
            'fragmentation_type': 'unknown', 
            'fragmentation_tag': 'unknown',
            'n_phospho_loss_spectra': 0, 
            'n_unfragmented_precursors': 0, 
            'zero_peak_spectra':0,
            'mz_calibration_offset_ppm': 0.0,
            'mz_calibration_sigma_ppm': 2.5,
            'mz_calibration_tolerance_ppm': 7.5,
            }

        self.metadata['files'][self.mzml_file]['spectra_stats'] = stats
        self.stats = stats

        #### Store the fragmentation types in a list for storing in the fragmentation_type_file
        scan_fragmentation_list = []

        #### Short circuit for testing
        if 0:
            table = pd.read_csv('residuals.tsv',sep='\t')
            self.assess_fiducial_peaks(table)
            sys.exit(0)

        #### Show information
        if self.verbose >= 1:
            eprint(f"INFO: Processing mzML file {self.mzml_file}")

        #### If the mzML is gzipped, then open with zlib, else a plain open
        match = re.search(r'\.gz$',self.mzml_file)
        if match:
            infile = gzip.open(self.mzml_file)
        else:
            infile = open(self.mzml_file, 'rb')

        #### Store all the spectra in a list for later potential multiprocessing
        spectra = []

        #### Read spectra from the file
        eprint(f"INFO: Reading spectra.. ", end='', flush=True)
        with mzml.read(infile) as reader:
            try:
                for spectrum in reader:

                    #### Look for a filter string and parse it
                    filter_string = None
                    if 'filter string' in spectrum['scanList']['scan'][0]:
                        filter_string = spectrum['scanList']['scan'][0]['filter string']
                        self.parse_filter_string(filter_string,stats)

                    #### If the ms level is greater than 2, fail
                    if spectrum['ms level'] > 2:
                        self.log_event('ERROR','MSnTooHigh',f"MS level is greater than we can handle at '{spectrum['ms level']}'")
                        break

                    #### If the ms level is 2, then examine it for information
                    if spectrum['ms level'] == 2 and 'm/z array' in spectrum:

                        # Explore the structure of the first spectrum
                        #auxiliary.print_tree(spectrum)
                        #sys.exit(0)

                        # If this is a 0 length spectrum, don't bother with it
                        if len(spectrum['m/z array']) == 0:
                            stats['zero_peak_spectra'] += 1

                        # Otherwise, prepare its structure and add it to the list to process
                        else:
                            charge_state = 1
                            if 'charge state' in spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]:
                                charge_state = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']
                            input_spectrum = {
                                'id': spectrum['id'],
                                'mz_array': spectrum['m/z array'],
                                'intensity_array': spectrum['intensity array'],
                                'precursor_mz': spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'],
                                'precursor_charge': charge_state,
                                'isolation_window': {
                                    'target_mz': spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z'],
                                    'lower_offset': spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window lower offset'],
                                    'upper_offset': spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window upper offset'] } }
                            spectra.append(input_spectrum)

                    #### Update counters and print progress if verbose
                    stats['n_spectra'] += 1
                    if self.verbose >= 1:
                        if stats['n_spectra']/500 == int(stats['n_spectra']/500):
                            eprint(f"{stats['n_spectra']}.. ", end='', flush=True)
                    if stats['n_spectra'] >= 30000:
                        print(f"\nWARNING: Exiting eary after {stats['n_spectra']} spectra for development purposes")
                        break
            except:
                self.log_event('ERROR','MzMLCorrupt',f"Pyteomics mzML reader threw an error reading mzML file! File may be corrupt. Check file '{self.mzml_file}'")

        # Close the file and print final timing information
        infile.close()
        if self.verbose >= 1: eprint("")
        t1 = timeit.default_timer()
        print(f"INFO: Read {stats['n_spectra']} spectra from {self.mzml_file} in {int(t1-t0)} sec ({stats['n_spectra']/(t1-t0)} spectra per sec)")

        #### If there are no spectra that were saved, then we're done with this file
        if stats['n_ms2_spectra'] == 0:
            self.log_event('ERROR','NoMS2Scans',f"This mzML file has no MS2 scans. Check file '{self.mzml_file}'")
            return


        # Perform a first pass to determine mass calibration
        options = { 'perform_mass_calibration': True, 'write_residuals_file': True }
        result = self.process_spectra(spectra, options)

        #### Print first pass timing information
        t2 = timeit.default_timer()
        print(f"INFO: First pass analysis of {stats['n_ms2_spectra']} spectra from {self.mzml_file} complete in {int((t2-t1)*100)/100.0} sec ({stats['n_ms2_spectra']/(t2-t1)} spectra per sec)")


        # Perform a second pass to do the de novo sequencing
        options = { 'denovo_test': True, 'run_parallel': False }
        result = self.process_spectra(spectra, options)

        #### Print first pass timing information
        t3 = timeit.default_timer()
        print(f"INFO: Second pass analysis of {stats['n_ms2_spectra']} spectra from {self.mzml_file} complete in {int((t3-t2)*100)/100.0} sec ({stats['n_ms2_spectra']/(t3-t2)} spectra per sec)")


        # Print out some stats about the run
        for stat_name,stat_value in stats.items():
            print(f"  {stat_name}={stat_value}")

        print(f"INFO: Done")



    ####################################################################################################
    def process_spectra(self, spectra, options):

        # Create a list to store all the generated ppm residuals
        residuals_list = []
        run_parallel = True

        # Set the run_mode based on the options
        if 'perform_mass_calibration' in options and options['perform_mass_calibration'] is True:
            self.run_mode = 'mass_calibration'
        if 'perform_denovo_sequencing' in options and options['perform_denovo_sequencing'] is True:
            self.run_mode = 'perform_denovo_sequencing'
        if 'denovo_test' in options and options['denovo_test'] is True:
            self.run_mode = 'denovo_test'
        if 'run_parallel' in options:
            run_parallel = options['run_parallel']

        # Now process the jobs in parallel
        if run_parallel:
            #n_threads = params.n_threads or multiprocessing.cpu_count()
            n_threads = multiprocessing.cpu_count()
            #n_threads = 1
            eprint(f"INFO: Processing spectra with n_threads={n_threads}", end='', flush=True)
            pool = multiprocessing.Pool(processes=n_threads)
            results = pool.map_async(self.process_spectrum, spectra)
            pool.close()
            pool.join()
            eprint("\nProcessing complete. Merging information.")

            # Coalesce the results
            results = results.get()
            for result in results:
                residuals_list.extend(result['residuals_list'])
                for key,value in result['msrun_attributes'].items():
                    self.stats[key] += value

        # Else if single threaded, just loop over the spectra
        else:
            n_spectra = len(spectra)
            eprint(f"INFO: Processing spectra in single-thread mode...")
            for spectrum in spectra:
                result = self.process_spectrum(spectrum)
                residuals_list.extend(result['residuals_list'])
                for key,value in result['msrun_attributes'].items():
                    self.stats[key] += value

        # Some custom analysis code
        if 'write_residuals_file' in options:
            print(f"INFO: Writing residuals file")
            column_name_list = [ 'type', 'scan_number', 'parent_mz', 'parent_intensity', 'child_mz', 'child_intensity', 'charge', 'delta', 'residual', 'residual_ppm', 'tag' ]
            table = pd.DataFrame(residuals_list, columns = column_name_list)
            table.to_csv('residuals.tsv',sep='\t')

        # Assess all the fiducial peaks
        if 'perform_mass_calibration' in options and options['perform_mass_calibration'] is True:
            self.perform_mass_calibration(table)


    ####################################################################################################
    def process_spectrum(self, input_spectrum):

        # Determine the scan number
        scan_number = -1
        match = re.search(r'scan=(\d+)',input_spectrum['id'])
        if match:
            scan_number = match.group(1)

        # Compute intensity profile: min, max, noise, signal, signal to noise ratio
        # FIXME This would be better and more accurate if we ignore everthing less than 200 and ignore the precursor
        intensity_profile = self.compute_intensity_profile(input_spectrum['intensity_array'])

        # Store information about the spectrum to be returned
        spectrum = {
            'id': input_spectrum['id'],
            'scan_number': scan_number,
            'precursor_mz': input_spectrum['precursor_mz'],
            'precursor_charge': input_spectrum['precursor_charge'],
            'isolation_window': input_spectrum['isolation_window'],
            'n_peaks': len(input_spectrum['mz_array']),
            'intensity_profile': intensity_profile,

            'mz_array': input_spectrum['mz_array'] ,
            'intensity_array': input_spectrum['intensity_array'],

            'charge_array': input_spectrum['mz_array'].astype(int) * 0,
            'isotope_array': input_spectrum['mz_array'].astype(int) * 0,
            'n_isotopes_array': input_spectrum['mz_array'].astype(int) * 0,
            'parent_peak': input_spectrum['mz_array'].astype(int) * 0,
            'isotope_delta': input_spectrum['mz_array'] * 0,
            'keep': input_spectrum['mz_array'].astype(int) * 0 + 1,

            'residuals_list': [],
            'attributes': {
                'has_unfragmented_precursor': 0,
                'n_identified_peptide_low_mass_ions': 0,
                'n_identified_reporter_ions': 0,
                'n_neutral_loss_peaks': 0,
                'n_phospho_loss_peaks': 0,
                'n_oxidation_loss_peaks': 0 },
            'msrun_attributes': {
                'n_unfragmented_precursors': 0,
                'n_phospho_loss_spectra': 0, },
            'peak_attributes': [],
            'peak_interpretation': [],
            'peak_comment': [],
            }


        # Populate some arrays for each peak
        for i_peak in range(spectrum['n_peaks']):
            spectrum['peak_attributes'].append( { 'n_neutral_losses': 0, 'is_neutral_loss': 0 } )
            spectrum['peak_interpretation'].append('')
            spectrum['peak_comment'].append('')

        # Dump out some data for inspection
        #print("\t".join( [ str(scan_number), str(spectrum['n_peaks']), '{:.1f}'.format(intensity_profile[2]), 
        #    '{:.1f}'.format(intensity_profile[3]), '{:.1f}'.format(intensity_profile[4]) ] ) )

        if self.run_mode == 'mass_calibration':
            self.identify_isotopes(spectrum)
            self.index_peaks(spectrum)
            self.identify_neutral_losses(spectrum)
            self.identify_reporter_ions(spectrum)
            self.identify_low_mass_ions(spectrum)

        if self.run_mode == 'perform_denovo_sequencing':
            self.identify_isotopes(spectrum)
            self.index_peaks(spectrum)
            self.identify_neutral_losses(spectrum)
            self.identify_reporter_ions(spectrum)
            self.identify_low_mass_ions(spectrum)
            self.clean(spectrum, deisotope=True, threshold=1, charge_reduce=True)
            self.interpret_de_novo(spectrum)

        if self.run_mode == 'denovo_test':
            interesting_scans = { '2680': 1, '3564': 1, '4814': 1, '4853': 1, '5009': 1, '2818': 1 }
            interesting_scans = { '3564': 1, '4814': 1, '4853': 1, '5009': 1, '2818': 1 }
            interesting_scans = { '3659': 1 }
            interesting_scans = { '3332': 1 }
            interesting_scans = { '3564': 1, '4814': 1, '3796': 1, '4853': 1, '5230': 1, '4336': 1, '5009': 1, '3345': 1, '3659': 1, '2818': 1, '2700': 1, '3651': 1, '3332': 1 }
            #interesting_scans = { '4814': 1 }
            interesting_scans = { '4336': 1 }

            if str(scan_number) in interesting_scans:
                print(f"Spectrum: {input_spectrum['id']}")
                self.index_peaks(spectrum)
                self.identify_isotopes(spectrum)
                self.identify_precursors(spectrum)
                self.identify_neutral_losses(spectrum)
                self.identify_reporter_ions(spectrum)
                self.identify_low_mass_ions(spectrum)
                self.clean(spectrum, deisotope=True, threshold=1, charge_reduce=True)
                self.index_peaks(spectrum)
                self.interpret_de_novo(spectrum)
                self.show_spectrum(spectrum)
                input("Enter to continue...")
                #sys.exit(0)

        # Free memory of stuff we don't need anymore
        attributes_to_remove = [ 'mz_array', 'intensity_array', 'mz_dict', 'processed', 'proc_mz_dict',
            'charge_array', 'isotope_array', 'n_isotopes_array', 'parent_peak', 'isotope_delta' ]
        for attr in attributes_to_remove:
            if attr in spectrum:
                del spectrum[attr]
        #spectrum['residuals_list'] = []

        return spectrum


    ####################################################################################################
    def identify_precursors(self, spectrum):

        # See if there is an unfragmented precursor here
        precursor_mz = spectrum['precursor_mz']
        n_peaks = spectrum['n_peaks']
        tolerance = 20

        match = self.find_closest_ion(spectrum,precursor_mz,tolerance,spectrum_version='input')
        if match:
            if len(spectrum['peak_interpretation'][match['peak_number']]) > 0:
                spectrum['peak_interpretation'][match['peak_number']] += ', '
            spectrum['peak_interpretation'][match['peak_number']] += "p/" + '{:.1f}'.format(match['delta_ppm'])
            spectrum['attributes']['has_unfragmented_precursor'] += 1
            spectrum['msrun_attributes']['n_unfragmented_precursors'] += 1
            spectrum['keep'][match['peak_number']] = 0

        # Otherwise add an artificial one at the end to gather neutral losses more easily
        else:
            spectrum['mz_array'] = np.hstack([spectrum['mz_array'],np.array([precursor_mz])])
            spectrum['intensity_array'] = np.hstack([spectrum['intensity_array'],np.array([spectrum['intensity_profile'][1]])])
            spectrum['peak_interpretation'].append('artificially added precursor')
            spectrum['peak_comment'].append('')
            spectrum['charge_array'] = np.hstack([spectrum['charge_array'],np.array([spectrum['precursor_charge']])])
            spectrum['isotope_array'] = np.hstack([spectrum['isotope_array'],np.array([0])])
            spectrum['n_isotopes_array'] = np.hstack([spectrum['n_isotopes_array'],np.array([0])])
            spectrum['isotope_delta'] = np.hstack([spectrum['isotope_delta'],np.array([0.0])])
            spectrum['isotope_delta'] = np.hstack([spectrum['isotope_delta'],np.array([0.0])])
            spectrum['keep'] = np.hstack([spectrum['keep'],np.array([0])])
            spectrum['peak_attributes'].append( { 'n_neutral_losses': 0, 'is_neutral_loss': 0 } )
            spectrum['n_peaks'] += 1

        # Loop over all peaks looking for peaks in the isolation window and exclude them
        lower_bound = spectrum['isolation_window']['target_mz'] - spectrum['isolation_window']['lower_offset']
        upper_bound = spectrum['isolation_window']['target_mz'] + spectrum['isolation_window']['upper_offset']
        other_precursor_count = 1
        for i_peak in range(n_peaks):
            if spectrum['mz_array'][i_peak] < lower_bound:
                continue
            if spectrum['mz_array'][i_peak] > upper_bound:
                break
            if spectrum['isotope_array'][i_peak] > 0:
                continue
            if len(spectrum['peak_interpretation'][i_peak]) > 0:
                spectrum['peak_interpretation'][i_peak] += ', '
            spectrum['peak_interpretation'][i_peak] += f"0#p{other_precursor_count}/0.00"
            spectrum['keep'][i_peak] = 0
            other_precursor_count += 1


    ####################################################################################################
    def identify_neutral_losses(self, spectrum):

        n_peaks = spectrum['n_peaks']
        mzs_list = []
        tolerance = 10

        # Loop over all peaks looking for neutral losses
        for i_peak in range(n_peaks):

            #### If this peak is already labeled as an isotope, then we can ignore it
            if spectrum['isotope_array'][i_peak]:
                continue

            mz = float(spectrum['mz_array'][i_peak])
            intensity = float(spectrum['intensity_array'][i_peak])
            charge = float( spectrum['charge_array'][i_peak] )

            for neutral_loss,neutral_loss_attrs in self.neutral_losses.items():
                delta_mass = neutral_loss_attrs['delta_mass']
                test_charge = charge or 1
                test_mz = mz - delta_mass / test_charge
                match = self.find_closest_ion(spectrum,test_mz,tolerance,spectrum_version='input')
                if match:
                    spectrum['peak_attributes'][i_peak]['n_neutral_losses'] += 1
                    spectrum['peak_attributes'][match['peak_number']]['is_neutral_loss'] += 1
                    spectrum['parent_peak'][match['peak_number']] = i_peak
                    spectrum['attributes']['n_neutral_loss_peaks'] += 1

                    if len(spectrum['peak_interpretation'][match['peak_number']]) > 0:
                        spectrum['peak_interpretation'][match['peak_number']] += ', '
                    spectrum['peak_interpretation'][match['peak_number']] += neutral_loss + f"({i_peak})/" + '{:.1f}'.format(match['delta_ppm'])
                    if neutral_loss == 'phosphoric acid':
                        spectrum['attributes']['n_phospho_loss_peaks'] += 1

                    # Add this match to the residuals list
                    # [ 'type', 'scan_number', 'parent_mz', 'parent_intensity', 'child_mz', 'child_intensity', 'charge', 'delta', 'residual', 'residual_ppm', 'tag' ]
                    #row = [ neutral_loss_attrs['formula'], spectrum['scan_number'], mz, intensity, match['mz'], spectrum['intensity_array'][match['peak_number']],
                    #    charge, mz - match['mz'], -1, match['delta_ppm'], 'loss' ]
                    #spectrum['residuals_list'].append(row)


        if spectrum['attributes']['n_phospho_loss_peaks'] > 1:
            spectrum['msrun_attributes']['n_phospho_loss_spectra'] += 1


    ####################################################################################################
    def identify_low_mass_ions(self, spectrum):

        tolerance = 20

        for low_mass_ion_name,low_mass_ion_mz in self.low_mass_ions.items():

            match = self.find_closest_ion(spectrum,low_mass_ion_mz,tolerance,spectrum_version='input')
            if match:
                if len(spectrum['peak_interpretation'][match['peak_number']]) > 0:
                    spectrum['peak_interpretation'][match['peak_number']] += ', '
                spectrum['peak_interpretation'][match['peak_number']] += f"0#{low_mass_ion_name}/" + '{:.1f}'.format(match['delta_ppm'])

                spectrum['attributes']['n_identified_peptide_low_mass_ions'] += 1

                # Add this match to the residuals list
                # [ 'type', 'scan_number', 'parent_mz', 'parent_intensity', 'child_mz', 'child_intensity', 'charge', 'delta', 'residual', 'residual_ppm', 'tag' ]
                row = [ low_mass_ion_name, spectrum['scan_number'], low_mass_ion_mz, spectrum['intensity_array'][match['peak_number']], match['mz'],
                    spectrum['intensity_array'][match['peak_number']], 1, low_mass_ion_mz - match['mz'], -1, match['delta_ppm'], low_mass_ion_name[:2] ]
                spectrum['residuals_list'].append(row)


    ####################################################################################################
    def identify_reporter_ions(self, spectrum):

        tolerance = 20

        for reporter_ion_name,reporter_ion_attributes in self.reporter_ions.items():

            match = self.find_closest_ion(spectrum,reporter_ion_attributes['mz'],tolerance,spectrum_version='input')
            if match:
                if len(spectrum['peak_interpretation'][match['peak_number']]) > 0:
                    spectrum['peak_interpretation'][match['peak_number']] += ', '
                spectrum['peak_interpretation'][match['peak_number']] += reporter_ion_name + "/" + '{:.1f}'.format(match['delta_ppm'])

                spectrum['attributes']['n_identified_reporter_ions'] += 1

                # Add this match to the residuals list
                # [ 'type', 'scan_number', 'parent_mz', 'parent_intensity', 'child_mz', 'child_intensity', 'charge', 'delta', 'residual', 'residual_ppm', 'tag' ]
                row = [ reporter_ion_name, spectrum['scan_number'], reporter_ion_attributes['mz'], spectrum['intensity_array'][match['peak_number']], match['mz'],
                    spectrum['intensity_array'][match['peak_number']], 1, reporter_ion_attributes['mz'] - match['mz'], -1, match['delta_ppm'], reporter_ion_name[:3] ]
                spectrum['residuals_list'].append(row)


    ####################################################################################################
    def identify_isotopes(self, spectrum):

        n_peaks = len(spectrum['mz_array'])
        max_intensity = spectrum['intensity_profile'][1]
        tolerance = self.stats['mz_calibration_tolerance_ppm']
        average_isotope_delta = 1.00291
        max_charge = 4

        # Loop through and identify isotopes
        i_peak = 0
        while i_peak < n_peaks - 1:
            mz = float(spectrum['mz_array'][i_peak])
            intensity = float( spectrum['intensity_array'][i_peak] )
            #print(f"Analyzing peak {i_peak} at {mz}")
            i_lookahead_peak = i_peak + 1
            i_isotope = 1
            i_charge = 0
            charge = max_charge
            done = False
            pursue_more_isotopes = False
            while not done:
                lookahead_mz = float(spectrum['mz_array'][i_lookahead_peak])
                lookahead_intensity = float(spectrum['intensity_array'][i_lookahead_peak])
                diff = lookahead_mz - mz
                delta = diff * charge - i_isotope * average_isotope_delta
                #delta = diff * charge - i_isotope - ( 1 - average_isotope_delta ) * charge
                delta_ppm = delta / mz * 1e6
                abs_delta_ppm = abs(delta_ppm)
                if abs_delta_ppm < tolerance:
                    #print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} matches charge {charge} with abs_delta_ppm={abs_delta_ppm}")
                    spectrum['charge_array'][i_peak] = charge
                    spectrum['n_isotopes_array'][i_peak] += 1
                    spectrum['charge_array'][i_lookahead_peak] = charge
                    spectrum['isotope_array'][i_lookahead_peak] += 1
                    spectrum['isotope_delta'][i_lookahead_peak] = delta
                    spectrum['parent_peak'][i_lookahead_peak] = i_peak

                    if len(spectrum['peak_interpretation'][i_lookahead_peak]) > 0:
                        spectrum['peak_interpretation'][i_lookahead_peak] += ', '
                    spectrum['peak_interpretation'][i_lookahead_peak] += f"i({i_peak})/" + '{:.1f}'.format(delta_ppm)

                    # Store the isotope measurements in the residuals list
                    # [ 'type', 'scan_number', 'parent_mz', 'parent_intensity', 'child_mz', 'child_intensity', 'charge', 'delta', 'residual', 'residual_ppm', 'tag' ]
                    #row = [ 'i', spectrum['scan_number'], mz, intensity, lookahead_mz, lookahead_intensity, charge, diff, delta, delta_ppm, 'i' ]
                    #spectrum['residuals_list'].append(row)

                    pursue_more_isotopes = True
                    done = True

                elif charge == max_charge and delta < 0:
                    #print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is too close even for charge {charge} with delta_ppm={delta_ppm}. Next lookahead peak")
                    charge = max_charge
                    i_lookahead_peak += 1
                    if i_lookahead_peak >= n_peaks:
                        done = True
                elif charge == 1 and delta > 0:
                    #print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is too far even for charge {charge} with delta_ppm={delta_ppm}. Move on to the next peak")
                    done = True
                elif charge == 1 and delta < 0:
                    if debug: print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is closer than the next charge {charge} isotope with delta_ppm={delta_ppm}. Next lookahead peak")
                    charge = max_charge
                    i_lookahead_peak += 1
                    if i_lookahead_peak >= n_peaks:
                        done = True
                else:
                    #print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is not charge {charge} isotope with delta_ppm={delta_ppm}")
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
                    lookahead_mz = float(spectrum['mz_array'][i_lookahead_peak])
                    lookahead_intensity = float(spectrum['intensity_array'][i_lookahead_peak])
                    diff = lookahead_mz - mz
                    delta = diff * charge - i_isotope * average_isotope_delta
                    delta_ppm = delta / mz * 1e6
                    abs_delta_ppm = abs(delta_ppm)
                    if abs_delta_ppm < tolerance:
                        #print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} matches charge {charge} isotope {i_isotope} with abs_delta_ppm={abs_delta_ppm}")
                        spectrum['charge_array'][i_lookahead_peak] = charge
                        spectrum['n_isotopes_array'][i_peak] += 1
                        spectrum['isotope_array'][i_lookahead_peak] += 1
                        spectrum['isotope_delta'][i_lookahead_peak] = delta
                        spectrum['parent_peak'][i_lookahead_peak] = i_peak

                        if len(spectrum['peak_interpretation'][i_lookahead_peak]) > 0:
                            spectrum['peak_interpretation'][i_lookahead_peak] += ', '
                        spectrum['peak_interpretation'][i_lookahead_peak] += f"i({i_peak}{i_isotope})/" + '{:.1f}'.format(delta_ppm)

                        # Store the isotope measurements in the residuals list
                        #row = [ f'i{i_isotope}', spectrum['scan_number'], mz, intensity, lookahead_mz, lookahead_intensity, charge, diff, delta, delta_ppm, 'i' ]
                        #spectrum['residuals_list'].append(row)
                        i_isotope += 1
                    elif delta < 0:
                        #print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is too close for isotope {i_isotope} with delta_ppm={delta_ppm}. Next lookahead peak")
                        pass
                    elif delta > 0:
                        #print(f"  Lookahead peak {i_lookahead_peak} at {lookahead_mz} is too far for isotope {i_isotope} with delta_ppm={delta_ppm}. Done looking.")
                        done = True
                    else:
                        print(f"  Cannot get here!!!!!!!!!!!!!")

                    i_lookahead_peak += 1
                    if i_lookahead_peak >= n_peaks:
                        done = True

            i_peak += 1


    ####################################################################################################
    def clean(self, spectrum, deisotope=False, threshold=0, charge_reduce=False, merge_losses=True):

        # FIXME: redo to this to annotate the origina spectrum with keeps and commands and origins and stuff and then go through and purge the keep=0s
        # Maybe keep a separate dict of flags like "probably_b_ion" and such with a pointer from the original peak number or something. Maybe better as flags in an array?

        n_peaks = spectrum['n_peaks']

        new_peaks_mz = []
        new_peaks_intensity = []
        new_peaks_keep = []
        new_peaks_peak_interpretation = []
        new_peaks_peak_comment = []

        for i_peak in range(0,n_peaks):
            mz = float(spectrum['mz_array'][i_peak])
            charge = float( spectrum['charge_array'][i_peak] )

            if deisotope:
                if spectrum['isotope_array'][i_peak]:
                    spectrum['keep'][i_peak] = 0

            if threshold:
                if spectrum['intensity_array'][i_peak] < threshold:
                    spectrum['keep'][i_peak] = 0

            if merge_losses:
                if spectrum['peak_attributes'][i_peak]['is_neutral_loss']:
                    parent_peak = spectrum['parent_peak'][i_peak]
                    spectrum['intensity_array'][parent_peak] += spectrum['intensity_array'][i_peak]
                    spectrum['peak_comment'][parent_peak] += f"merged "  + spectrum['peak_interpretation'][i_peak]
                    spectrum['keep'][i_peak] = 0

            # If we want to charge reduce and the charge is greater than 1, then reduce it
            elif charge_reduce and spectrum['keep'][i_peak]:
                if charge > 1:
                    mz = mz * charge - (charge - 1) * masses['proton']
                    spectrum['peak_comment'][i_peak] += f"charge reduced {charge}->1"
                    spectrum['keep'][i_peak] = 0
                    new_peaks_mz.append(mz)
                    new_peaks_intensity.append(spectrum['intensity_array'][i_peak])
                    new_peaks_keep.append(1)
                    new_peaks_peak_interpretation.append('')
                    new_peaks_peak_comment.append(f"charge{charge}->1")

        combined_peaks_mz = list(spectrum['mz_array'])
        combined_peaks_intensity = list(spectrum['intensity_array'])
        combined_peaks_keep = list(spectrum['keep'])
        combined_peaks_peak_interpretation = list(spectrum['peak_interpretation'])
        combined_peaks_peak_comment = list(spectrum['peak_comment'])

        combined_peaks_mz.extend(new_peaks_mz)
        combined_peaks_intensity.extend(new_peaks_intensity)
        combined_peaks_keep.extend(new_peaks_keep)
        combined_peaks_peak_interpretation.extend(new_peaks_peak_interpretation)
        combined_peaks_peak_comment.extend(new_peaks_peak_comment)

        column_name_list = [ 'mz', 'intensity', 'keep', 'interpretation', 'comment' ]
        zipped_list =  list(zip(combined_peaks_mz, combined_peaks_intensity, combined_peaks_keep, combined_peaks_peak_interpretation, combined_peaks_peak_comment))
        table = pd.DataFrame(zipped_list, columns = column_name_list)
        table = table[table['keep']==1]
        table.sort_values(by=['mz'], ascending=True, inplace=True, ignore_index=True )

        spectrum['processed'] = {
            'mz_array': table['mz'],
            'intensity_array': table['intensity'],
            'interpretation_array': table['interpretation'],
            'comment_array': table['comment'],
        }



    ####################################################################################################
    def show_spectrum(self, spectrum):

        print(f"Dump of spectrum id {spectrum['id']}")
        print(f"scan={spectrum['scan_number']}")
        print(f"n_peaks={spectrum['n_peaks']}")
        print(f"precursor_mz={spectrum['precursor_mz']}")
        print(f"precursor_charge={spectrum['precursor_charge']}")
        print(f"intensity_profile=",spectrum['intensity_profile'])
        if 'attributes' in spectrum:
            print("Other attributes:")
            for attribute,attribute_value in spectrum['attributes'].items():
                print(f"  {attribute}={attribute_value}")

        print( "i_peak  mz              inten   charge  n_isos  is_iso  delta   n_loss  is_loss interpretation")
        print( "------- --------------- ------- ------- ------- ------- ------- ------- ------- -----------------------")
        #print("4       129.0659        24.9    0       0       0       0.0000  0       0       b2(GA)/0.6")

        rescale = 2
        i = 0
        for i in range(len(spectrum['mz_array'])):
            mz = float(spectrum['mz_array'][i])
            if rescale == 1:
                intensity = int(spectrum['intensity_array'][i] / 300 )
                row_data = [ str(i), '{:.4f}'.format(mz), '{:d}'.format(intensity) ]
            elif rescale == 2:
                intensity = float(spectrum['intensity_array'][i] / spectrum['intensity_profile'][1] * 100)
                row_data = [ str(i), '{:.4f}'.format(mz), '{:.1f}'.format(intensity) ]
            else:
                intensity = float(spectrum['intensity_array'][i])
                row_data = [ str(i), '{:.4f}'.format(mz), '{:.1f}'.format(intensity) ]
            if 'charge_array' in spectrum:
                row_data.extend( [ str(spectrum['charge_array'][i]), str(spectrum['n_isotopes_array'][i]), str(spectrum['isotope_array'][i]), 
                    '{:.4f}'.format(spectrum['isotope_delta'][i]) ] )
            if 'peak_attributes' in spectrum:
                row_data.extend( [ str(spectrum['peak_attributes'][i]['n_neutral_losses']), str(spectrum['peak_attributes'][i]['is_neutral_loss']) ] )
            if 'peak_interpretation' in spectrum:
                row_data.extend( [ spectrum['peak_interpretation'][i] + '  --  ' + spectrum['peak_comment'][i] ] )
            print("\t".join(row_data) )

        i = 0
        if 'processed' in spectrum:
            print("============================================\n")
            for i in range(len(spectrum['processed']['mz_array'])):
                mz = float(spectrum['processed']['mz_array'][i])
                if rescale == 1:
                    intensity = int(spectrum['processed']['intensity_array'][i] / 300 )
                    row_data = [ str(i), '{:.4f}'.format(mz), '{:d}'.format(intensity) ]
                elif rescale == 2:
                    intensity = float(spectrum['processed']['intensity_array'][i] / spectrum['intensity_profile'][1] * 100)
                    row_data = [ str(i), '{:.4f}'.format(mz), '{:.1f}'.format(intensity) ]
                else:
                    intensity = float(spectrum['processed']['intensity_array'][i])
                    row_data = [ str(i), '{:.4f}'.format(mz), '{:.1f}'.format(intensity) ]
                row_data.extend( [ spectrum['processed']['interpretation_array'][i] + '  --  ' + spectrum['processed']['comment_array'][i] ] )
                print("\t".join(row_data) )


    ####################################################################################################
    def compute_intensity_profile(self, intensities):

        intensity_array = numpy.asarray(intensities)
        sorted_array = numpy.sort(intensity_array)
        if len(intensities) >= 10:
            noise = sorted_array[3]
            signal = sorted_array[-3]
        else:
            index = int(len(intensities)/5)
            noise = sorted_array[index]
            signal = sorted_array[-1*index]

        # Prevent division by 0
        if noise == 0:
            noise = 0.01
        signal_to_noise = signal / noise * 3.0
        return [sorted_array[0], sorted_array[-1], noise, signal, signal_to_noise]


    ####################################################################################################
    def interpret_de_novo(self, spectrum):

        n_peaks = len(spectrum['processed']['mz_array'])
        spectrum['interpretations'] = {}

        self.read_de_novo(spectrum, start='y')


    ####################################################################################################
    def read_de_novo(self, spectrum, start=None):

        n_peaks = len(spectrum['processed']['intensity_array'])
        spectrum['proc_interpretations'] = {}
        precursor_mz = spectrum['precursor_mz']
        precursor_charge = spectrum['precursor_charge']
        precursor_mw = ( precursor_mz - masses['proton'] ) * precursor_charge
        tolerance = 5
        verbose = 0

        #### Set up a structure to hold intermediate sequencing results for later stitching together
        partial_sequencing_index = {}

        self.index_peaks(spectrum)

        fragmentation_type = 'HCD'

        # Create a list for all possible sequences to emerge from here
        sequences_list = []
        finished_sequences_list = []

        ion_series_list = [ 'b', 'y']
        terminal_modifications = {
            'nterm': { 'n-': 0.0, '[Acetyl]-': 42.010565 },
            #'nterm': { 'n-': 0.0, '[Acetyl]-': 42.010565, '[Carbamidomethyl]-': 57.021464 },   # [Carbamidomethyl]- is the same as having a G!
            #'nterm': { 'n-': 0.0, '[TMT6]-': 229.162932, '[iTRAQ4]-': 144.102063 },
            'cterm': { '-c': 0.0 },
            #'cterm': { '-c': 0.0, '[AlphaAmidation]':-58.005479, '[Methylthio]-': 45.987721 },
        }
        terminal_reversal = {
            'nterm': 'cterm',
            'cterm': 'nterm',
        }

        # Set up the list of sequences to start with
        for ion_series in ion_series_list:
            terminus_type = ion_series_attr[ion_series]['terminus_type']
            for terminal_modification in terminal_modifications[terminus_type]:
                sequence = {
                    'terminus_type': terminus_type,
                    'working_position': 0,
                    'working_mass': ion_series_attr[ion_series]['mass'] + terminal_modifications[terminus_type][terminal_modification] + masses['proton'],
                    'complement_mass': precursor_mw - ion_series_attr[ion_series]['mass'] - ion_series_attr[ion_series_attr[ion_series]['complement']]['mass'],
                    'residues': [ terminal_modification ],
                    'preliminary_score': 0.0,
                    'status': 'GO',
                }
                sequences_list.append(sequence)

        # Main loop to try to find the next possible residues
        done = False
        while not done:
            new_sequences_list = []
            for sequence in sequences_list:

                if verbose: print(f"Looking to extend " + ''.join(sequence['residues']))
                found_match = False

                # If we're starting from the cterm and we might be nearing the nterm, look for completion with nterm residues
                aa_masses = self.aa_masses
                if sequence['terminus_type'] == 'cterm' and sequence['working_position'] != 0:
                     aa_masses = self.nterm_aa_masses

                #### For each of the potential terminal modifications on the opposite end that we're starting from
                for terminal_modification_label,terminal_modification_mass in terminal_modifications[terminal_reversal[sequence['terminus_type']]].items():

                    #### And for each potential residue that might be the last one, see if we can complete the sequence
                    for residue,residue_mass in aa_masses.items():
                        # If this residue is the missing piece, declare victory
                        unexplained_mass_ppm = (sequence['complement_mass'] - residue_mass - terminal_modification_mass) / precursor_charge / precursor_mz * 1e6
                        if verbose: print(f"For residue {residue}, we still have {unexplained_mass_ppm} ppm to go")
                        if abs(unexplained_mass_ppm) < tolerance:
                            if verbose: print("Complete!")
                            new_sequence = copy.deepcopy(sequence)
                            new_sequence['working_position'] += 1
                            new_sequence['working_mass'] += residue_mass + terminal_modification_mass
                            new_sequence['complement_mass'] -= residue_mass + terminal_modification_mass
                            new_sequence['preliminary_score'] += 1.0
                            new_sequence['residues'].append(terminal_modification_label + residue)
                            new_sequence['status'] = 'COMPLETE'
                            finished_sequences_list.append(new_sequence)
                            found_match = True
                            continue

                #### If we didn't manange to complete all the way to the other terminus, see which residues might take us another step
                if not found_match:

                    #### Set the usual residues, except if we're just starting out from the nterm, then use possible nterm residues
                    aa_masses = self.aa_masses
                    if sequence['terminus_type'] == 'nterm' and sequence['working_position'] == 0:
                        aa_masses = self.nterm_aa_masses

                    #### For each potential residue that might be the next one, see if we can take one more step
                    for residue,residue_mass in aa_masses.items():
                        lookahead_mz = sequence['working_mass'] + residue_mass
                        if verbose: print(f"  Look for {residue} at {lookahead_mz}...")
                        match = self.find_closest_ion(spectrum,lookahead_mz,tolerance)
                        if match:
                            if verbose: print(f"  Found match for {residue} at {match['mz']} with delta_ppm {match['delta_ppm']} ppm")
                            new_sequence = copy.deepcopy(sequence)
                            new_sequence['working_position'] += 1
                            new_sequence['working_mass'] = match['mz']
                            new_sequence['complement_mass'] -= residue_mass
                            new_sequence['preliminary_score'] += match['score']
                            new_sequence['residues'].append(residue)
                            new_sequences_list.append(new_sequence)
                            found_match = True

                            #### Store this partial sequencing result in the partial_sequencing_index for possible later stitching
                            integer_index_mass = int(new_sequence['working_mass'])
                            if integer_index_mass not in partial_sequencing_index:
                                partial_sequencing_index[integer_index_mass] = []
                            partial_sequencing_index[integer_index_mass].append(new_sequence)

                #### If we didn't find a next step to take, or always at the nterm, try looking ahead two steps
                #### There's a danger here that if a single step leads us astray, then we don't have the opportunity to recover the correct chain with a double hope if we never try it
                if not found_match or ( sequence['terminus_type'] == 'nterm' and sequence['working_position'] == 0 ):

                    # If we're starting from the cterm and we might be nearing the nterm, look for completion with nterm residues
                    double_aa_masses = self.double_aa_masses
                    if sequence['terminus_type'] == 'cterm' and sequence['working_position'] != 0:
                        double_aa_masses = self.nterm_double_aa_masses

                    #### For each of the potential terminal modifications on the opposite end that we're starting from
                    for terminal_modification_label,terminal_modification_mass in terminal_modifications[terminal_reversal[sequence['terminus_type']]].items():

                        #### And for each potential residue pair that might be the last one, see if we can complete the sequence
                        #### FIXME. Maybe this should be done for a triplet eventually as well
                        for residue,residue_mass in double_aa_masses.items():
                            # If this residue is the missing piece, declare victory
                            unexplained_mass_ppm = (sequence['complement_mass'] - residue_mass - terminal_modification_mass) / precursor_charge / precursor_mz * 1e6
                            if verbose: print(f"For residue {terminal_modification_label}({residue}), we still have {unexplained_mass_ppm} ppm to go")
                            if abs(unexplained_mass_ppm) < tolerance:
                                if verbose: print("Complete!")
                                new_sequence = copy.deepcopy(sequence)
                                new_sequence['working_position'] += 2
                                new_sequence['working_mass'] += residue_mass + terminal_modification_mass
                                new_sequence['complement_mass'] -= residue_mass + terminal_modification_mass
                                new_sequence['preliminary_score'] += 1.0
                                new_sequence['residues'].append(f"{terminal_modification_label}({residue})")
                                new_sequence['status'] = 'COMPLETE'
                                finished_sequences_list.append(new_sequence)
                                found_match = True

                    #### If we didn't manange to complete all the way to the other terminus, see which residue pair might take us another step
                    if not found_match:

                        #### Set the usual residue pairs, except if we're just starting out from the nterm, then use possible nterm residues
                        double_aa_masses = self.double_aa_masses
                        if sequence['terminus_type'] == 'nterm' and sequence['working_position'] == 0:
                            double_aa_masses = self.nterm_double_aa_masses

                        #### For each potential residue that might be the next one, see if we can take one more step
                        for residue,residue_mass in double_aa_masses.items():
                            lookahead_mz = sequence['working_mass'] + residue_mass
                            if verbose: print(f"  Look for {residue} at {lookahead_mz}...")
                            match = self.find_closest_ion(spectrum,lookahead_mz,tolerance)
                            if match:
                                if verbose: print(f"  Found match for {residue} with delta_ppm {match['delta_ppm']} ppm")
                                new_sequence = copy.deepcopy(sequence)
                                new_sequence['working_position'] += 2
                                new_sequence['working_mass'] = match['mz']
                                new_sequence['complement_mass'] -= residue_mass
                                new_sequence['preliminary_score'] += match['score']
                                new_sequence['residues'].append(f"({residue})")
                                new_sequences_list.append(new_sequence)
                                found_match = True

                #### If we haven't found anywhere to continue, then we're just as a dead end here with nowhere farther to go
                if not found_match:
                    sequence['status'] = 'STOPPED'
                    finished_sequences_list.append(sequence)

            #### If there are sequences that need more work, copy them back to the primary sequences list
            if len(new_sequences_list) > 0:
                if verbose: print("Copying new_sequences_list to sequences_list")
                sequences_list = new_sequences_list
            #### Or else if there are no new sequences then we're done
            else:
                done = True

            # Relief valve for now
            if len(sequences_list) > 1000:
                print("WARNING: Exceeding 1000 test sequences. Giving up for now..")
                return


        #return

        # Print out the results
        for terminus_type in [ 'cterm', 'nterm' ]:
            print(f"{terminus_type} sequencing results:")
            sequences_list = []
            preliminary_scores_list = []
            for sequence in finished_sequences_list:
                if sequence['terminus_type'] == terminus_type:
                    sequence['preliminary_score'] /= ( precursor_mw / 1000 )
                    if terminus_type == 'cterm':
                        residues = sequence['residues']
                        residues.reverse()
                        complement_mass = '{:.4f}'.format(sequence['complement_mass'])
                        if sequence['status'] == 'COMPLETE':
                            sequences_list.append(''.join(residues))
                            preliminary_scores_list.append(sequence['preliminary_score'])
                        else:
                            sequences_list.append(f"n-(+{complement_mass})" + ''.join(residues))
                            preliminary_scores_list.append(sequence['preliminary_score'])
                            #### See if there's something with the complement mass
                            integer_index_mass = int(sequence['complement_mass'])
                            if integer_index_mass in partial_sequencing_index:
                                for potential_sequence_match in partial_sequencing_index[integer_index_mass]:
                                    complement_mass = '{:.4f}'.format(potential_sequence_match['complement_mass'])
                                    potential_sequence_residues = sequence['residues']
                                    potential_sequence_residues.reverse()
                                    print("    " + f"n-(+{complement_mass})" + ''.join(potential_sequence_residues))
                    else:
                        residues = sequence['residues']
                        complement_mass = '{:.4f}'.format(sequence['complement_mass'])
                        if sequence['status'] == 'COMPLETE':
                            sequences_list.append(''.join(residues) + "-c")
                            preliminary_scores_list.append(sequence['preliminary_score'])
                        else:
                            sequences_list.append(''.join(residues) + f"(+{complement_mass})-c")
                            preliminary_scores_list.append(sequence['preliminary_score'])

            #### Sort the final list by reverse score
            column_name_list = [ 'score', 'sequence' ]
            zipped_list =  list(zip(preliminary_scores_list, sequences_list))
            table = pd.DataFrame(zipped_list, columns = column_name_list)
            table.sort_values(by=['score'], ascending=False, inplace=True, ignore_index=True )

            scores_list = table['score']
            sequences_list = table['sequence']
            i_sequence = 0
            for score in scores_list:
                print("\t".join( [ '{:.3f}'.format(score),sequences_list[i_sequence] ] ))
                i_sequence += 1
                if i_sequence > 100:
                    break


    ####################################################################################################
    # Find the closest ion in the spectrum
    def find_closest_ion(self,spectrum,lookahead_mz,tolerance,spectrum_version='processed'):

        # Determine which spectrum to look in, either 'input' or 'processed'
        lookup_dict = 'proc_mz_dict'
        if spectrum_version == 'input':
            lookup_dict = 'mz_dict'

        # Override the input tolerance with the reference tolerance
        tolerance = self.stats['mz_calibration_tolerance_ppm']

        # Get the integer mass as a dict key
        int_lookahead_mz = int(lookahead_mz)
        if int_lookahead_mz not in spectrum[lookup_dict]:
            return

        # Loop over all the peaks in this bin
        # FIXME: This just takes the first one within the tolerance, not necessarily the best one if there's more than one within tolerance!
        for peak_data in spectrum[lookup_dict][int_lookahead_mz]:
            peak_mz = peak_data['mz']
            peak_intensity = peak_data['intensity']
            delta = peak_mz - lookahead_mz
            delta_ppm = delta / lookahead_mz * 1e6
            abs_delta_ppm = abs(delta_ppm)
            if abs_delta_ppm < tolerance:

                # Hack a delta score. This is basically flat at 1.0 between 0 and sigma and then linear down to 0.0 at tolerance
                sigma = self.stats['mz_calibration_sigma_ppm']
                delta_score = -1.0/(tolerance-sigma)*abs_delta_ppm + 1 + sigma/(tolerance-sigma)
                #delta_score = ( tolerance - abs_delta_ppm ) / ( tolerance / 2.0 )
                if delta_score > 1.0:
                    delta_score = 1.0

                # Hack an intensity score
                minimum_intensity = 1500
                maximum_intensity = 15000
                intensity_score = ( peak_intensity - minimum_intensity ) / ( maximum_intensity - minimum_intensity )
                if intensity_score > 1.0:
                    intensity_score = 1.0
                if intensity_score < 0.1:
                    intensity_score = 0.1
                match = { 'mz': peak_mz, 'delta_ppm': delta_ppm, 'score': delta_score * intensity_score, 'peak_number': peak_data['peak_number'] }
                return match


    ####################################################################################################
    # Put all the peaks into a dict by integer mass to make lookups faster
    def index_peaks(self,spectrum):

        # If there isn't yet an index dict for the input spectrum, then create the index dict
        if 'mz_dict' not in spectrum:
            spectrum['mz_dict'] = {}
            i_peak = 0
            for peak_mz in spectrum['mz_array']:
                int_peak_mz = int(peak_mz)
                if int_peak_mz not in spectrum['mz_dict']:
                    spectrum['mz_dict'][int_peak_mz] = []
                peak_info = { 'mz': peak_mz, 'intensity': spectrum['intensity_array'][i_peak], 'peak_number': i_peak }
                spectrum['mz_dict'][int_peak_mz].append(peak_info)
                i_peak += 1

        # If there is a processed spectrum available, then create or rebuild the index dict
        if 'processed' in spectrum:
            spectrum['proc_mz_dict'] = {}
            i_peak = 0
            for peak_mz in spectrum['processed']['mz_array']:
                int_peak_mz = int(peak_mz)
                if int_peak_mz not in spectrum['proc_mz_dict']:
                    spectrum['proc_mz_dict'][int_peak_mz] = []
                peak_info = { 'mz': peak_mz, 'intensity': spectrum['processed']['intensity_array'][i_peak], 'peak_number': i_peak }
                spectrum['proc_mz_dict'][int_peak_mz].append(peak_info)
                i_peak += 1


    ####################################################################################################
    #### Parse the filter string
    def parse_filter_string(self,filter_string,stats):

        fragmentation_tag = '??'

        mass_accuracy = '??'
        match = re.search(r'^(\S+)',filter_string)
        if match.group(0) == 'FTMS':
            mass_accuracy = 'HR'
        elif match.group(0) == 'ITMS':
            mass_accuracy = 'LR'
        else:
            self.log_event('ERROR','UnrecognizedDetector',f"Unrecognized detector '{match.group(0)}' in filter string {filter_string}")
            return

        ms_level = 0
        match = re.search(r' ms ',filter_string)
        if match:
            ms_level = 1
        else:
            match = re.search(r' ms(\d) ',filter_string)
            if match:
                ms_level = int(match.group(1))
            else:
                self.log_event('ERROR','UnrecognizedMSLevel',f"Unrecognized MS level in filter string {filter_string}")
                return

        # See if there's fragmentation of a known type
        have_cid = 0
        have_etd = 0
        have_hcd = 0
        have_fragmentation = 0
        match = re.search(r'@hcd',filter_string)
        if match: have_hcd = 1
        match = re.search(r'@cid',filter_string)
        if match: have_cid = 1
        match = re.search(r'@etd',filter_string)
        if match: have_etd = 1
        match = re.search(r'@',filter_string)
        if match: have_fragmentation = 1
        dissociation_sum = have_cid + have_etd + have_hcd
        if have_fragmentation > dissociation_sum:
            self.log_event('ERROR','UnrecognizedFragmentation',f"Unrecognized string after @ in filter string '{filter_string}'")
            return

        #### If this is an MS1 scan, learn what we can from it
        if ms_level == 0:
            if 'n_ms0_spectra' not in stats:
                stats['n_ms0_spectra'] = 0
            stats['n_ms0_spectra'] += 1
            self.log_event('ERROR','UnknownMSLevel',f"Unable to determine MS level in filter string '{filter_string}'")
            return

        #### If this is an MS1 scan, learn what we can from it
        if ms_level == 1:
            stats['n_ms1_spectra'] += 1
            spectrum_type = 'MS1'
            if mass_accuracy == 'HR':
                if stats['high_accuracy_precursors'] == 'unknown':
                    stats['high_accuracy_precursors'] = 'true'
                elif stats['high_accuracy_precursors'] == 'true' or stats['high_accuracy_precursors'] == 'multiple':
                    pass
                elif stats['high_accuracy_precursors'] == 'false':
                    stats['high_accuracy_precursors'] = 'multiple'
                else:
                    self.log_event('ERROR','InternalStateError',f"high_accuracy_precursors has a strange state '{stats['high_accuracy_precursors']}'")
                    return
            elif mass_accuracy == 'LR':
                if stats['high_accuracy_precursors'] == 'unknown':
                    stats['high_accuracy_precursors'] = 'false'
                elif stats['high_accuracy_precursors'] == 'false' or stats['high_accuracy_precursors'] == 'multiple':
                    pass
                elif stats['high_accuracy_precursors'] == 'true':
                    stats['high_accuracy_precursors'] = 'multiple'
                else:
                    self.log_event('ERROR','InternalStateError',f"high_accuracy_precursors has a strange state '{stats['high_accuracy_precursors']}'")
                    return
            else:
                self.log_event('ERROR','UnknownPrecursorScanType',f"Unknown precursor scan type '{match.group(0)}'")
                return

        #### Else it's an MSn (n>1) scan so learn what we can from that
        elif ms_level == 2:
            dissociation_sum = have_cid + have_etd + have_hcd
            stats['n_ms2_spectra'] += 1
            if dissociation_sum == 1:
                if have_hcd:
                    stats[f"n_{mass_accuracy}_HCD_spectra"] += 1
                    spectrum_type = f"{mass_accuracy}_HCD"
                    fragmentation_tag = f"{mass_accuracy} HCD"
                elif have_cid:
                    stats[f"n_{mass_accuracy}_IT_CID_spectra"] += 1
                    spectrum_type = f"{mass_accuracy}_IT_CID"
                    fragmentation_tag = f"{mass_accuracy} IT CID"
                elif have_etd:
                    stats[f"n_{mass_accuracy}_IT_ETD_spectra"] += 1
                    spectrum_type = f"{mass_accuracy}_IT_ETD"
                    fragmentation_tag = f"{mass_accuracy} IT ETD"
                else:
                    self.log_event('ERROR','UnrecognizedFragmentation',f"Did not find @hcd, @cid, @etd in filter string '{filter_string}'")
                    return

            elif dissociation_sum == 2:
                if have_etd and have_hcd:
                    stats[f"n_{mass_accuracy}_EThcD_spectra"] += 1
                    spectrum_type = f"{mass_accuracy}_EThcD"
                    fragmentation_tag = f"{mass_accuracy} EThcD"
                elif have_etd and have_cid:
                    stats[f"n_{mass_accuracy}_ETciD_spectra"] += 1
                    spectrum_type = f"{mass_accuracy}_ETciD"
                    fragmentation_tag = f"{mass_accuracy} ETciD"
                else:
                    self.log_event('ERROR','UnrecognizedFragmentation',f"Did not find known fragmentation combination in filter string '{filter_string}'")
                    return

            else:
                self.log_event('ERROR','CantParseFilterStr',f"Unable to determine dissociation type in filter string '{filter_string}'")
                return

            #### Record the fragmentation type for this run
            if spectrum_type != stats['fragmentation_type']:
                if stats['fragmentation_type'] == 'unknown':
                    stats['fragmentation_type'] = spectrum_type
                elif stats['fragmentation_type'] == 'multiple':
                    pass
                else:
                    stats['fragmentation_type'] = 'multiple'
                    self.log_event('ERROR','MultipleFragTypes',f"There are multiple fragmentation types in this MS run. Split them.")


        # If we have ms_level > 2, then make a note, but this is uncharted waters
        else:
            if 'n_ms3+_spectra' not in stats:
                stats['n_ms3+_spectra'] = 0
            stats['n_ms3+_spectra'] += 1

        # Store the fragmentation_tag for use by the caller
        stats['fragmentation_tag'] = fragmentation_tag


    ####################################################################################################
    #### Add spectrum
    #def add_spectrum(self,spectrum,spectrum_type,precursor_mz):
    def add_spectrum(self,params):

        stats, precursor_mz, peaklist = params
        spectrum_type = stats['fragmentation_type']

        if 'n_spectra_examined' not in self.metadata['files'][self.mzml_file]['spectra_stats']:
            self.metadata['files'][self.mzml_file]['spectra_stats']['n_spectra_examined'] = 0
        self.metadata['files'][self.mzml_file]['spectra_stats']['n_spectra_examined'] += 1

        destination = f"lowend_{spectrum_type}"
        #destination2 = f"neutral_loss_{spectrum_type}"

        if destination not in self.composite:
            if spectrum_type == 'HR_HCD':
                #print(f"INFO: Creating a composite spectrum {destination}")
                minimum = 100
                maximum = 400
                binsize = 0.0002
                array_size = int( (maximum - minimum ) / binsize ) + 1
                self.composite[destination] = { 'minimum': minimum, 'maximum': maximum, 'binsize': binsize }
                self.composite[destination]['intensities'] = numpy.zeros(array_size,dtype=numpy.float32)
                self.composite[destination]['n_peaks'] = numpy.zeros(array_size,dtype=numpy.int32)

                #print(f"INFO: Creating a composite spectrum {destination2}")
                # minimum = 0
                # maximum = 165
                # binsize = 0.001
                # array_size = int( (maximum - minimum ) / binsize ) + 1
                # self.composite[destination2] = { 'minimum': minimum, 'maximum': maximum, 'binsize': binsize }
                # self.composite[destination2]['intensities'] = numpy.zeros(array_size,dtype=numpy.float32)
                # self.composite[destination2]['n_peaks'] = numpy.zeros(array_size,dtype=numpy.int32)

            elif spectrum_type == 'LR_IT_CID':
                #print(f"INFO: Creating a composite spectrum {destination}")
                minimum = 100
                maximum = 460
                binsize = 0.05
                array_size = int( (maximum - minimum ) / binsize ) + 1
                self.composite[destination] = { 'minimum': minimum, 'maximum': maximum, 'binsize': binsize }
                self.composite[destination]['intensities'] = numpy.zeros(array_size,dtype=numpy.float32)
                self.composite[destination]['n_peaks'] = numpy.zeros(array_size,dtype=numpy.int32)

            elif spectrum_type == 'LR_IT_ETD':
                #print(f"INFO: Creating a composite spectrum {destination}")
                minimum = 100
                maximum = 400
                binsize = 0.05
                array_size = int( (maximum - minimum ) / binsize ) + 1
                self.composite[destination] = { 'minimum': minimum, 'maximum': maximum, 'binsize': binsize }
                self.composite[destination]['intensities'] = numpy.zeros(array_size,dtype=numpy.float32)
                self.composite[destination]['n_peaks'] = numpy.zeros(array_size,dtype=numpy.int32)

            else:
                #print(f"ERROR: Unrecognized spectrum type {spectrum_type}")
                return

        #### Convert the peak lists to numpy arrays
        intensity_array = numpy.asarray(peaklist['intensity_array'])
        mz_array = numpy.asarray(peaklist['mz_array'])

        #### Limit the numpy arrays to the region we're interested in
        intensity_array = intensity_array[ ( mz_array > self.composite[destination]['minimum'] ) & ( mz_array < self.composite[destination]['maximum'] ) ]
        mz_array = mz_array[ ( mz_array > self.composite[destination]['minimum'] ) & ( mz_array < self.composite[destination]['maximum'] ) ]

        #### Compute their bin locations and store n_peaks and intensities
        bin_array = (( mz_array - self.composite[destination]['minimum'] ) / self.composite[destination]['binsize']).astype(int)
        self.composite[destination]['n_peaks'][bin_array] += 1
        self.composite[destination]['intensities'][bin_array] += intensity_array

        #### The stuff below was trying to compute neutral losses. It works, but it's too slow. Need a Cython implementation I guess

        # sorted_mz_intensity = sorted(mz_intensity, key=lambda pair: pair[2], reverse=True)
        # length = 30
        # if len(mz_intensity) < 30:
            # length = len(mz_intensity)
        # clipped_mz_intensity = sorted_mz_intensity[0:(length-1)]
        # mz_intensity = sorted(clipped_mz_intensity, key=lambda pair: pair[1], reverse=False)
        # i = 0
        # while i < length-1:
            # mz_intensity[i][0] = i
            # i += 1
        # sorted_mz_intensity = sorted(mz_intensity, key=lambda pair: pair[2], reverse=True)


        # imajor_peak = 0
        # while imajor_peak < 20 and imajor_peak < len(mz_intensity):
            # ipeak,mz,major_intensity = sorted_mz_intensity[imajor_peak]
            # #print(ipeak,mz,major_intensity)
            # if major_intensity == 0.0:
                # major_intensity = 1.0
            # imz = mz
            # while ipeak > 0 and mz - imz < self.composite[destination2]['maximum'] and imz > 132 and abs(imz - precursor_mz) > 5:
                # bin = int( ( mz - imz ) / self.composite[destination2]['binsize'])
                # #print("  - ",ipeak,imz,bin,intensity)
                # self.composite[destination2]['n_peaks'][bin] += 1
                # self.composite[destination2]['intensities'][bin] += intensity
                # ipeak -= 1
                # imz = mz_intensity[ipeak][1]
                # intensity = mz_intensity[ipeak][2]
            # imajor_peak += 1


    ####################################################################################################
    #### Read header
    #### Open the mzML file and read line by line into a text buffer that we will XML parse
    def read_header(self):

        #### If the mzML is gzipped, then open with zlib, else a plain open
        match = re.search(r'\.gz$',self.mzml_file)
        if match:
            infile = gzip.open(self.mzml_file)
        else:
            infile = open(self.mzml_file)

        #### set up a text buffer to hold the mzML header
        buffer = ''
        counter = 0

        #### Read line by line
        for line in infile:

            if not isinstance(line, str):
                line = str(line, 'utf-8', 'ignore')

            #### Completely skip the <indexedmzML> tag if present
            match = re.search('<indexedmzML ',line)
            if match: continue

            #### Look for first tag after the header and end when found
            match = re.search('<run ',line)
            if match: break
            if counter > 0: buffer += line
            counter += 1

        #### Close file and process the XML in the buffer
        infile.close()

        #### Finish the XML by closing the tags
        buffer += '  </mzML>\n'

        #### Get the root of the XML and the namespace
        xmlroot = etree.fromstring(buffer)
        namespace = xmlroot.nsmap
        #print(namespace)
        if None in namespace:
            namespace = '{'+namespace[None]+'}'
        else:
            namespace = ''

        #### Create a reference of instruments we know
        instrument_by_category = {
            'pureHCD': [
                'MS:1001911|Q Exactive',
                'MS:1002526|Q Exactive Plus',
                'MS:1002523|Q Exactive HF',
                'MS:1002877|Q Exactive HF-X'
             ],
            'ion_trap': [
                'MS:1000447|LTQ',
                'MS:1000638|LTQ XL ETD',
                'MS:1000854|LTQ XL',
                'MS:1000855|LTQ Velos',
                'MS:1000856|LTQ Velos ETD',
                'MS:1000167|LCQ Advantage',
                'MS:1000168|LCQ Classic',
                'MS:1000169|LCQ Deca XP Plus',
                'MS:1000554|LCQ Deca',
                'MS:1000578|LCQ Fleet'
             ],
            'variable': [ 
                'MS:1000448|LTQ FT', 
                'MS:1000557|LTQ FT Ultra', 
                'MS:1000449|LTQ Orbitrap',
                'MS:1000555|LTQ Orbitrap Discovery', 
                'MS:1000556|LTQ Orbitrap XL',
                'MS:1000639|LTQ Orbitrap XL ETD',
                'MS:1001910|LTQ Orbitrap Elite', 
                'MS:1001742|LTQ Orbitrap Velos', 
                'MS:1002835|LTQ Orbitrap Classic', 
                'MS:1002416|Orbitrap Fusion', 
                'MS:1002417|Orbitrap Fusion ETD', 
                'MS:1002732|Orbitrap Fusion Lumos',
                'MS:1003028|Orbitrap Exploris 480',
                'MS:1000483|Thermo Fisher Scientific instrument model'
             ] }

        #### Restructure it into a dict by PSI-MS identifier
        instrument_attributes = {}
        for instrument_category in instrument_by_category:
            for instrument_string in instrument_by_category[instrument_category]:
                accession,name = instrument_string.split('|')
                instrument_attributes[accession] = { 'category': instrument_category, 'accession': accession, 'name': name }

        #### Get all the CV params in the header and look for ones we know about
        cv_params = xmlroot.findall(f'.//{namespace}cvParam')
        found_instrument = 0
        for cv_param in cv_params:
            accession = cv_param.get('accession')
            if accession in instrument_attributes:
                #### Store the attributes about this instrument model
                model_data = { 'accession': accession, 'name': instrument_attributes[accession]['name'], 'category': instrument_attributes[accession]['category'] }
                self.metadata['files'][self.mzml_file]['instrument_model'] = model_data
                found_instrument = 1

        #### If none are an instrument we know about about, ask for help
        if not found_instrument:
            self.log_event('ERROR','UnrecogInstr',f"Did not recognize the instrument. Please teach me about this instrument.")
            return


    ####################################################################################################
    #### Assess all the fiducial peaks
    def assess_fiducial_peaks(self, table):

        tolerance = self.user_parameters['first_pass_tolerance']
        min = -1 * tolerance
        max = tolerance
        binsize = 0.2
        n_bins = int( (max - min) / binsize )

        for ion_name,ion in self.reporter_ions.items():
            #print(f"  Assessing data for peak {reporter_ion_name}")
            self.reference_peak_attrs[ion_name] = {}
            table_selected = table[table['type']==ion_name]

            residual_ppm = table_selected['residual_ppm']
            y,x = numpy.histogram(residual_ppm,range=[min,max],bins=n_bins)
            x=x[0:(len(y))]+binsize/2.0

            #### Look for and fit the largest peak
            peak = self.find_peak(x,y,0.0)
            detection_result = '-'
            if peak['assessment']['is_found']:
                detection_result = 'mz: ' + '{:6.2f}'.format(peak['fit']['mz']) + f"   intensity: {int(peak['fit']['intensity'])}"
            print(f"{ion_name}: {detection_result}")
            self.reference_peak_attrs[ion_name]['peak'] = peak

        for ion_name,ion in self.low_mass_ions.items():
            #print(f"  Assessing data for peak {reporter_ion_name}")
            self.reference_peak_attrs[ion_name] = {}
            table_selected = table[table['type']==ion_name]

            residual_ppm = table_selected['residual_ppm']
            y,x = numpy.histogram(residual_ppm,range=[min,max],bins=n_bins)
            x=x[0:(len(y))]+binsize/2.0

            #### Look for and fit the largest peak
            peak = self.find_peak(x,y,0.0)
            detection_result = '-'
            if peak['assessment']['is_found']:
                detection_result = 'mz: ' + '{:6.2f}'.format(peak['fit']['mz']) + f"   intensity: {int(peak['fit']['intensity'])}"
            if detection_result != '-':
                print(f"{ion_name}: {detection_result}")
            self.reference_peak_attrs[ion_name]['peak'] = peak

        #self.metadata['files'][self.mzml_file]['ROIs'] = ROIs


    ####################################################################################################
    #### Perform a mass calibration based on expected low mass ions
    def perform_mass_calibration(self, table):

        # For now, manually select a subset of peaks to use as the calibrator
        selected_calibration_ions = { 'y1(K)': 1, 'y1(K)-H20': 1, 'y1(R)': 1 }

        # Build a table of peak data to use for calibration, combining the selected ions
        i_loop = 0
        for ion_name,ion in self.low_mass_ions.items():
            if ion_name in selected_calibration_ions:
                if i_loop == 0:
                    table_selected = table[table['type']==ion_name]
                else:
                    table_selected.extend(table[table['type']==ion_name])

        # Set up parameters for the histogram
        tolerance = self.user_parameters['first_pass_tolerance']
        min = -1 * tolerance
        max = +1 * tolerance
        binsize = 0.2
        n_bins = int( (max - min) / binsize )

        # Build the histogram of residuals
        residual_ppm = table_selected['residual_ppm']
        y,x = numpy.histogram(residual_ppm,range=[min,max],bins=n_bins)
        x=x[0:(len(y))]+binsize/2.0

        # Look for and fit the largest peak in the histogram
        peak = self.find_peak(x,y,0.0)

        # If a decent peak is found, update the calibration parameters for subsequence use
        if peak['assessment']['is_found']:
            self.stats['mz_calibration_offset_ppm'] = peak['fit']['mz']
            self.stats['mz_calibration_sigma_ppm'] = peak['fit']['sigma']
            self.stats['mz_calibration_tolerance_ppm'] = peak['fit']['sigma'] * 4


    ####################################################################################################
    #### Peak finding routine
    def find_peak(self,x,y,reference_x):

        show_interactive_plots = 0

        if show_interactive_plots:
            plt.step(x,y)
            plt.title("Input region")
            plt.show()

        peak = { 'mode_bin': { 'ibin': 0, 'mz': 0, 'intensity': 0 } }
        for ibin in range(0,len(x)-1):
            if y[ibin] > peak['mode_bin']['intensity']:
                peak['mode_bin']['ibin'] = ibin
                peak['mode_bin']['mz'] = x[ibin]
                peak['mode_bin']['intensity'] = y[ibin]

        #### If there were 0 peaks found, return a verdict now
        if peak['mode_bin']['ibin'] == 0:
            peak['assessment'] = { 'is_found': False, 'fraction': 0.0 }
            return(peak)

        #### Grow the extent to get most of the signal
        ibin = peak['mode_bin']['ibin']
        iextent = 1
        done = 0
        prev_n_spectra = peak['mode_bin']['intensity']
        peak['extended'] = { 'extent': 0, 'intensity': float(peak['mode_bin']['intensity']) }
        #print(spec[composite_type]['n_peaks'][ibin-10:ibin+10])
        while not done:
            prev_intensity = peak['extended']['intensity']
            new_intensity = prev_intensity + y[ibin-iextent] + y[ibin+iextent]
            peak['extended'] = { 'extent': iextent, 'intensity': float(new_intensity) }
            if new_intensity/prev_intensity < 1.02:
                #print('Extent reached by intensity')
                done = 1
            if iextent > 15:
                #print('Extent reached by max extent')
                done = 1

            # Increment and check for reaching the edge of the window
            iextent += 1
            if ibin-iextent < 0 or ibin+iextent >= len(x):
                done = 1


        if peak['extended']['intensity'] > 100:
            extent = peak['extended']['extent'] * 2
            min_bin = int(ibin - extent)
            original_min_bin = min_bin
            if min_bin < 0:
                min_bin = 0
            max_bin = int(ibin + extent)
            if max_bin >= len(x):
                max_bin = len(x) - 1
            window_x = x[min_bin:max_bin]
            window_y = y[min_bin:max_bin]

            center = original_min_bin - min_bin + extent
            binsize = x[1] - x[0]
            #print([original_min_bin,min_bin,max_bin,original_min_bin + extent,extent,len(window_x)])
            if show_interactive_plots:
                plt.step(window_x,window_y)
                plt.plot(window_x-binsize/2,gaussian_function(window_x,window_y[center],window_x[center],binsize*3),'ro:')
                plt.show()

            try:
                popt,pcov = curve_fit(gaussian_function,window_x,window_y,p0=[window_y[center],window_x[center],binsize*3])
                peak['fit'] = { 'mz': popt[1], 'intensity': popt[0], 'sigma': popt[2], 'delta_mz': popt[1]-reference_x }
                peak['assessment'] = { 'is_found': True, 'fraction': 0.0, 'comment': 'Peak found and fit' }
            except:
                peak['assessment'] = { 'is_found': False, 'fraction': 0.0, 'comment': 'Gaussian fit failed to converge' }

            if show_interactive_plots:
                plt.step(window_x,window_y)
                plt.plot(window_x-binsize/2,gaussian_function(window_x,*popt),'ro:')
                title = "Fit results: mz = %.4f, sigma = %.4f" % (popt[1], popt[2])
                print(peak)
                plt.title(title)
                plt.show()

        else:
            peak['assessment'] = { 'is_found': False, 'fraction': 0.0, 'comment': 'Too few peaks found in ROI' }
            #print(peak)

        return(peak)


    ####################################################################################################
    #### Assess Regions of Interest to make calls
    def assess_ROIs(self):
        results = { 'labeling': { 'scores': { 'TMT': 0, 'TMT6': 0, 'TMT10': 0, 'iTRAQ': 0, 'iTRAQ4': 0, 'iTRAQ8': 0 } } }

        #### Determine what the denominator of MS2 spectra should be
        n_ms2_spectra = self.metadata['files'][self.mzml_file]['spectra_stats']['n_ms2_spectra']
        if ( 'n_HCD_spectra' in self.metadata['files'][self.mzml_file]['spectra_stats'] and
                self.metadata['files'][self.mzml_file]['spectra_stats']['n_HCD_spectra'] > 0 and 
                self.metadata['files'][self.mzml_file]['spectra_stats']['n_HCD_spectra'] <
                self.metadata['files'][self.mzml_file]['spectra_stats']['n_ms2_spectra'] ):
            n_ms2_spectra = self.metadata['files'][self.mzml_file]['spectra_stats']['n_HCD_spectra']

        #### Determine the fragmentation_type
        results['fragmentation_type'] = self.metadata['files'][self.mzml_file]['spectra_stats']['fragmentation_type']
        self.metadata['files'][self.mzml_file]['summary'] = results

        #### If no ROIs were computed, cannot continue
        if 'ROIs' not in self.metadata['files'][self.mzml_file]:
            return

        for ROI in self.metadata['files'][self.mzml_file]['ROIs']:
            ROIobj = self.metadata['files'][self.mzml_file]['ROIs'][ROI]
            if ROIobj['peak']['assessment']['is_found']:
                type = ROIobj['type']
                n_spectra = ROIobj['peak']['extended']['n_spectra']
                ROIobj['peak']['assessment']['fraction'] = n_spectra / n_ms2_spectra
                results['labeling']['scores'][type] += ROIobj['peak']['assessment']['fraction']

        #### Make the call for TMT or iTRAQ
        if results['labeling']['scores']['TMT'] > results['labeling']['scores']['iTRAQ']:
            if results['labeling']['scores']['TMT'] > 2:
                results['labeling']['call'] = 'TMT'
            elif results['labeling']['scores']['TMT'] < 1:
                results['labeling']['call'] = 'none'
            else:
                results['labeling']['call'] = 'ambiguous'
        else:
            if results['labeling']['scores']['iTRAQ'] > 1:
                if results['labeling']['scores']['iTRAQ4'] > results['labeling']['scores']['iTRAQ8']:
                    results['labeling']['call'] = 'iTRAQ4'
                else:
                    results['labeling']['call'] = 'iTRAQ8'
            elif results['labeling']['scores']['iTRAQ'] < 0.5:
                results['labeling']['call'] = 'none'
            else:
                results['labeling']['call'] = 'ambiguous'



    ####################################################################################################
    #### Log an event
    def log_event(self, status, code, message):

        category = 'UNKNOWN'
        if status == 'WARNING':
            category = 'warnings'
        elif status == 'ERROR':
            category = 'errors'
        else:
            eprint(f"FATAL ERROR: Unrecognized event status '{status}'")
            eprint(sys.exc_info())
            sys.exit()

        #### Record the event
        full_message = f"{status}: [{code}]: {message}"
        self.metadata['problems'][category]['count'] += 1
        self.metadata['problems'][category]['list'].append(full_message)
        if code not in self.metadata['problems'][category]['codes']:
            self.metadata['problems'][category]['codes'][code] = 1
        else:
            self.metadata['problems'][category]['codes'][code] += 1

        #### If this is an error, also update the overall state
        if status == 'ERROR':
            self.metadata['state']['status'] = status
            self.metadata['state']['code'] = code
            self.metadata['state']['message'] = message
            if self.verbose >= 1:
                eprint(full_message)


####################################################################################################
#### Gaussian function used for curve fitting
def gaussian_function(x,a,x0,sigma):
    return a * numpy.exp( -(x-x0)**2 / (2*sigma**2) )


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Creates an index for an MSP spectral library file')
    argparser.add_argument('--use_cache', action='store', help='Set to true to use the cached file instead of regenerating')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('files', type=str, nargs='+', help='Filenames of one or more mzML files to read')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None: verbose = 1

    #### Loop over all the files to ensure that they are really there before starting work
    for file in params.files:
        if not os.path.isfile(file):
            print(f"ERROR: File '{file}' not found or not a file")
            return

    #### Loop over all the files processing them
    for file in params.files:

        #### Process the mzML file
        processor = SpectrumProcessor(file, verbose=verbose)
        processor.user_parameters['first_pass_tolerance'] = 20
        processor.read_header()
        if not params.use_cache:
            processor.read_spectra()

#### For command line usage
if __name__ == "__main__": main()
