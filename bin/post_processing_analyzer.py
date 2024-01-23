#!/usr/bin/env python3

import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

import os
import argparse
import multiprocessing
import copy
import csv

sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/../lib")
from simple_pepxml_reader import SimplePepXmlReader
from qualscore_mzML_reader import QualscoreMzMlReader


####################################################################################################
#### Read spectra
def read_mzMLs(search_results, mzMLs_path=None, verbose=None):

    #### Show progress information
    if verbose is not None and verbose >= 0:
        eprint(f"INFO: Reading spectra from mzMLs")
        progress_intro = False

    spectrum_data = []

    for mzml_root_name in search_results.mzmls:

        #### Try to find the mzML file
        paths = [ '.' ]
        if mzMLs_path is not None:
            paths.append(mzMLs_path)
        found_file = False
        tried_paths = []
        for path in paths:
            for suffix in [ '.mzML', '.mzML.gz' ]:
                mzml_filename = path + '/' + mzml_root_name + suffix
                tried_paths.append(mzml_filename)
                if os.path.isfile(mzml_filename):
                    found_file = True
                    break
        if not found_file:
            print(f"ERROR: Unable to find mzML file to load spectra from. Tried these:")
            for path in tried_paths:
                print(f"  - {path}")
            return

        #### Read the spectra
        reader = QualscoreMzMlReader(mzml_filename, msrun_name=mzml_root_name, verbose=verbose)
        reader.read_header()
        reader.read_spectra()
        spectrum_list = reader.spectra

        print(f"INFO: spectrum_list has {len(spectrum_list)} items")
        spectrum_data.extend(spectrum_list)

    return spectrum_data


####################################################################################################
#### Read spectra
def prepare_spectrum_data(spectrum_data, search_results, verbose):

    #### Show progress information
    if verbose >= 0:
        eprint(f"INFO: Computing spectrum metrics")
        progress_intro = False

    pxd_identifier = 'PXD999001'
    merged_data = []

    #### Index the search_results
    counter = 0
    search_results_index = {}
    print(f"INFO: Processing {len(search_results.psms)} PSMs")
    for psm in search_results.psms:
        msrun_name = psm[0]
        scan_number = int(psm[1])
        if msrun_name not in search_results_index:
            search_results_index[msrun_name] = {}
        search_results_index[msrun_name][scan_number] = counter
        counter += 1

    print(f"INFO: Processing {len(spectrum_data)} raw spectra")
    for raw_spectrum in spectrum_data:
        msrun_name = raw_spectrum['msrun_name']
        scan_number = int(raw_spectrum['scan_number'])
        charge_state = int(raw_spectrum['charge_state'])
        precursor_mz = float(raw_spectrum['precursor_mz'])
        minimum_intensity = float(raw_spectrum['minimum intensity'])
        weighted_snr = float(raw_spectrum['weighted_snr'])
        reporter_ions_median_snr = float(raw_spectrum['reporter_ions_median_snr'])
        quality_score = float(raw_spectrum['quality_score'])
        probability = 0
        peptidoform = 0
        if scan_number in search_results_index[msrun_name]:
            index = search_results_index[msrun_name][scan_number]
            probability = float(search_results.psms[index][3])
            peptidoform = search_results.psms[index][4]
        n_peaks = len(raw_spectrum['mzs'])

        usi = f"mzspec:{pxd_identifier}:{msrun_name}:scan:{scan_number}:{peptidoform}/{charge_state}"

        row = [ msrun_name, scan_number, charge_state, precursor_mz, probability , peptidoform, n_peaks, minimum_intensity, weighted_snr, reporter_ions_median_snr, quality_score, usi ]
        merged_data.append(row)

    columns = [ 'msrun_name', 'scan', 'charge', 'precursor_mz', 'probability', 'peptidoform', 'n_peaks', 'minimum_intensity', 'weighted_snr', 'reporter_ions_median_snr', 'quality_score', 'usi' ]
    with open('qualscore_metrics.tsv', 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        outfile.write("\t".join(columns) + "\n")
        writer.writerows(merged_data)

    return merged_data


####################################################################################################
#### Main function for command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Computes the quality and identifyability of all MS2 spectra')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--pepxml_file', action='store', default=None, help='Filename of the PeptideProphet or iProphet PepXML file from first pass processing')
    argparser.add_argument('--mzMLs_path', action='store', default=None, help='Path of the mzML files if they will not be found in the expected locations')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    #### Ensure that the PepXML file is really there before starting work
    input_file = params.pepxml_file
    if not os.path.isfile(input_file):
        print(f"ERROR: File '{input_file}' not found or not a file")
        return

    #### Load the relevant search result information from the pepXML file
    search_results = SimplePepXmlReader()
    search_results.read(input_file, verbose=verbose)

    spectrum_data = read_mzMLs(search_results, mzMLs_path=params.mzMLs_path, verbose=verbose)

    prepare_spectrum_data(spectrum_data, search_results, verbose=verbose)



if __name__ == "__main__": main()
