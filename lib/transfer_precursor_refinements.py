#!/usr/bin/env python3

import sys
import os
import argparse
import os.path
import timeit
import re
import json
import numpy
import pickle
import gzip

from pyteomics import mzml, auxiliary
from psims.transform.mzml import MzMLTransformer

def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


####################################################################################################
class CustomMzMLTransformer(MzMLTransformer):
    def __init__(self, input_stream, output_stream, transform=None, transform_description=None,
                 sort_by_scan_time=False, scan_information=None, log_file_handle=None):
        super().__init__(input_stream, output_stream, transform=transform, transform_description=transform_description,
                 sort_by_scan_time=sort_by_scan_time)

        self.scan_information = scan_information
        self.log_file_handle = log_file_handle

    def format_spectrum(self, spectrum):
        new_spectrum = super().format_spectrum(spectrum)
        ms_level = None
        for param in new_spectrum['params']:
            if param['name'] == 'MS:1000511':
                ms_level = param['value']
        if ms_level > 1:
            #print(new_spectrum)
            nativeId = spectrum['id']
            scan_number = -1
            match = re.search(r'scan=(\d+)',nativeId)
            if match:
                scan_number = int(match.group(1))
            else:
                eprint(f"ERROR: Unable to extract scan number from {nativeId}")

            #print(f"scan={scan_number}")
            previous_precursor_mz = new_spectrum['precursor_information'][0]['mz']
            previous_charge = new_spectrum['precursor_information'][0]['charge']
            new_precursor = self.scan_information[scan_number]['precursor_mz']
            new_charge = self.scan_information[scan_number]['charge_state']

            new_spectrum['precursor_information'][0]['mz'] = new_precursor
            new_spectrum['precursor_information'][0]['charge'] = new_charge

            print(f"{scan_number:6d}\t{previous_precursor_mz:10.4f}\t{previous_charge:2d}\t{new_precursor:10.4f}\t" +
                f"{new_charge:2d}\t{float(new_precursor)-float(previous_precursor_mz):10.4f}\t{int(new_charge)-int(previous_charge):2d}", file=self.log_file_handle )
            #print(self.scan_information[scan_number])

        return new_spectrum




####################################################################################################
#### Monocle CSV Reader class reads a Monocle CSV and extracts the precursor m/z values and charge values
class MonocleCsvReader:


    ####################################################################################################
    #### Constructor
    def __init__(self, csv_file, verbose=None):

        self.csv_file = csv_file

        #### General stats
        self.stats = {}

        #### Information about each scan
        self.scans = {}

        #### Set verbosity
        if verbose is None: verbose = 0
        self.verbose = verbose


    ####################################################################################################
    #### Read spectra
    def read_spectra(self):

        #### Set up information
        t0 = timeit.default_timer()
        stats = {
            'n_spectra': 0,
            'n_ms0_spectra': 0,
            'n_ms1_spectra': 0,
            'n_ms2_spectra': 0,
            'n_ms3_spectra': 0,
            'n_length_zero_spectra': 0
        }
        self.stats = stats

        #### Show information
        if self.verbose >= 1:
            eprint(f"INFO: Reading Monocle CSV file {self.mzml_file}")

        #### Read spectra from the file
        with open(self.csv_file) as infile:
            for line in infile:

                columns = line.strip().split(',')

                #### Extract the MS level of the spectrum
                scan_number = columns[0]
                ms_level = columns[1]
                ms_level_stat = f"n_ms{ms_level}_spectra"
                stats[ms_level_stat] += 1

                precursor_mz = float(columns[2])
                charge_state = int(columns[4])

                charge_stat = f"n_charge_{charge_state}_precursors"
                if charge_stat not in stats:
                    stats[charge_stat] = 0
                stats[charge_stat] += 1

                #### Store the supplied precursor_ms and charge state in the cache
                self.scans[scan_number] = {'precursor_mz': precursor_mz, 'charge_state': charge_state }

                #### Update counters
                stats['n_spectra'] += 1

        infile.close()

        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"INFO: Read {stats['n_spectra']} spectra from {self.csv_file} in {int(t1-t0)} sec ({stats['n_spectra']/(t1-t0)} spectra per sec)")




####################################################################################################
#### Monocle mzML Reader class reads an mzML and extracts the precursor m/z values and charge values
class MonocleMzmlReader:


    ####################################################################################################
    #### Constructor
    def __init__(self, mzml_file, verbose=None):

        self.mzml_file = mzml_file

        #### General stats
        self.stats = {}

        #### Information about each scan
        self.scans = {}

        #### Set verbosity
        if verbose is None: verbose = 0
        self.verbose = verbose


    ####################################################################################################
    #### Read spectra
    def read_spectra(self):

        #### Set up information
        t0 = timeit.default_timer()
        stats = {
            'n_spectra': 0,
            'n_ms0_spectra': 0,
            'n_ms1_spectra': 0,
            'n_ms2_spectra': 0,
            'n_ms3_spectra': 0,
            'n_length_zero_spectra': 0
        }
        self.stats = stats

        #### Show information
        if self.verbose >= 1:
            eprint(f"INFO: Reading Monocle mzML file {self.mzml_file}", end='', flush=True)
            progress_intro = False

        #### If the mzML is gzipped, then open with zlib, else a plain open
        if self.mzml_file.endswith('.gz'):
            infile = gzip.open(self.mzml_file)
        else:
            infile = open(self.mzml_file, 'rb')

        #### Read spectra from the file
        with mzml.read(infile) as reader:
            #try:
            if True:
                for spectrum in reader:

                    #### Debugging. Print the data structure of the first spectrum
                    #if stats['n_spectra'] == 0:
                    #    auxiliary.print_tree(spectrum)
                    #    print(spectrum)
                    #    return

                    #### Extract the MS level of the spectrum
                    ms_level = spectrum['ms level']
                    ms_level_stat = f"n_ms{ms_level}_spectra"
                    stats[ms_level_stat] += 1

                    #### Save the spectrum in our cache
                    nativeId = spectrum['id']
                    scan_number = -1
                    match = re.search(r'scan=(\d+)',nativeId)
                    if match:
                        scan_number = int(match.group(1))
                    else:
                        eprint(f"ERROR: Unable to extract scan number from {nativeId}")

                    #### If the ms level is 2, then examine it for information
                    if ms_level == 2 and 'm/z array' in spectrum:
                        #### Monocle outputs the wrong terms
                        #precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                        #charge_state = int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
                        precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['isolation window target m/z']
                        charge_state = int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['possible charge state'])

                        charge_stat = f"n_charge_{charge_state}_precursors"
                        if charge_stat not in stats:
                            stats[charge_stat] = 0
                        stats[charge_stat] += 1

                        #### Store the supplied precursor_ms and charge state in the cache
                        self.scans[scan_number] = {'precursor_mz': precursor_mz, 'charge_state': charge_state }

                        #### Check for zero length and very sparse spectra
                        if len(spectrum['m/z array']) == 0:
                            stats['n_length_zero_spectra'] += 1

                    #### Update counters and print progress
                    stats['n_spectra'] += 1
                    if self.verbose >= 1:
                        if stats['n_spectra']/500 == int(stats['n_spectra']/500):
                            if not progress_intro:
                                #eprint("INFO: Reading spectra.. ", end='')
                                progress_intro = True
                            #eprint(f"{stats['n_spectra']}.. ", end='', flush=True)
                            eprint(".", end='', flush=True)
            #except:
            else:
                self.log_event('ERROR','MzMLCorrupt',f"Pyteomics threw an error reading mzML file! File may be corrupt. Check file '{self.mzml_file}'")

        infile.close()
        if self.verbose >= 1:
            eprint("")

        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"INFO: Read {stats['n_spectra']} spectra from {self.mzml_file} in {int(t1-t0)} sec ({stats['n_spectra']/(t1-t0)} spectra per sec)")


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Read and write an mzML file')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--input_filename', type=str, action='store', required=True, help='Name of the input ThermoRawFileParser mzML file')
    argparser.add_argument('--monocle_filename', type=str, action='store', required=True, help='Name of the Monocle csv or mzML output from which to get precursor information')
    argparser.add_argument('--output_filename', type=str, action='store', required=True, help='Name of the output mzML based on the TRFP mzML but with Monocle precursor information')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    #### Check that the filenames are defined and exist
    if not os.path.isfile(params.input_filename):
        print(f"ERROR: --input_filename '{params.input_filename}' not found or not a file")
        return
    if not os.path.isfile(params.monocle_filename):
        print(f"ERROR: --monocle_filename '{params.monocle_filename}' not found or not a file")
        return

    #### Ensure the output is not the same as the input
    if params.output_filename == params.input_filename:
        print(f"ERROR: --output_filename '{params.input_filename}' may not be the same as the --input_filename")
        return
    if params.output_filename == params.monocle_filename:
        print(f"ERROR: --monocle_filename '{params.monocle_filename}' may not be the same as the --input_filename")
        return

    #### If the Monocle file is an mzML, read and extract information from that
    if params.monocle_filename.endswith('.mzML') or params.monocle_filename.endswith('.mzML.gz'):
        reader = MonocleMzmlReader(params.monocle_filename, verbose=verbose)
        reader.read_spectra()
    elif params.monocle_filename.endswith('.csv'):
        reader = MonocleCsvReader(params.monocle_filename, verbose=verbose)
        reader.read_spectra()
    else:
        print(f"ERROR: Expected Monocle file to have extensions .csv or .mzML or .mzML.gz but '{params.monocle_filename}' does not")
        return

    #### Open in the input file
    print(f"INFO: Opening mzML file '{params.input_filename}' for reading")
    if params.input_filename.endswith('.gz'):
        infile = gzip.open(params.input_filename)
    else:
        infile = open(params.input_filename, 'rb')

    #### Open the output file
    try:
        outfile = open(params.output_filename, 'wb')
    except:
        print(f"ERROR: --output_filename '{params.output_filename}' is not writable. Permissions problem?")
        return

    log_file_handle = open(params.output_filename+'.log', 'w')

    print(f"INFO: Opening {params.output_filename} for writing")
    CustomMzMLTransformer(infile, outfile, scan_information=reader.scans, log_file_handle=log_file_handle,
        transform_description='Transfer Monocle precursor refinements to well-formed mzML').write()


#### For command line usage
if __name__ == "__main__": main()
