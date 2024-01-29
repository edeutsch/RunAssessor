#!/usr/bin/env python3

import sys

from pytz import NonExistentTimeError
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

import os
import argparse
import os.path
import timeit
import re
import numpy
import gzip
from lxml import etree
from multiprocessing.pool import ThreadPool

#### Import technical modules and pyteomics
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy import exp
from pyteomics import mzml, auxiliary

#### Import the metadata handler
from metadata_handler import MetadataHandler
from spectrum import Spectrum


####################################################################################################
#### mzML Assessor class
class QualscoreMzMlReader:


    ####################################################################################################
    #### Constructor
    def __init__(self, mzml_file=None, msrun_name=None, verbose=None):

        self.mzml_file = mzml_file
        self.msrun_name = msrun_name

        #### Create a place to store the spectra for analysis
        self.spectra = []

        #### Create a storage area for referenceable param groups, which are not really handled by pyteomics
        self.referenceable_param_group_list = {}

        #### Set verbosity
        if verbose is None:
            verbose = 0
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
            'n_length_zero_spectra': 0,
            'n_length_lt10_spectra': 0,
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
            'n_HR_QTOF_spectra': 0,
            'high_accuracy_precursors': 'unknown',
            'fragmentation_type': 'unknown',
            'fragmentation_tag': 'unknown'
        }

        #### Show information
        if self.verbose >= 1:
            eprint(f"INFO: Reading mzML file {self.mzml_file}")
            progress_intro = False

        #### If the mzML is gzipped, then open with zlib, else a plain open
        if self.mzml_file.endswith('.gz'):
            infile = gzip.open(self.mzml_file)
        else:
            infile = open(self.mzml_file, 'rb')

        #### Try storing all spectra for multithreaded processing
        #spectra = []

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

                    #### Extract the spectrum id and scan number
                    id = spectrum['id']
                    match = re.search(r'scan=(\d+)', id)
                    if match:
                        scan_number = match.group(1)
                    else:
                        scan_number = 0

                    #### Extract the MS level of the spectrum
                    try:
                        ms_level = spectrum['ms level']
                    #### If no MS level is found, see if it can be extracted from a referenceable param group
                    except:
                        ms_level = 0
                        if 'ref' in spectrum:
                            group_id = spectrum['ref']
                            if group_id in self.referenceable_param_group_list and 'ms level' in self.referenceable_param_group_list[group_id]:
                                ms_level = int(self.referenceable_param_group_list[group_id]['ms level'])

                    #### Look for a filter string and parse it
                    if 'filter string' in spectrum['scanList']['scan'][0]:
                        filter_string = spectrum['scanList']['scan'][0]['filter string']
                        self.parse_filter_string(filter_string,stats)

                    #### There's only a filter string for Thermo data, so for others, record a subset of information
                    else:
                        stats[f"n_ms{ms_level}_spectra"] += 1
                        if self.metadata['files'][self.mzml_file]['instrument_model']['category'] == 'QTOF':
                            stats['high_accuracy_precursors'] = 'true'
                            stats['fragmentation_type'] = 'HR_QTOF'
                            stats['fragmentation_tag'] = 'HR QTOF'
                            if ms_level > 1:
                                stats['n_HR_QTOF_spectra'] += 1

                    #### If the ms level is greater than 2, fail
                    if ms_level > 4:
                        self.log_event('ERROR','MSnTooHigh',f"MS level is greater than we can handle at '{spectrum['ms level']}'")
                        break

                    #### If the ms level is 2, then examine it for information
                    if ms_level == 2 and 'm/z array' in spectrum:

                        precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                        #### Try to get the charge information
                        try:
                            charge_state = int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
                        except:
                            charge_state = 'unknown'
                        charge_stat = f"n_charge_{charge_state}_precursors"
                        if charge_stat not in stats:
                            stats[charge_stat] = 0
                        stats[charge_stat] += 1

                        #### Check for zero length and very sparse spectra
                        if len(spectrum['m/z array']) == 0:
                            stats['n_length_zero_spectra'] += 1
                        if len(spectrum['m/z array']) < 10:
                            stats['n_length_lt10_spectra'] += 1

                        #### Transform the mzML spectrum into a standardized Spectrum object
                        new_spectrum_obj = Spectrum()
                        new_spectrum_obj.fill(
                            mzs=spectrum['m/z array'],
                            intensities=spectrum['intensity array'],
                            precursor_mz=precursor_mz,
                            charge_state=charge_state
                        )
                        new_spectrum_obj.compute_spectrum_metrics()
                        #self.spectra.append(new_spectrum_obj)

                        new_spectrum_dict = {
                            'msrun_name': self.msrun_name,
                            'scan_number': scan_number,
                            'charge_state': charge_state,
                            'precursor_mz': precursor_mz,
                            'mzs': spectrum['m/z array'],
                            'intensities': spectrum['intensity array'],
                            'minimum intensity': new_spectrum_obj.attributes['minimum intensity'],
                            'weighted_snr': new_spectrum_obj.metrics['weighted_snr'],
                            'quality_score': new_spectrum_obj.metrics['quality_score'],
                        }
                        if 'reporter_ions' in new_spectrum_obj.metrics['stats_per_bin']:
                            new_spectrum_dict['reporter_ions_median_snr'] = new_spectrum_obj.metrics['stats_per_bin']['reporter_ions']['median_snr']
                            new_spectrum_dict['reporter_ion_intensities'] = new_spectrum_obj.metrics['stats_per_bin']['reporter_ions']['intensities']
                        else:
                            new_spectrum_dict['reporter_ions_median_snr'] = 0.0

                        self.spectra.append(new_spectrum_dict)


                    #### Update counters and print progress
                    stats['n_spectra'] += 1
                    if self.verbose >= 1:
                        if stats['n_spectra']/5000 == int(stats['n_spectra']/5000):
                            if not progress_intro:
                                eprint("INFO: Reading spectra.. ", end='')
                                progress_intro = True
                            eprint(f"{stats['n_spectra']}.. ", end='', flush=True)
                            #eprint(".", end='', flush=True)
            #except:
            else:
                self.log_event('ERROR','MzMLCorrupt',f"Pyteomics threw an error reading mzML file! File may be corrupt. Check file '{self.mzml_file}'")

        infile.close()
        if self.verbose >= 1: eprint("")

        #### If there are no spectra that were saved, then we're done with this file
        if stats['n_ms2_spectra'] == 0:
            self.log_event('ERROR','NoMS2Scans',f"This mzML file has no MS2 scans. Check file '{self.mzml_file}'")


        #### Now that we've read anything, try to process the data multithreaded
        #n_threads = 8
        #if self.verbose >= 1: eprint(f"INFO: processing {len(spectra)} spectra multithreaded with {n_threads} threads")
        #pool = ThreadPool(n_threads)
        #results = pool.imap_unordered(self.add_spectrum, spectra)
        #results = pool.map(self.add_spectrum, spectra)
        #pool.close()
        #pool.join()
        #if self.verbose >= 1: eprint(f"INFO: Done.")

        #### Write the composite spectrum to a file
        #with open('zzcomposite_spectrum.pickle','wb') as outfile:
        #    pickle.dump(self.composite,outfile)

        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"INFO: Read {stats['n_spectra']} spectra from {self.mzml_file} in {int(t1-t0)} sec ({stats['n_spectra']/(t1-t0)} spectra per sec)")


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

        #### Keep a dict of random recognized things
        recognized_things = {}

        #### Read line by line
        for line in infile:

            if not isinstance(line, str):
                line = str(line, 'utf-8', 'ignore')

            #### Completely skip the <indexedmzML> tag if present
            if '<indexedmzML ' in line:
                continue

            #### Look for first tag after the header and end when found
            if '<run ' in line:
                break

            #### Look for some known artificts that are useful to record
            if 'softwareRef="tdf2mzml"' in line:
                recognized_things['converted with tdf2mzml'] = True
                recognized_things['likely a timsTOF'] = True

            #### Update the buffer and counter
            if counter > 0:
                buffer += line
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
                'MS:1000649|Exactive',
                'MS:1001911|Q Exactive',
                'MS:1002526|Exactive Plus',
                'MS:1002634|Q Exactive Plus',
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
                'MS:1003029|Orbitrap Eclipse',
                'MS:1000483|Thermo Fisher Scientific instrument model'
            ],
            'QTOF': [
                'MS:1000126|Waters instrument model',
                'MS:1000122|Bruker Daltonics instrument model',
                'MS:1000121|AB SCIEX instrument model'
            ]
        }

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
                model_data = {
                    'accession': accession,
                    'name': instrument_attributes[accession]['name'],
                    'category': instrument_attributes[accession]['category']
                }
                if 'likely a timsTOF' in recognized_things:
                    model_data['inferred_name'] = 'timsTOF'
                #self.metadata['files'][self.mzml_file]['instrument_model'] = model_data
                found_instrument = 1

        #### If none are an instrument we know about about, ask for help
        if not found_instrument:
            self.log_event('ERROR','UnrecogInstr',f"Did not recognize the instrument. Please teach me about this instrument.")
            model_data = {
                'accession': None,
                'name': 'unknown',
                'category': 'unknown'
            }
            #self.metadata['files'][self.mzml_file]['instrument_model'] = model_data
            return


        #### Read and store the referenceable param groups
        referenceable_param_group_list = xmlroot.findall(f'.//{namespace}referenceableParamGroup')
        for referenceable_param_group in referenceable_param_group_list:
            group_id = referenceable_param_group.get('id')
            self.referenceable_param_group_list[group_id] = {}
            for subelement in list(referenceable_param_group):
                name = subelement.get('name')
                value = subelement.get('value')
                self.referenceable_param_group_list[group_id][name] = value



####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Creates an index for an MSP spectral library file')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('files', type=str, nargs='+', help='Filenames of one or more mzML files to read')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    #### Loop over all the files to ensure that they are really there before starting work
    for file in params.files:
        if not os.path.isfile(file):
            print(f"ERROR: File '{file}' not found or not a file")
            return

    #### Loop over all the files processing them
    for file in params.files:

        #### Assess the mzML file
        reader = QualscoreMzMlReader(file, verbose=verbose)
        reader.read_header()
        reader.read_spectra()


#### For command line usage
if __name__ == "__main__": main()
