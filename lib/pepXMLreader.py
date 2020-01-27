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
import gzip

#### Import technical modules and pyteomics
from pyteomics import pepxml, auxiliary


####################################################################################################
#### PepXml Reader class
class PepXmlReader:


    ####################################################################################################
    #### Constructor
    def __init__(self, pepxml_file, verbose=None):
        self.pepxml_file = pepxml_file

        #### Store the provided metadata or create a template to work in
        if metadata is None:
            self.metadata = { 'files': { mzml_file: {} } }
        else:
            self.metadata = metadata
            metadata['files'][mzml_file] = {}

        #### Create a place to store our composite spectra for analysis
        self.psms = []

        #### Set verbosity
        if verbose is None: verbose = 0
        self.verbose = verbose


    ####################################################################################################
    #### Read psms
    def read_psms(self):

        #### Set up information
        t0 = timeit.default_timer()
        stats = { 'n_psms': 0 }

        #### Show information
        if self.verbose >= 1:
            eprint(f"INFO: Reading pepXML file {self.pepxml_file}")
            progress_intro = False

        #### If the pepXML is gzipped, then open with zlib, else a plain open
        match = re.search('\.gz$',self.pepxml_file)
        if match:
            infile = gzip.open(self.pepxml_file)
        else:
            infile = open(self.pepxml_file, 'rb')

        #### Read psms from the file
        with pepxml.read(infile) as reader:
            for psm in reader:

                #### Testing. Print the data structure of the first spectrum
                if stats['n_psms'] == 0:
                    auxiliary.print_tree(psm)
                    sys.exit(10)

                #### Update counters and print progress
                stats['n_psms'] += 1
                if self.verbose >= 1:
                    if stats['n_psms']/1000 == int(stats['n_psms']/1000):
                        if not progress_intro:
                            eprint("INFO: Reading psms.. ", end='')
                            progress_intro = True
                        eprint(f"{stats['n_psms']}.. ", end='', flush=True)

        infile.close()
        if self.verbose >= 1: eprint("")

        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"INFO: Read {stats['n_psms']} psms from {self.mzml_file}")
        print(f"INFO: Elapsed time: {t1-t0}")
        print(f"INFO: Processed {stats['n_psms']/(t1-t0)} psms per second")


    ####################################################################################################
    #### Read header
    #### Open the mzML file and read line by line into a text buffer that we will XML parse
    def read_header(self):

        #### If the mzML is gzipped, then open with zlib, else a plain open
        match = re.search('\.gz$',self.mzml_file)
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
        namespace = namespace[None]

        #### Create a reference of instruments we know
        instrument_by_category = { 'pureHCD': [ 'MS:1001911|Q Exactive' ],
            'variable': [ 'MS:1001910|LTQ Orbitrap Elite', 'MS:1001742|LTQ Orbitrap Velos', 'MS:1000556|LTQ Orbitrap XL',
                'MS:1000555|LTQ Orbitrap Discovery' ] }

        #### Restructure it into a dict by PSI-MS identifier
        instrument_attributes = {}
        for instrument_category in instrument_by_category:
            for instrument_string in instrument_by_category[instrument_category]:
                accession,name = instrument_string.split('|')
                instrument_attributes[accession] = { 'category': instrument_category, 'accession': accession, 'name': name }

        #### Get all the CV params in the header and look for ones we know about
        cv_params = xmlroot.findall('.//{http://psi.hupo.org/ms/mzml}cvParam')
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
            print("ERROR: Did not recognize the instrument. Please teach me about this instrument.")


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
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Reads a pepXML file')
    argparser.add_argument('--min_psm_probability', action='store', default=0.90, help='Minimum PSM (PeptideProphet or iProphet) probability to accept')
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

        #### Read the pepXML file
        reader = PepXmlReader(file, verbose=verbose)
        #reader.read_header()
        reader.read_psms()


#### For command line usage
if __name__ == "__main__": main()
