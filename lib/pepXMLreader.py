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
    def __init__(self, pepxml_file, reference_file=None, verbose=None):

        self.pepxml_file = pepxml_file
        self.reference_file = reference_file
        self.reference_peptides = {}

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

        #### Print a header
        print("\t".join(['scan', 'pool', 'PepProProb', 'iProProb', 'charge', 'PTMProProb', 'isSeqInDataset', 'isSeqInPool', 'IsPepformInPool', 'SameNoPhosphos', 'HasAPhospho', 
            'sequence', 'RefPepform','CalcPepform', 'PTMProProbsString', 'USI' ]))

        #### Read psms from the file
        with pepxml.read(infile) as reader:
            for psm in reader:

                peptideprophet_probability = None
                iprophet_probability = None
                keep = False
                mean_best_probability = -1
                peptide_str = 'xx'
                pool = '?'
                msrun_name = '?'

                #print(psm)
                sequence = psm['search_hit'][0]['peptide']
                charge = psm['assumed_charge']
                spectrum_name = psm['spectrum']
                match = re.search("_(pool\d)_",spectrum_name)
                if match:
                    pool = match.group(1)
                match = re.match("(.+)\.\d+\.\d+\.\d+$",spectrum_name)
                if match:
                    msrun_name = match.group(1)

                for analysis_result in psm['search_hit'][0]['analysis_result']:
                    if analysis_result['analysis'] == 'peptideprophet':
                        peptideprophet_probability = analysis_result['peptideprophet_result']['probability']
                    if analysis_result['analysis'] == 'interprophet':
                        iprophet_probability = analysis_result['interprophet_result']['probability']
                    if analysis_result['analysis'] == 'ptmprophet':
                        #print(analysis_result)
                        #print(analysis_result['ptmprophet_result']['ptm'][0:3])
                        if analysis_result['ptmprophet_result']['ptm'][0:3] == 'STY':
                            peptide_str = analysis_result['ptmprophet_result']['ptm_peptide']
                            mean_best_probability = analysis_result['ptmprophet_result']['parameter']['mean_best_prob']
                if iprophet_probability is not None and iprophet_probability >= 0.90:
                    keep = True
                if mean_best_probability < 0.9:
                    keep = False

                #### Generate a peptidoform in proper notation
                peptidoform = '??????'
                phospho_peptidoform = '??????'
                has_alanine_phospho = 'N'
                n_phosphos = 0
                #print(psm['search_hit'][0])
                if 'modifications' in psm['search_hit'][0]:
                    #print(psm['search_hit'][0]['modifications'])
                    residues = list(sequence)
                    phospho_residues = list(sequence)
                    nterm = ''
                    for modification in psm['search_hit'][0]['modifications']:
                        offset = modification['position']
#                        if 'variable' in modification:
#                            if abs( modification['variable'] - 79.966 ) < 0.01:
#                                residues[offset-1] += '[Phospho]'
#                            elif abs( modification['variable'] - 15.9949 ) < 0.01:
#                                residues[offset-1] += '[Hydroxylation]'
#                            elif abs( modification['variable'] - 0.984 ) < 0.01:
#                                residues[offset-1] += '[Deamidation]'
#                            elif abs( modification['variable'] - (-17.026) ) < 0.01:
#                                residues[offset-1] += '[Pyro-glu]'
#                            elif abs( modification['variable'] - (-18.010) ) < 0.01:
#                                residues[offset-1] += '[Pyro_glu]'
#                            else:
#                                print(f"ERROR: Unable to translate {modification}")
                        if 'mass' in modification:

                            #### Phospho-only peptidoforms
                            if abs( modification['mass'] - 181.01401 ) < 0.01:
                                phospho_residues[offset-1] += '[Phospho]'
                                n_phosphos += 1
                            elif abs( modification['mass'] - 243.0297 ) < 0.01:
                                phospho_residues[offset-1] += '[Phospho]'
                                n_phosphos += 1
                            elif abs( modification['mass'] - 166.998359 ) < 0.01:
                                phospho_residues[offset-1] += '[Phospho]'
                                n_phosphos += 1
                            elif abs( modification['mass'] - 151.003445 ) < 0.01:
                                phospho_residues[offset-1] += '[Phospho]'
                                n_phosphos += 1
                                has_alanine_phospho = 'Y'

                            #### All-mod peptidforms
                            if abs( modification['mass'] - 160.030649 ) < 0.01:
                                #residues[offset-1] += '[Carbamidomethyl]'
                                pass
                            elif abs( modification['mass'] - 181.01401 ) < 0.01:
                                residues[offset-1] += '[Phospho]'
                            elif abs( modification['mass'] - 243.0297 ) < 0.01:
                                residues[offset-1] += '[Phospho]'
                            elif abs( modification['mass'] - 166.998359 ) < 0.01:
                                residues[offset-1] += '[Phospho]'
                            elif abs( modification['mass'] - 151.003445 ) < 0.01:
                                residues[offset-1] += '[Phospho]'
                            elif abs( modification['mass'] - 147.0354 ) < 0.01:
                                residues[offset-1] += '[Hydroxylation]'
                            elif abs( modification['mass'] - 202.074213 ) < 0.01:
                                residues[offset-1] += '[Hydroxylation]'
                            elif abs( modification['mass'] - 113.047664 ) < 0.01:
                                residues[offset-1] += '[Hydroxylation]'
                            elif abs( modification['mass'] - 115.026943 ) < 0.01:
                                residues[offset-1] += '[Deamidation]'
                            elif abs( modification['mass'] - 129.042594 ) < 0.01:
                                residues[offset-1] += '[Deamidation]'
                            elif abs( modification['mass'] - 111.032029 ) < 0.01:
                                residues[offset-1] += '[Pyro-glu]'
                            elif abs( modification['mass'] - 111.032029 ) < 0.01:
                                residues[offset-1] += '[Pyro_glu]'
                            elif abs( modification['mass'] - 143.0041 ) < 0.01:
                                #residues[offset-1] += '[Carbamidomethyl]'
                                residues[offset-1] += '[Pyro_glu]'
                            elif abs( modification['mass'] - 43.018425 ) < 0.01:
                                nterm = '[Acetyl]-'
                            else:
                                print(f"ERROR: Unable to translate {modification}")
                    peptidoform = nterm + ''.join(residues)
                    phospho_peptidoform = ''.join(phospho_residues)

                #print(peptidoform)
                is_sequence_in_dataset = 'N'
                if sequence in self.reference_peptides['by_sequence']:
                    is_sequence_in_dataset = 'Y'

                is_sequence_in_pool = 'N'
                pool_sequence = f"{pool}-{sequence}"
                if pool_sequence in self.reference_peptides['by_pool_sequence']:
                    is_sequence_in_pool = 'Y'

                same_number_of_phosmods = 'N'
                if pool_sequence in self.reference_peptides['by_pool_sequence']:
                    ref_n_mods = int(self.reference_peptides['by_pool_sequence'][pool_sequence][0]['n_mods'])
                    #eprint(f"{n_phosphos}, {ref_n_mods}, {type(n_phosphos)}, {type(ref_n_mods)}")
                    same_number_of_phosmods = f"{n_phosphos},{ref_n_mods}"
                    if n_phosphos == ref_n_mods:
                        same_number_of_phosmods = 'Y'

                is_peptidoform_in_pool = 'N'
                pool_peptidoform = f"{pool}-{phospho_peptidoform}"
                if pool_peptidoform in self.reference_peptides['by_pool_peptidoform']:
                    is_peptidoform_in_pool = 'Y'

                #### Find the reference peptide
                reference_peptidoform = '---------------------'
                if is_sequence_in_pool == 'Y':
                    reference_peptidoform = self.reference_peptides['by_pool_sequence'][pool_sequence][0]['peptidoform']

                if keep:
                    usi = f"mzspec:PXD007058:{msrun_name}:scan:{psm['start_scan']}:{peptidoform}/{charge}"
                    row = [ str(psm['start_scan']), pool, str(peptideprophet_probability), str(iprophet_probability), 
                        str(charge), str(mean_best_probability), is_sequence_in_dataset, is_sequence_in_pool, is_peptidoform_in_pool, same_number_of_phosmods, has_alanine_phospho, 
                        sequence, reference_peptidoform, peptidoform, peptide_str, usi ]
                    print("\t".join(row))

                #### Testing. Print the data structure of the first spectrum
                #if stats['n_psms'] >1000:
                    #auxiliary.print_tree(psm)
                    #sys.exit(10)

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
        print(f"INFO: Read {stats['n_psms']} psms from {self.pepxml_file}")
        print(f"INFO: Elapsed time: {t1-t0}")
        print(f"INFO: Processed {stats['n_psms']/(t1-t0)} psms per second")


    ####################################################################################################
    #### Read reference file
    def read_reference_file(self):

        #### Show information
        if self.verbose >= 1:
            eprint(f"INFO: Reading reference file {self.reference_file}")

        iline = 0
        reference_peptides = { 'by_sequence': {}, 'by_pool_sequence': {}, 'by_pool_peptidoform': {}, 'by_pool_sequence_mods': {} }
        self.reference_peptides = reference_peptides

        #### Read psms from the file
        with open(self.reference_file) as infile:
            for line in infile:
                line = line.rstrip()
                columns = line.split(",")
                if iline == 0:
                    if columns[0] != 'pool':
                        print(f"ERROR: Expected first column of head to be 'pool', but instead found this first line:\n{line}")
                        sys.exit(6)
                    iline += 1
                    continue

                pool = columns[0]
                sequence = columns[1]
                mods = columns[5]

                pool_sequence = f"{pool}-{sequence}"
                pool_sequence_mods = f"{pool}-{sequence}-{mods}"
                offsets = mods.split("&")

                #### Create the peptidoform
                residues = list(sequence)
                for offset in offsets:
                    ioffset = int(offset)
                    residues[ioffset-1] += '[Phospho]'
                peptidoform = ''.join(residues)
                pool_peptidoform = f"{pool}-{peptidoform}"

                #print(f"{pool}, {sequence}, {mods}, {offsets}, {peptidoform}")
                peptide = { 'pool': pool, 'sequence': sequence, 'n_mods': columns[4], 'offsets': offsets, 'peptidoform': peptidoform }

                #### Add this entry to the by_sequence lookup
                if sequence not in reference_peptides['by_sequence']:
                    reference_peptides['by_sequence'][sequence] = []
                reference_peptides['by_sequence'][sequence].append(peptide)

                #### Add this entry to the by_pool_sequence lookup
                if pool_sequence not in reference_peptides['by_pool_sequence']:
                    reference_peptides['by_pool_sequence'][pool_sequence] = []
                reference_peptides['by_pool_sequence'][pool_sequence].append(peptide)

                #### Add this entry to the by_pool_peptidoform lookup
                if pool_peptidoform not in reference_peptides['by_pool_peptidoform']:
                    reference_peptides['by_pool_peptidoform'][pool_peptidoform] = []
                reference_peptides['by_pool_peptidoform'][pool_peptidoform].append(peptide)

                #### Add this entry to the by_pool_sequence_mods lookup
                if pool_sequence_mods not in reference_peptides['by_pool_sequence_mods']:
                    reference_peptides['by_pool_sequence_mods'][pool_sequence_mods] = []
                reference_peptides['by_pool_sequence_mods'][pool_sequence_mods].append(peptide)

                iline += 1

                #### Exit for debugging
                #if iline > 10:
                #    sys.exit()


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
    argparser.add_argument('--reference_file', action='store', default=None, help='Ferries et al peptide reference tsv file')
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

    reader = PepXmlReader(file, params.reference_file, verbose=verbose)

    #### If a reference file was given, read it
    if params.reference_file is not None:
        reader.read_reference_file()

    #### Loop over all the files processing them
    for file in params.files:

        #### Read the pepXML file
        reader.pepxml_file = file
        reader.read_psms()


#### For command line usage
if __name__ == "__main__": main()
