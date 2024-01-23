#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import timeit
import re
import gzip

#### Import technical modules and pyteomics
from pyteomics import pepxml, auxiliary


mass_lookup_table = {
    'n': { 43.018425: 'Acetyl', 305.214971: 'TMTpro' },
    'C': { 143.0041: 'Pyro-Carbamidomethyl', 160.030649: 'Carbamidomethyl', 69.021464: 'Cys->Dha' },
    'H': { 153.053827: 'Oxidation' },
    'E': { 111.032029: 'Glu->pyro-Glu' },
    'K': { 432.302109: 'TMTpro' },
    'M': { 147.0354: 'Oxidation', 163.030314: 'Dioxidation' },
    'N': { 115.026943: 'Deamidated' },
    'Q': { 111.032029: 'Gln->pyro-Glu', 129.042594: 'Deamidated' },
    'R': { 460.308257: 'TMTpro' },
    'S': { 166.998359: 'Phospho', 391.239174: 'TMTpro' },
    'T': { 181.01401: 'Phospho', 405.254825: 'TMTpro'  },
    'Y': { 243.0297: 'Phospho' },
    'W': { 202.074228: 'Oxidation', 218.069142: 'Dioxidation', 190.074228: 'Trp->Kynurenin' },
}

def mass_lookup(total_mass, residue):

    if residue in mass_lookup_table:
        for mass, mod_name in mass_lookup_table[residue].items():
            if abs( total_mass - mass ) < 0.01:
                return(mod_name)
    return



####################################################################################################
#### PepXml Reader class
class SimplePepXmlReader:


    ####################################################################################################
    #### Constructor
    def __init__(self, pepxml_file=None, verbose=None):

        #### Set verbosity
        if verbose is None:
            verbose = 0
        self.verbose = verbose

        #### Store parameters
        self.pepxml_file = pepxml_file

        #### Create a place to store the input mzML files and spectrum stats
        self.mzmls = {}
        self.psms = []


    ####################################################################################################
    #### Read the minimal information from the PepXML file
    def read(self, filename=None, verbose=None):

        #### Set parameters
        if filename is not None:
            self.pepxml_file = filename
        if verbose is not None:
            self.verbose = verbose

        #### Show progress information
        if self.verbose >= 0:
            eprint(f"INFO: Reading pepXML file {self.pepxml_file}")

        self.read_psms()



    ####################################################################################################
    #### Read header
    #### Open the PepXML file and read line by line into a text buffer that we will XML parse
    def read_header(self):

        #### If the mzML is gzipped, then open with zlib, else a plain open
        if self.pepxml_file.endswith('.gz'):
            infile = gzip.open(self.pepxml_file)
        else:
            infile = open(self.pepxml_file)

        #### set up a text buffer to hold the mzML header
        buffer = ''
        counter = 0

        #### Read line by line
        for line in infile:

            if not isinstance(line, str):
                line = str(line, 'utf-8', 'ignore')

            #### Look for first tag after the header and end when found
            match = re.search(r'<spectrum_query ',line)
            if match:
                break
            if counter > 0: buffer += line
            counter += 1

        #### Close file and process the XML in the buffer
        infile.close()

        #### Finish the XML by closing the tags
        buffer += '  </msms_run_summary>\n'
        buffer += '</msms_pipeline_analysis>\n'

        #### Get the root of the XML and the namespace
        xmlroot = etree.fromstring(buffer)
        namespace = xmlroot.nsmap
        namespace = namespace[None]

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


    ####################################################################################################
    #### Read psms
    def read_psms(self):

        #### Set up information
        t0 = timeit.default_timer()
        stats = { 'n_psms': 0 }
        self.psms = []

        #### Show information
        if self.verbose >= 1:
            eprint(f"INFO: Reading PSMs in pepXML file {self.pepxml_file}")
            progress_intro = False

        #### If the pepXML is gzipped, then open with zlib, else a plain open
        if self.pepxml_file.endswith('.gz'):
            infile = gzip.open(self.pepxml_file)
        else:
            infile = open(self.pepxml_file, 'rb')

        #### Read psms from the file
        with pepxml.read(infile) as reader:
            for psm in reader:

                peptideprophet_probability = None
                keep = False
                msrun_name = '?'

                sequence = psm['search_hit'][0]['peptide']
                charge = psm['assumed_charge']
                peptide_length = len(sequence)

                spectrum_name = psm['spectrum']
                match = re.match(r"(.+)\.\d+\.\d+\.\d+$",spectrum_name)
                if match:
                    msrun_name = match.group(1)

                for analysis_result in psm['search_hit'][0]['analysis_result']:
                    if analysis_result['analysis'] == 'peptideprophet':
                        peptideprophet_probability = analysis_result['peptideprophet_result']['probability']

                #### Generate a peptidoform in proper notation
                peptidoform = '??????'
                if 'modifications' in psm['search_hit'][0]:
                    residues = list(sequence)
                    nterm = ''

                    for modification in psm['search_hit'][0]['modifications']:
                        offset = modification['position']
                        total_mass = modification['mass']
                        if offset == 0:
                            mod_name = mass_lookup(total_mass, 'n')
                            if mod_name is not None:
                                nterm = f"[{mod_name}]-"
                            else:
                                print(f"ERROR: Unable to lookup nterm with mass {total_mass}")
                                exit()
                        elif 'mass' in modification:
                            residue = residues[offset-1]
                            mod_name = mass_lookup(total_mass, residue)
                            if mod_name is not None:
                                residues[offset-1] += f"[{mod_name}]"
                            else:
                                print(f"ERROR: Unable to lookup residue {residue} with total mass {total_mass}")
                                exit()

                    peptidoform = nterm + ''.join(residues)

                keep = True

                if keep:
                    row = [ msrun_name, str(psm['start_scan']), str(charge), str(peptideprophet_probability), peptidoform, peptide_length ]
                    self.psms.append(row)
                    self.mzmls[msrun_name] = True

                #### Testing. Print the data structure of a PSM
                #if stats['n_psms'] < 10:
                #    print("\t".join(row))
                #if stats['n_psms'] >1000:
                    #auxiliary.print_tree(psm)
                    #sys.exit(10)

                #### Update counters and print progress
                stats['n_psms'] += 1
                if self.verbose >= 1:
                    if stats['n_psms']/5000 == int(stats['n_psms']/5000):
                        if not progress_intro:
                            eprint("INFO: Reading psms.. ", end='')
                            progress_intro = True
                        eprint(f"{stats['n_psms']}.. ", end='', flush=True)

        infile.close()
        if self.verbose >= 1:
            eprint('')

        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"INFO: Read {stats['n_psms']} psms from {self.pepxml_file}")

        if self.verbose >= 1:
            print(f"INFO: Elapsed time: {t1-t0}")
            print(f"INFO: Processed {stats['n_psms']/(t1-t0)} psms per second")
            print(f"INFO: Found data from {len(self.mzmls)} mzML files")


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Reads a PepXML file')
    argparser.add_argument('--min_psm_probability', action='store', default=0.90, help='Minimum PSM (PeptideProphet or iProphet) probability to accept')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('pepxml_file', type=str, help='Filename of the PeptideProphet or iProphet PepXML file from first pass')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 0

    #### Ensure that the PepXML file is really there before starting work
    input_file = params.pepxml_file
    if not os.path.isfile(input_file):
        print(f"ERROR: File '{input_file}' not found or not a file")
        return

    #### Load the relevant search result information from the pepXML file
    search_results = SimplePepXmlReader()
    search_results.read(input_file, verbose=verbose)


#### For command line usage
if __name__ == "__main__": main()
