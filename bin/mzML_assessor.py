#!/usr/bin/env python3

#### Define eprint() as printing to stderr
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../lib")
from mzML_assessor import MzMLAssessor
from metadata_handler import MetadataHandler


####################################################################################################
#### Main function for command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Creates an index for an MSP spectral library file')
    argparser.add_argument('--use_cache', action='store', help='Set to true to use the cached file instead of regenerating')
    argparser.add_argument('--refresh', action='count', default=0, help='If set, existing metadata for a file will be overwritten rather than skipping the file')
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

    #### Load the current metadata file in order to update it
    study = MetadataHandler(verbose=verbose)
    study.read_or_create()

    #### Loop over all the files processing them
    for file in params.files:

        #### If this file is already in the metadata store, purge it and generate new
        if file in study.metadata['files']:
            if params.refresh is not None and params.refresh > 0:
                study.metadata['files'][file] = None
            else:
                if verbose >= 1: eprint(f"INFO: Already have results for '{file}'. Skipping..")
                continue

        #### Assess the mzML file
        assessor = MzMLAssessor(file, metadata=study.metadata, verbose=verbose)
        assessor.read_header()
        if not params.use_cache:
            assessor.read_spectra()
        assessor.assess_composite_spectra()
        assessor.assess_ROIs()

    #### Infer parameters based on the latest data
    study.infer_search_criteria()

    #### Write out our state of mind
    study.store()

if __name__ == "__main__": main()
