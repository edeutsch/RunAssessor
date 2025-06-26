#!/usr/bin/env python3

import sys
import os
import argparse
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../lib")
from metadata_handler import MetadataHandler


####################################################################################################
#### Main function for command-line usage
def main():

    #### Parse command line arguments
    argparser = argparse.ArgumentParser(description='Study metadata processor')
    argparser.add_argument('--metadata_file', action='store', help='Specify a metadata file if different than the default study_metadata.json')
    argparser.add_argument('--write_sdrf_file', action='count', help='If set, then write an SDRF file based on the metadata file')
    argparser.add_argument('--verbose', action='count' )
    params = argparser.parse_args()

    #### Create metadata handler object and read or create the metadata structure
    metadata = MetadataHandler(params.metadata_file, params.verbose)
    metadata.read_or_create()

    #### Try to infer search criteria based on the information we have available
    metadata.infer_search_criteria()

    #### Store the metadata structure
    metadata.store()

    #### If selected, also write an SDRF file
    if params.write_sdrf_file:
        key_value_file = metadata.find_txt_file()
        metadata.read_txt_file(key_value_file)
        metadata.generate_sdrf_table()  
        sdrf_filename = metadata.infer_sdrf_filename()
        metadata.write_sdrf_file(sdrf_filename) 



if __name__ == "__main__": main()
