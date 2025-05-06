#!/usr/bin/env python3

import sys
import os
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../lib")
from metadata_handler import MetadataHandler


####################################################################################################
#### Main function for command-line usage
def main():

    #### Parse command line arguments
    argparser = argparse.ArgumentParser(description='handler for study metadata files')
    argparser.add_argument('--metadata_file', action='store', help='Specify a metadata file if different than the default')
    argparser.add_argument('--verbose', action='count' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    params = argparser.parse_args()

    #### Create metadata handler object and read or create the metadata structure
    metadata = MetadataHandler(params.metadata_file, params.verbose)
    metadata.read_or_create()

    #### Try to infer search criteria based on the information we have available
    metadata.infer_search_criteria()

    #### Store the metadata structure
    metadata.store()

if __name__ == "__main__": main()
