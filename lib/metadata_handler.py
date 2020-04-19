#!/usr/bin/env python3

#### Define eprint() as printing to stderr
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import timeit
import re
import json



####################################################################################################
#### Metadata handler class
class MetadataHandler:


    ####################################################################################################
    #### Constructor
    def __init__(self, metadata_file=None, verbose=None):

        #### Check if a metadata_file was provided and if not, set to default
        if metadata_file is None or metadata_file == '':
            metadata_file = 'study_metadata.json'
        self.metadata_file = metadata_file

        #### Check if a verbose was provided and if not, set to default
        if verbose is None:
            verbose = 1
        self.verbose = verbose


    ####################################################################################################
    #### Read or create a study metadata file
    def read_or_create(self):

        file = self.metadata_file

        #### If this is a directory, then just append the default filename
        if os.path.isdir(file):
            file += '/study_metadata.json'

        #### If the specified (or inferred) file does not exist, we should create it
        if not os.path.isfile(file):
            if self.verbose >= 1: eprint(f"INFO: File '{file}' not found or not a file. Will create it.")
            self.create_template()
            with open(file,'w') as outfile:
                json.dump(self.metadata,outfile,sort_keys=True,indent=2)
            return

        #### If there is such a file, read it
        with open(file,'r') as infile:
            self.metadata = json.load(infile)

        #### Verify that this is the right kind of JSON
        if 'knowledge' not in self.metadata or 'state' not in self.metadata:
            eprint(f"ERROR: File '{file}' is JSON, but does not appear to be a study metadata file")

        #### Verify the version number of the file
        if 'version' not in self.metadata['state']:
            eprint(f"ERROR: File '{file}' does not have version information. Something is wrong")
        if self.metadata['state']['version'] < 'v0.1':
            eprint(f"ERROR: Version {self.metadata['state']['version']} in file '{file}' is too old and is not supported anymore")

        #### If verbose, print some information
        if self.verbose >= 1:
            eprint(f"INFO: Study metadata file '{file}' is loaded")
            eprint(f"INFO: Status is {self.metadata['state']['status']}")
            if self.metadata['state']['status'] != 'OK':
                eprint(f"INFO: Code is '{self.metadata['state']['code']}'")
                eprint(f"INFO: Message is '{self.metadata['state']['message']}'")

        #### If we got this far, then everything seems okay
        return


    ####################################################################################################
    #### Create a blank template metadata object
    def create_template(self):

        self.metadata = {
            'files': {},
            'knowledge': {},
            'search_criteria': {},
            'problems': {
                'warnings': {
                    'count': 0,
                    'list': [],
                    'codes': {} },
                'errors': {
                    'count': 0,
                    'list': [],
                    'codes': {} }
            },
            'state': { 'status': 'OK', 'message': 'No problems detected', 'code': 'OK', 'version': 'v0.1' }
        }


    ####################################################################################################
    #### Store study metadata file
    def store(self):

        if self.verbose >= 1: eprint(f"INFO: Writing study metadata file '{self.metadata_file}'")
        with open(self.metadata_file,'w') as outfile:
            json.dump(self.metadata,outfile,sort_keys=True,indent=2)


    ####################################################################################################
    #### Infer search criteria based on available information
    def infer_search_criteria(self):

        #### Short handle for the search criteria information
        criteria = self.metadata['search_criteria']
        knowledge = self.metadata['knowledge']
        verbose = self.verbose
        if verbose >= 1: eprint("INFO: Inferring search criteria from the available information")

        # Since we will recompute the numbers, zero out all the counters
        knowledge['instrument_model'] = 'unknown'
        if 'instrument_models' in knowledge:
            for instrument_model in knowledge['instrument_models']:
                knowledge['instrument_models'][instrument_model] = 0

        criteria['fragmentation_type'] = 'unknown'
        if 'fragmentation_types' in criteria:
            for fragmentation_type in criteria['fragmentation_types']:
                criteria['fragmentation_types'][fragmentation_type] = 0

        #### Loop over all the files and decide on search criteria
        for file in self.metadata['files']:
            fileinfo = self.metadata['files'][file]

            #### Assemble consensus instrument
            if 'instrument_model' in fileinfo:
                file_instrument = fileinfo['instrument_model']['name'] or 'unknown'
                if 'instrument_model' not in knowledge: knowledge['instrument_model'] = file_instrument
                if 'instrument_models' not in knowledge: knowledge['instrument_models'] = {}
                if file_instrument not in knowledge['instrument_models']: knowledge['instrument_models'][file_instrument] = 0
                knowledge['instrument_models'][file_instrument] += 1
                if knowledge['instrument_model'] == 'unknown':
                    knowledge['instrument_model'] = file_instrument
                elif knowledge['instrument_model'] == file_instrument:
                    pass
                else:
                    knowledge['instrument_model'] = 'multiple'
                    self.log_event('ERROR','MultipleInstruments',f"There are multiple instruments in this group of MS runs. Unless they are highly similar, it may be best to separate them.")
            else:
                knowledge['instrument_model'] = 'unknown'
                self.log_event('ERROR','MissingInstrument',f"The instrument was not determined. This should be handled better.")

            #### Assemble fragmentation type
            if 'fragmentation_type' in fileinfo['spectra_stats']:
                fragmentation_type = fileinfo['spectra_stats']['fragmentation_type'] or 'unknown'
                if 'fragmentation_type' not in criteria: criteria['fragmentation_type'] = fragmentation_type
                if 'fragmentation_types' not in criteria: criteria['fragmentation_types'] = {}
                if fragmentation_type not in criteria['fragmentation_types']: criteria['fragmentation_types'][fragmentation_type] = 0
                criteria['fragmentation_types'][fragmentation_type] += 1
                if criteria['fragmentation_type'] == 'unknown':
                    criteria['fragmentation_type'] = fragmentation_type
                elif criteria['fragmentation_type'] == fragmentation_type:
                    pass
                else:
                    criteria['fragmentation_type'] = 'multiple'
                    self.log_event('ERROR','MultipleFragTypes',f"There are multiple fragmentation types in this group of MS runs. Unless your search engine can handle this, you should split them.")
            else:
                criteria['fragmentation_type'] = 'unknown'
                self.log_event('ERROR','MissingFragType',f"The instrufragmentation type was not determined. This should be handled better.")

            #### Assemble high_accuracy_precursors
            if 'high_accuracy_precursors' not in criteria: criteria['high_accuracy_precursors'] = 'unknown'
            high_accuracy_precursors = fileinfo['spectra_stats']['high_accuracy_precursors']
            if criteria['high_accuracy_precursors'] == 'unknown':
                criteria['high_accuracy_precursors'] = high_accuracy_precursors
            elif criteria['high_accuracy_precursors'] == high_accuracy_precursors:
                pass
            else:
                criteria['high_accuracy_precursors'] = 'multiple'
                self.log_event('ERROR','MultipleFragTypes',f"There are multiple precursor accuracy types in this group of MS runs. They should probably be separated.")

            #### Assemble labeling
            if 'labeling' not in criteria: criteria['labeling'] = 'unknown'
            labeling = fileinfo['summary']['labeling']['call']
            if criteria['labeling'] == 'unknown':
                if labeling != 'ambiguous':
                    criteria['labeling'] = labeling
            elif criteria['labeling'] == labeling:
                pass
            else:
                if labeling != 'ambiguous':
                    criteria['labeling'] = 'multiple'
                    self.log_event('ERROR','MultipleLabelingTypes',f"There are multiple labeling types in this MS run group. Split them.")


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
