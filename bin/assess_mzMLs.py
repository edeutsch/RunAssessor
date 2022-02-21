#!/usr/bin/env python3

#### Define eprint() as printing to stderr
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import multiprocessing
import copy

sys.path.append(os.path.dirname(os.path.realpath(__file__))+"/../lib")
from mzML_assessor import MzMLAssessor
from metadata_handler import MetadataHandler


####################################################################################################
#### Process one mzML
def process_job(job):

    #### Assess the mzML file
    assessor = MzMLAssessor(job['filename'], metadata=job['metadata'], verbose=job['verbose'])
    assessor.read_header()
    assessor.read_spectra(write_fragmentation_type_file=job['write_fragmentation_type_files'])
    assessor.assess_composite_spectra()
    if assessor.metadata['state']['status'] != 'ERROR':
        assessor.assess_ROIs()
    return assessor


####################################################################################################
#### Main function for command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Creates an index for an MSP spectral library file')
    argparser.add_argument('--use_cache', action='store', help='Set to true to use the cached file instead of regenerating')
    argparser.add_argument('--write_fragmentation_type_files', action='count', help='If set, write a fragmentation_type file for each mzML')
    argparser.add_argument('--n_threads', action='store', type=int, help='Set the number of files to process in parallel (defaults to number of cores)')
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

    # Set up the list of jobs to process
    jobs = []

    #### Loop over all the files and put them in a list of jobs for parallel processing
    for file in params.files:

        #### If this file is already in the metadata store, purge it and generate new
        if file in study.metadata['files']:
            if params.refresh is not None and params.refresh > 0:
                study.metadata['files'][file] = None
            else:
                if verbose >= 1: eprint(f"INFO: Already have results for '{file}'. Skipping..")
                continue

        job = { 'filename': file, 'metadata': copy.deepcopy(study.metadata), 'verbose': verbose,
            'write_fragmentation_type_files': params.write_fragmentation_type_files }
        jobs.append(job)

    #### Now process the jobs in parallel
    n_threads = params.n_threads or multiprocessing.cpu_count()
    eprint(f"Processing files with n_threads={n_threads} (one mzML per thread)", end='', flush=True)
    pool = multiprocessing.Pool(processes=n_threads)
    results = pool.map_async(process_job, jobs)
    #results = pool.map(process_job, jobs)
    pool.close()
    pool.join()
    eprint("")

    # Coalesce the results
    results = results.get()
    for result in results:
        for filename in result.metadata['files']:
            study.metadata['files'][filename] = result.metadata['files'][filename]
        study.metadata['state'] = result.metadata['state']

    #### Infer parameters based on the latest data
    study.infer_search_criteria()

    #### Write out our state of mind
    study.store()

if __name__ == "__main__": main()
