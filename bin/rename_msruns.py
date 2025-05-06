#!/usr/bin/env python3

import sys
import os
import argparse
import re
import traceback
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


####################################################################################################
#### get a list of files to process
def fix_msrun_names(verbose=0):

    #### Get a list of current files in the cwd
    files = os.listdir('.')
    rawfile_list = []
    for filename in files:
        if filename.endswith('.raw') or filename.endswith('.RAW'):
            rawfile_list.append(filename)

    #### If there are no .raw or .RAW files, then we don't know how to proceed
    #### FIXME. Extend to .wiff and .wiff.scan and other kinds of files
    if len(rawfile_list) == 0:
        print(f"INFO: Did not find any RAW files in the current working directory. Nothing to do")
        return

    #### Check to see if there is already a map file here
    map_file = 'renamed_msruns.tsv'
    map_file_mode = 'w'
    if os.path.exists(map_file):
        print(f"WARNING: There is already a '{map_file}' in this directory. rename_msruns may already have been run?")
        map_file_mode = 'a'

    is_map_file_open = False
    n_renamed_files = 0

    #### Check the filenames and fix if necessary
    for filename in rawfile_list:

        #### Find and strip the filename extension
        match = re.search(r'\.raw$',filename,flags=re.IGNORECASE)
        if match:
            suffix = match.group(0)
            filename = re.sub(r'\.raw$','',filename,flags=re.IGNORECASE)
        else:
            print("ERROR: File does not end with .RAW")
            return

        #### Substitute all non A-Za-z0-9_- characters with a _
        new_filename = re.sub(r'[^\w-]','_',filename)
        if new_filename != filename:

            #### If the map file is not open yet, open it
            if not is_map_file_open:
                try:
                    outfile = open(map_file,map_file_mode)
                except Exception as error:
                    exception_type, exception_value, exception_traceback = sys.exc_info()
                    print(f"ERROR: Unable to open {map_file}: {error}: {repr(traceback.format_exception(exception_type, exception_value, exception_traceback))}")
                    return
                is_map_file_open = True

            #### Rename the file to the new one
            print(f"INFO: Rename {filename}{suffix} to {new_filename}{suffix}")
            try:
                os.rename(f"{filename}{suffix}", f"{new_filename}{suffix}")
            except Exception as error:
                exception_type, exception_value, exception_traceback = sys.exc_info()
                print(f"ERROR: Unable to rename file {filename}{suffix}: {error}: {repr(traceback.format_exception(exception_type, exception_value, exception_traceback))}")
                outfile.close()
                return

            #### Write out the mapping in to the map_file
            outfile.write(f"{filename}\t{new_filename}\n")
            n_renamed_files += 1


    #### Inform of the final result and close the map_file
    if n_renamed_files == 0:
        print(f"INFO: No files needed renaming. Nothing done.")
    else:
        print(f"INFO: Renamed {n_renamed_files} files, and map file {map_file} written")
        outfile.close()


####################################################################################################
#### Main function for command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Renames raw files in the current directory to replace weird characters with underscore and write a lookup table of changes')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None: verbose = 0

    #### Get the current list of msruns
    fix_msrun_names(verbose=verbose)


if __name__ == "__main__": main()
