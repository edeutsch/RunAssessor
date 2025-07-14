import os
import sys
import shutil
import argparse

def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


def update_expected_output(files):
    large_files_dir = os.path.dirname(os.path.abspath(__file__))
    tests_dir = os.path.dirname(large_files_dir)
    expected_output_dir = os.path.join(tests_dir, "expected_output")

    total_files = len(files)

    for i, file in enumerate(files):

        current_path = os.path.join(large_files_dir, file)
        destination_path = os.path.join(expected_output_dir, file)

        continuation_message = ""
        if i < total_files - 1:
            continuation_message = "... continuing to next file"

        if not os.path.exists(current_path):
            eprint(f"ERROR: '{file}' not found in 'tests/large_files' directory{continuation_message}")
            continue
        
        try:
            if os.path.exists(destination_path):
                shutil.copy2(current_path, expected_output_dir)
                eprint(f"INFO: Overwriting original expected '{file}' in 'tests/expected_output' directory...")
                eprint(f"INFO: New '{file}' successfully copied into 'tests/expected_output' directory{continuation_message}")
            else:
                shutil.copy2(current_path, expected_output_dir)
                eprint(f"INFO: New '{file}' successfully copied into 'tests/expected_output' directory{continuation_message}")
        except (FileExistsError, FileNotFoundError, OSError, PermissionError, shutil.Error) as e:
            eprint(f"ERROR: {e}: Unable to copy '{file}' into 'tests/expected_output' directory{continuation_message}")
            continue
    


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description="Copy over one or more study_metadata.json files from 'tests/large_files' to 'tests/expected_output'"
    )

    argparser.add_argument(
        'files',
        type=str,
        nargs='+',
        help="List of file names of one or more study_metadata.json files to copy from 'tests/large_files' to tests/expected_output'"
    )

    params = argparser.parse_args()

    update_expected_output(files=params.files)
# To run the update_expected_outputs.py in command line:
# cd tests\large_files
# python update_expected_outputs.py <json file 1> <json file 2> ...