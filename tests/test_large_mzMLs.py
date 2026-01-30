import sys
import os
import subprocess
import json
from deepdiff import DeepDiff
import pytest
import glob

def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


# Set up to run assess_mzMLs.py on the large test files '*.mzML.gz'
tests_dir = os.path.dirname(os.path.abspath(__file__))
print(tests_dir)
project_root = os.path.dirname(tests_dir)
print(project_root)

bin_dir = os.path.join(project_root, "bin")
assess_mzMLs_file = os.path.join(bin_dir, "assess_mzMLs.py")

tests_large_files_dir = os.path.join(tests_dir, "large_files")
input_mzML_file = [     # Located in tests/large_files directory which is the current working directory for subprocess
    "P50.mzML.gz",
    "QEX03_210305_PTMscan_wt_IP_63.mzML.gz",
    "06CPTAC_BCprospective_W_BI_20161116_BA_f17.mzML.gz",
    "Q20181210_06.mzML.gz",
    "140924_11.mzML.gz",
    "ALS081215_10.mzML.gz",
    "Orbi_CID_HCD_iTRAQ8plex.mzML.gz",
    "QExactive_TMT6.mzML.gz",
    "ADPhosTMT3R2.mzML.gz",
    "20250613_PRM34_AST_EV1109_G20_PRM_Pool10_iRT25_HeLaS3-P1.mzML.gz",
    "20171001_QE3_nLC7_DBJ_SA_LFQphos_Tech_Rep_03.mzML.gz",
    "ea02315_OrbiterFAIMS_CHX_HEK293T_F6.mzML.gz",
    "CRC_phospho_iTRAQ_20.mzML.gz",
    "20180914_QE8_nLC0_BDA_SA_DIA_Skin_Dermal_T_cells_CD4_MT_200000_2.mzML.gz",
    "EXP22038_2022ms0310aX49_A_BA7_1_9837.mzML.gz",
    "WS1_11.mzML.gz",
    "iTRAQ_Phos_T1.mzML.gz",
    "C_i012_PO4_SCX_15_ETD_qMODE.mzML.gz",
    "20150730_LWP_HBE_1.mzML.gz"
]

# The current output of 'assess_mzMLs.py' is '<root_name>_study_metadata.json' in the tests/large_files directory
# The expected output is '<root_name>_study_metadata.json' in the tests/expected_output directory
tests_expected_output_dir = os.path.join(tests_dir, "expected_output")

if tests_large_files_dir not in sys.path:
    sys.path.insert(0, tests_large_files_dir)
from large_test_file_downloader import download_large_test_files


@pytest.mark.slow
def test_large_files_present():
    download_large_test_files(urls=[
        "https://peptideatlas.org/data/RunAssessor/testfiles/P50.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/QEX03_210305_PTMscan_wt_IP_63.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/06CPTAC_BCprospective_W_BI_20161116_BA_f17.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/Q20181210_06.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/140924_11.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/ALS081215_10.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/Orbi_CID_HCD_iTRAQ8plex.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/QExactive_TMT6.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/ADPhosTMT3R2.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/20250613_PRM34_AST_EV1109_G20_PRM_Pool10_iRT25_HeLaS3-P1.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/20171001_QE3_nLC7_DBJ_SA_LFQphos_Tech_Rep_03.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/ea02315_OrbiterFAIMS_CHX_HEK293T_F6.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/CRC_phospho_iTRAQ_20.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/20180914_QE8_nLC0_BDA_SA_DIA_Skin_Dermal_T_cells_CD4_MT_200000_2.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/EXP22038_2022ms0310aX49_A_BA7_1_9837.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/WS1_11.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/iTRAQ_Phos_T1.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/C_i012_PO4_SCX_15_ETD_qMODE.mzML.gz",
        "https://peptideatlas.org/data/RunAssessor/testfiles/20150730_LWP_HBE_1.mzML.gz"
    ]
    )
    for input_file in input_mzML_file:
        assert os.path.exists(os.path.join(tests_large_files_dir, input_file))


# Mark following test as slow and skip it unless specify "--runslow" in command line
@pytest.mark.slow
def test_compare_study_metadatas():     # This test takes ~ 1 hour

    keep_output = False
    halt_if_different = True

    study_metadata_files = [
        os.path.join(tests_large_files_dir, "*_study_metadata.json"),
        os.path.join(tests_large_files_dir, "*_study_metadata.file_summary.tsv"),
        os.path.join(tests_large_files_dir, "*_study_metadata.txt")
]
    for file_type in study_metadata_files:
        for file_path in glob.glob(file_type):
            if os.path.exists(file_path):
                os.remove(file_path)

    for i, input_file in enumerate(input_mzML_file):
        
        file_root_name = input_file.split('.')[0]
        current_study_metadata = os.path.join(tests_large_files_dir, f"{file_root_name}_study_metadata.json")
        metadata_filepath = f"{file_root_name}_study_metadata.json"

        output = subprocess.run(
            [sys.executable, assess_mzMLs_file, input_file, '--verbose', '--metadata_filepath', metadata_filepath],
            cwd=tests_large_files_dir,
            capture_output=True,
            text=True,
            check=True
        )

        # Testing verbosity
        stdout_lines = []
        for line in output.stdout.splitlines():
            if line != '':
                stdout_lines.append(line)
        print(stdout_lines)

        stderr_lines = []
        for line in output.stderr.splitlines():
            if line != '':
                stderr_lines.append(line)
        eprint(stderr_lines)

        assert len(stdout_lines) == 0
        assert len(stderr_lines) >= 11

        # Make sure the current generated '<root_name>_study_metadata.json' files exactly match the expected output files
        file_root = input_file.replace('.mzML.gz','')
        file_root = file_root.replace('.mzML','')
        expected_study_metadata = file_root + '_study_metadata.json'
        expected_study_metadata = os.path.join(tests_expected_output_dir, expected_study_metadata)
        with open(current_study_metadata, 'r') as study_metadata, open(expected_study_metadata, 'r') as expected_study_metadata:
            current = json.load(study_metadata)
            expected = json.load(expected_study_metadata)
        
        diff = DeepDiff(expected, current)
        if halt_if_different:
            assert not diff, \
                f"Current generated study_metadata.json for '{input_file}' does not match expected output." \
                f"\n\nDifferences found:\n{diff.pretty()}\n"
        else:
            print(f"Current generated study_metadata.json for '{input_file}' does not match expected output." \
                f"\n\nDifferences found:\n{diff.pretty()}\n")

        if not keep_output:
            study_metadata_files_glob = os.path.join(tests_large_files_dir, f"{file_root}_study_metadata.*")
            for file_path in glob.glob(study_metadata_files_glob):
                if os.path.exists(file_path):
                    os.remove(file_path)
