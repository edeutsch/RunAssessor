import os
import sys
import requests
import time

import argparse
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


def download_large_test_files(urls, output_filenames=None):
    dir = os.path.dirname(os.path.abspath(__file__))
    
    # Download large test file from each url in the list of urls
    for i, url in enumerate(urls):
        output_filename = None
        if output_filenames:
            output_filename = output_filenames[i]
        if output_filename is None:
            file_name = url.split('/')[-1]

            if '?' in file_name:
                file_name = file_name.split('?')[0]
            if not file_name:
                eprint(f"ERROR: Could not determine file name from URL: '{url}'. Skipping this download.")
                continue
        else:
            file_name = output_filename

        file_path = os.path.join(dir, file_name)
    
        if not os.path.exists(file_path):
            eprint(f"INFO: Downloading large test file '{file_name}' from URL: '{url}'")
            eprint(f"INFO: Saving '{file_name}' to file path: '{file_path}'")

            response = requests.get(url, stream=True)
            if response.status_code == 200:
                with open(file_path, "wb") as f:
                    # Tracking progress on file download
                    total_size = int(response.headers.get('content-length', 0))
                    block_size = 8192
                    downloaded_size = 0
                    start_time = time.time()

                    # Downloading file chunk by chunk
                    for chunk in response.iter_content(chunk_size=block_size):
                        if chunk:
                            f.write(chunk)
                            downloaded_size += len(chunk)
                            if total_size > 0:
                                progress = (downloaded_size / total_size) * 100
                                elapsed_time = time.time() - start_time
                                speed = (downloaded_size / (1024 * 1024)) / elapsed_time if elapsed_time > 0 else 0
                                sys.stderr.write(f"\rINFO: Downloading: {progress:.2f}% ({downloaded_size / (1024*1024):.2f} MB / {total_size / (1024*1024):.2f} MB) | Speed: {speed:.2f} MB/s")
                                sys.stderr.flush()
                    sys.stderr.write("\n")

                eprint(f"INFO: File downloaded successfully and saved as '{file_name}' to '{file_path}' in {elapsed_time:.2f} seconds")
            else:
                eprint(f"ERROR: Failed to download file. HTTP status code: {response.status_code} for URL: '{url}'")
        else:
            eprint(f"INFO: Large test file '{file_name}' already found. Skipping download.")
        
        sys.stderr.write("\n")

    # Check to make sure the number of URLs is the same as the number of output_filenames
    if output_filenames is not None and len(urls) != len(output_filenames):
        eprint("ERROR: Number of urls does not match number of output_filenames")
        raise ValueError("Number of urls does not match number of output_filenames")



if __name__ == "__main__":
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
# To run the large_test_file_downloader.py in command line,
# cd tests\large_files
# python large_test_file_downloader.py