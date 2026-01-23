# RunAssessor
The RunAssessor package provides functionality for assessing a set of input mass spectrometry proteomics mzML files, extracting all information that can
be gleaned from the files prior to any attempt to identify the peptides therein.


## Installation
Clone the repo and install dependencies:
```bash
git clone https://github.com/edeutsch/RunAssessor.git
cd RunAssesor
pip install -r requirements.txt
```

## Quick Start

RunAssessor processes mass spectrometry data in the `.mzML` format. `.mzML.gz` files are also supported.
```bash
python ~/RunAssessor/bin/assess_mzMLs.py ~/file_path/mass_spec_data.mzML.gz
```

RunAssessor can also accept multiple `.mzML` or '.mzML.gz' files at a time 
```bash
python ~/RunAssessor/bin/assess_mzMLs.py --write_sdrf_files --verbose *.mzML
```

### Running a quick test

In `RunAssessor/tests/data/` there is a small `.mzML.gz` file that can be used to quickly try RunAssessor
```bash
cd tests/data
python ../../bin/assess_mzMLs.py chlo_6_tiny.mzML.gz
```

There is also a set of larger `.mzML.gz` files of different kinds that are used to test RunAssessor in `RunAssessor/tests/large_files/`.
You can try RunAssessor on those files. Files in `large_files` will generate a more substantial JSON and summary file compared to `chlo_6_tiny.mzML.gz`.
```bash
cd tests/large_files
python large_test_file_downloader.py
python ../../bin/assess_mzMLs.py *.mzML.gz   # Or single file of your choice
```


# Usage of the main RunAssessor programs

## assess_mzMLs.py

The assess_mzMLs.py program is the main RunAssessor program that accepts a set of input mzML files along with various parameters and processes
the mzML files accordingly. This is generally the only program that needs to be executed.

### Command-line Arguments

| Argument | Description | Default | Example |
|----------|-------------|---------|---------|
| `--metadata_filepath` | Filepath of the metadata file. | `study_metadata.json` | `--metadata_filepath ~/file_path/my_metadata.json` |
| `--preserve_existing_metadata` | If set, read the existing metadata file, preserve its contents, and build on it (default is to overwrite current file). | Not set | `--preserve_existing_metadata` |
| `--n_threads` | Number of files to read in parallel. | Number of CPU cores | `--n_threads 4` |
| `--write_fragmentation_type_files` | If set, writes a fragmentation type file for each mzML input. | Not set | `--write_fragmentation_type_files` |
| `--verbose` | If set, prints more detailed processing information. Use multiple times to increase verbosity. | Not set | `--verbose` |
| `--version` | Prints the version of the script/tool. | Not set | `--version` |
| `--write_sdrf_file` | If set, writes an SDRF file based on the metadata text file. | Not set | `--write_sdrf_file` |
| `--include_sdrf_provenance` | If set, then all information written to the SDRF file will have provenance tags (this is experimental non-standard SDRF, not for routine use). | Not set | `--include_sdrf_provenance` |
| `--write_pdfs` | If set, generates a PDF of delta time and ppm graphs from precursor stats, and a PDF of neutral loss windows from composite intensities. | Not set | `--write_pdfs` |
| `--write_ions` | If set, writes a TSV file of the ion three-sigma table, which lists all fragment ions found. | Not set | `--write_ions` |


## assess_metadata.py

The assess_metadata.py program reads an existing study_metadata.json file and recomputes recommended parameters and potentially outputs to
different formats. This program is generally not routinely used, but can be helpful if, for example, a large set of input files are analyzed
and a study_metadata.json file is written, but the --write_sdrf_file was not provided to assess_mzMLs.py. This tool can take the previous analysis
and write out the SDRF file based on the existing JSON file, without reprocessing all mzML files again.

### Command-line Arguments
| Argument | Description | Default | Example |
|----------|-------------|---------|---------|
| `--metadata_file` | Specify a metadata file if different than the default. | `study_metadata.json` | `--metadata_file ~/file_path/custom_metadata.json` |
| `--write_sdrf_file` | If set, writes an SDRF file based on the metadata file. Can be passed multiple times. | Not set | `--write_sdrf_file` |
| `--verbose` | If set, enables verbose output. Use multiple times to increase verbosity. | Not set | `--verbose` |


## rename_msruns.py

Renames raw files in current directory to replace special characters in filenames with underscores and writes a mapping document of changes.
Some special characters in filename can cause problems in downstream processing tools (for example: spaces, parentheses, ampersands, and more). This program
fixes the filenames in a repeatable way and records how the filenames have been adjusted in an output file. This is not needed for normal RunAssessor operation,
but may be useful for large-scale reprocessing efforts when certain tools have trouble with special characters in filenames.




