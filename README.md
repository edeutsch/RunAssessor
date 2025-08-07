# RunAssessor
Assesses an MS Run mzML file for important information de novo


## Installation
Clone the repo and install dependencies:
```bash
git clone https://github.com/edeutsch/RunAssessor.git
cd RunAssesor
pip install -r requirements.txt
```

## Quick Start

Code processes mass spectrometry data in the `.mzML` format. `.mzML.gz` files are also supported.
```bash
python bin/assess_mzMLs.py ~/file_path/mass_spec_data.mzML.gz
```

Code can also accept multiple `.mzML` files at a time 
```bash
python bin/assess_mzMLs.py ~/file_path/*.mzML # can add .gz to end
```

# Usage

## assess_mzMLs.py

### Command-line Arguments

| Argument | Description | Default | Example |
|----------|-------------|---------|---------|
| `--metadata_filepath` | Filepath of the metadata file. | `study_metadata.json` | `--metadata_filepath ~/file_path/my_metadata.json` |
| `--preserve_existing_metadata` | If set, read the existing metadata file, preserve its contents, and build on it (default is to overwrite current file). | Not set | `--preserve_existing_metadata` |
| `--n_threads` | Number of files to read in parallel. | Number of CPU cores | `--n_threads 4` |
| `--write_fragmentation_type_files` | If set, writes a fragmentation type file for each mzML input. | Not set | `--write_fragmentation_type_files` |
| `--verbose` | If set, prints more detailed processing information. Use multiple times to increase verbosity. | Not set | `--verbose` |
| `--version` | Prints the version of the script/tool. | 0.8 | `--version` |
| `--write_sdrf_file` | If set, writes an SDRF file based on the metadata text file. | Not set | `--write_sdrf_file` |
| `--write_pdfs` | If set, generates a PDF of delta time and ppm graphs from precursor stats, and a PDF of neutral loss windows from composite intensities. | Not set | `--write_pdfs` |


## assess_metadata.py

### Comand-line Arguments
| Argument | Description | Default | Example |
|----------|-------------|---------|---------|
| `--metadata_file` | Specify a metadata file if different than the default. | `study_metadata.json` | `--metadata_file ~/file_path/custom_metadata.json` |
| `--write_sdrf_file` | If set, writes an SDRF file based on the metadata file. Can be passed multiple times. | Not set | `--write_sdrf_file` |
| `--verbose` | If set, enables verbose output. Use multiple times to increase verbosity. | Not set | `--verbose` |


## rename_msruns.py

## Interpretting JSON file: `study_metadata.json`
