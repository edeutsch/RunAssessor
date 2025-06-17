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
python bin/assess_mzMLs.py ~/file_path/*.mzML
```

# Usage

## assess_metadata.py

### Command-line Arguments

| Argument | Description | Default | Example |
|----------|-------------|---------|---------|
| `--metadata_filepath` | Filepath of the metadata file. | `study_metadata.json` | `--metadata_filepath ./data/my_metadata.json` |
| `--preserve_existing_metadata` | If set, read the existing metadata file, preserve its contents, and build on it (default is to overwrite current file). | Not set | `--preserve_existing_metadata` |
| `--n_threads` | Number of files to read in parallel. | Number of CPU cores | `--n_threads 4` |
| `--write_fragmentation_type_files` | If set, writes a fragmentation type file for each mzML input. | Not set | `--write_fragmentation_type_files` |
| `--verbose` | If set, prints more detailed processing information. Use multiple times to increase verbosity. | Not set | `--verbose` |
| `--version` | Prints the version of the script/tool. | 0.8 | `--version` |


## assess_mzMLs.py

## rename_msruns.py

## Interpretting JSON file: `study_metadata.json`
