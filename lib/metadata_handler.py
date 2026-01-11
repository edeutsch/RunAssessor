#!/usr/bin/env python3

import sys
import os
import argparse
import os.path
import json
import shutil
import pandas as pd
import numpy
import datetime
numpy.seterr(invalid='ignore')

import csv
import math
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)


####################################################################################################
#### Metadata handler class
class MetadataHandler:


    ####################################################################################################
    #### Constructor
    def __init__(self, metadata_filepath=None, verbose=None, write_ions=None, frag_units=None):

        default_metadata_filepath = 'study_metadata.json'

        self.metadata = {}
        self.default_filename = False
        #### Check if a metadata_filepath was provided and if not, set to default
        if metadata_filepath is None or metadata_filepath == '':
            metadata_filepath = default_metadata_filepath
            self.default_filename = True

        #### If this is a directory, then just append the default filename
        if os.path.isdir(metadata_filepath):
            metadata_filepath += '/' + default_metadata_filepath

        #### If metadata_filepath does not end with .json, +add a .json 
        if not metadata_filepath.endswith(".json"):
            metadata_filepath += ".json"

        self.metadata_filepath = metadata_filepath
        self.sdrf_hints = {}
        self.sdrf_table_column_titles = []
        self.sdrf_table_rows = []
        self.metadata = {}


        #### Check if a verbose was provided and if not, set to default
        if verbose is None:
            verbose = 0
        self.verbose = verbose

        if write_ions != None:
            self.write_ions = True
        else:
            self.write_ions = False
        
        #### Set up type of units for ms/ms tolerances
        # frag_units can only be mz, ppm, or None
        self.fragmentation_units = frag_units 

        #### Values used to calculate tolerance recomendations
        self.ppm_error = 5
        self.mz_error = 0.3


    ####################################################################################################
    #### Create a study metadata file
    def create(self):

        file = self.metadata_filepath

        #### Comment if the file already exists
        if os.path.isfile(file):
            if self.verbose >= 1:
                eprint(f"INFO: Study metadata file '{file}' already exists. It is being overwritten.")

        if self.verbose >= 1:
            eprint(f"INFO: Creating study metadata file '{file}'")
        self.create_template()
        txt_file = self.find_txt_file()
        self.read_txt_file(txt_file)
        try:
            with open(file, 'w') as outfile:
                json.dump(self.metadata,outfile, sort_keys=True, indent=2)
        except:
            eprint(f"ERROR: Cannot open file '{file}' for writing")
            return

        return 'OK'



    ####################################################################################################
    #### Read or create a study metadata file
    def read_or_create(self):

        file = self.metadata_filepath

        #### If the specified (or inferred) file does not exist, we should create it
        if not os.path.isfile(file):
            if self.verbose >= 1:
                eprint(f"INFO: Study metadata file '{file}' not found or not a file. Create it.")
            result = self.create()
            return result

        #### If there is such a file, read it
        try:
            infile = open(file, 'r')
            try:
                self.metadata = json.load(infile)
            except:
                eprint(f"ERROR: Cannot parse JSON from file '{file}'")
                return
        except:
            eprint(f"ERROR: Cannot open file '{file}' for reading")
            return

        #### Verify that this is the right kind of JSON
        n_errors = 0
        if 'knowledge' not in self.metadata or 'state' not in self.metadata:
            eprint(f"ERROR: File '{file}' is JSON, but does not appear to be a study metadata file")
            n_errors += 1

        #### Verify the version number of the file
        if 'version' not in self.metadata['state']:
            eprint(f"ERROR: File '{file}' does not have version information. Something is wrong")
            n_errors += 1
        if self.metadata['state']['version'] < 'v0.1':
            eprint(f"ERROR: Version {self.metadata['state']['version']} in file '{file}' is too old and is not supported anymore")
            n_errors += 1

        if n_errors > 0:
            return

        #### If verbose, print some information
        if self.verbose >= 1:
            eprint(f"INFO: Study metadata file '{file}' is loaded")
            eprint(f"INFO: Status is {self.metadata['state']['status']}")
            if self.metadata['state']['status'] != 'OK':
                eprint(f"INFO: Code is '{self.metadata['state']['code']}'")
                eprint(f"INFO: Message is '{self.metadata['state']['message']}'")
                return

        #### If we got this far, then everything seems okay
        return 'OK'



    ####################################################################################################
    #### Find and set the key-value study metadata text file
    def find_txt_file(self):

        file = self.metadata_filepath
        #### Replace .json with .txt
        if file.endswith('.json'):
            file = file.replace('.json', '.txt')
        else:
            if self.verbose >= 1:
                eprint(f"INFO: Study metadata file '{file}' does not end in .json, so cannot look for corresponding .txt file")
            return

        #### If the specified key-value file does not exist, no table will be generated. However, if none is specified, a template file will be used
        if not os.path.isfile(file) and not self.default_filename:
            eprint(f"INFO: Looked for but did not find study metadata key-value hints file '{file}'. File not found or not a file. SDRF table will not be generated")
            return None

        elif os.path.isfile(file):
            eprint(f"INFO: {file} found")

        #### Loads Template as key-value pair file
        else:
            if self.verbose >= 1:
                eprint(f"INFO: No study_metadata key-value hints file defined")
                file = os.path.dirname(os.path.abspath(__file__)) + "/study_metadata_template.txt"
                eprint(f"INFO: Using Template key-value hints {file}")

            file = os.path.dirname(os.path.abspath(__file__)) + "/study_metadata_template.txt"
            self.copy_template()

        return file


    ####################################################################################################
    #### Read the key-value study metadata text file
    def read_txt_file(self, txt_file):

        file = txt_file
        if file == None: return
        #### Check to see if file can be read
        try:
            infile = open(file, 'r')
        except:
            eprint(f"ERROR: Cannot open file '{file}' for reading")
            return

        sdrf_hints = { 'key_order': [], 'keys': {} }

        #### Parse the file
        n_errors = 0
        i_line = 0
        for line in infile:
            i_line += 1
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith('#'):
                continue
            position = line.find('=')
            if position < 1:
                eprint(f"ERROR: File '{file}' has malformed key-value pair at line {i_line}: {line}")
                n_errors += 1
            key = line[0:position]
            value = line[position+1:]
            if key not in sdrf_hints['keys']:
                sdrf_hints['keys'][key] = []
                sdrf_hints['key_order'].append(key)
            sdrf_hints['keys'][key].append(value)

        self.sdrf_hints = sdrf_hints

        if n_errors > 0:
            self.metadata['state']['status'] = 'ERROR'
            self.metadata['state']['code'] = 'TxtParseError'
            self.metadata['state']['message'] = f"{n_errors} errors parsing sdrf txt hints file"

        #### If verbose, print some information
        if self.verbose >= 1:
            eprint(f"INFO: Study metadata SDRF hints text file '{file}' is loaded")
            eprint(f"INFO: Status is {self.metadata['state']['status']}")
            if self.metadata['state']['status'] != 'OK':
                eprint(f"INFO: Code is '{self.metadata['state']['code']}'")
                eprint(f"INFO: Message is '{self.metadata['state']['message']}'")
                return

        #### If we got this far, then everything seems okay
        #eprint(json.dumps(sdrf_hints, indent=2, sort_keys=2))
        return 'OK'



    ####################################################################################################
    #### Create a blank template metadata object
    def create_template(self):

        self.metadata = {
            'files': {},
            'knowledge': {},
            'search_criteria': {},
            'spectra_stats': {},
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

        return self.metadata



    ####################################################################################################
    #### Store study metadata file
    def store(self):

        if self.verbose >= 1:
            eprint(f"INFO: Writing study metadata file '{self.metadata_filepath}'")
        with open(self.metadata_filepath, 'w') as outfile:
            json.dump(self.metadata,outfile, sort_keys=True, indent=2)


    ####################################################################################################
    #### Infers the sdrf file name and returns it if possible

    def infer_sdrf_filename(self):
        filename = self.metadata_filepath

        #### Replace .json with .sdrf.tsv
        if filename.endswith('.json'):
            filename = filename.replace('.json', '.sdrf.tsv')
            return filename
        else:
            if self.verbose >= 1:
                eprint(f"INFO: Study metadata file '{filename}' does not end in .json, so cannot create corresponding .sdrf.tsv file")
            return



    ####################################################################################################
    #### Write SDRF file
    def write_sdrf_file(self, filename):

        if len(self.sdrf_hints.get("keys", {})) == 0:
            return
        if self.verbose >= 1:
            eprint(f"INFO: Writing SDRF file '{filename}'")
        with open(filename, 'w') as outfile:
            print("\t".join(self.sdrf_table_column_titles), file=outfile)
            for row in self.sdrf_table_rows:
                print("\t".join(row), file=outfile)

    ####################################################################################################
    #### Creates a study_metadata_template.txt in the directory the user is in
    def copy_template(self):
        eprint("INFO: Creating a study_metadata_template.txt in users directory")
        file = os.path.dirname(os.path.abspath(__file__)) + "/study_metadata_template.txt"
        destination_file = os.getcwd() + "/study_metadata.txt"
        shutil.copy2(file, destination_file)
        if self.verbose >= 1:
            eprint("INFO: study_metadata_template.txt copied over as study_metadata.txt")


    ####################################################################################################
    #### Infer search criteria based on available information
    def infer_search_criteria(self):

        #### Creates a dictionary to store sigma values for all ions in each file
        ion_three_sigma_table = {}
        ### Creates a variable to store 3sigma values across all values
        all_3sigma_values_away = {"Status":False, "Highest": [], "Lowest":[], "precursors":[], "precursor_found": False}

        ### Store information about phosphoenrichment in files
        phosphoenrichment = {'has_phospho':0, 'no_phospho':0}

        ### Creates a list to hold information to generate a table with info about all files
        info = []
        

        #### Short handle for the search criteria information
        criteria = self.metadata['search_criteria']
        knowledge = self.metadata['knowledge']
        spectra_stats = self.metadata['spectra_stats']
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

        #### Loop over all the files twice to extract the start timestamps and compute the relative start times
        earliest_epoch_seconds = 0
        for file in self.metadata['files']:
            fileinfo = self.metadata['files'][file]
            if fileinfo.get('start_timestamp', None) and len(fileinfo['start_timestamp']) > 10:
                try:
                    start_timestamp = datetime.datetime.strptime(fileinfo['start_timestamp'], '%Y-%m-%dT%H:%M:%SZ')
                    fileinfo['delta_start_days'] = start_timestamp.timestamp()
                    if earliest_epoch_seconds == 0 or fileinfo['delta_start_days'] < earliest_epoch_seconds:
                        earliest_epoch_seconds = fileinfo['delta_start_days']
                except:
                    self.log_event('ERROR','Unparsable Timestamp',f"The run timestamp '{fileinfo['start_timestamp']}' cannot be parsed.")
        for file in self.metadata['files']:
            fileinfo = self.metadata['files'][file]
            if 'delta_start_days' in fileinfo:
                fileinfo['delta_start_days'] = round( ( fileinfo['delta_start_days'] - earliest_epoch_seconds ) / (60 * 60 * 24), 5)

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
                self.log_event('ERROR','MissingFragType',f"The fragmentation type was not determined. This should be handled better.")

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
            labeling = 'unknown'

            if 'summary' in fileinfo and 'call' in fileinfo['summary']['combined summary']:
                labeling = fileinfo['summary']['combined summary']['call']
                if criteria['labeling'] == 'unknown':
                    if labeling != 'ambiguous':
                        criteria['labeling'] = labeling
                elif criteria['labeling'] == labeling:
                    pass
                else:
                    if labeling != 'ambiguous':
                        criteria['labeling'] = 'multiple'
                        self.log_event('ERROR','MultipleLabelingTypes',f"There are multiple labeling types in this MS run group. Split them.")


            #### Review the acquisition types
            if 'acquisition_type' not in criteria:
                criteria['acquisition_type'] = 'unknown'
            acquisition_type = fileinfo['spectra_stats']['acquisition_type']
            if criteria['acquisition_type'] == 'unknown':
                criteria['acquisition_type'] = acquisition_type
            elif criteria['acquisition_type'] == acquisition_type:
                pass
            else:
                criteria['acquisition_type'] = 'multiple'
                self.log_event('ERROR','MultipleAcquisitionTypes',f"There are multiple acquisition types in this group of MS runs. They should probably be separated.")


            #### Roll up the spectra counts
            if 'spectra_stats' in fileinfo:
                for key,value in fileinfo['spectra_stats'].items():
                    if key.startswith('n_'):
                        if key not in spectra_stats:
                            spectra_stats[key] = 0
                        spectra_stats[key] += value

            #### Colect sigma values

            try:
                if self.fragmentation_units in {"mz", "ppm"}:
                    units = self.fragmentation_units
                else:
                    frag_type = self.metadata['files'][file]['spectra_stats']['fragmentation_type']
                    if frag_type.startswith('HR'):
                        units = "ppm"
                    elif frag_type.startswith('LR'):
                        units = "mz"
                    else:
                        units = None

                frag_tol = fileinfo['summary']['combined summary']['fragmentation tolerance']

                if units and 'warning' not in frag_tol:
                    if units == "ppm":
                        lower_key = 'fragment_tolerance_ppm_lower'
                        upper_key = 'fragment_tolerance_ppm_upper'
                    else:
                        lower_key = 'lower_m/z'
                        upper_key = 'upper_m/z'

                    all_3sigma_values_away['Lowest'].append(frag_tol[lower_key])
                    all_3sigma_values_away['Highest'].append(frag_tol[upper_key])
                    all_3sigma_values_away['Status'] = True

                #### precursor tolerance (independent of units)
                precursor_tol = fileinfo['summary']['combined summary']['recommended precursor tolerance (ppm)']
                if precursor_tol is not None:
                    all_3sigma_values_away['precursors'].append(precursor_tol)
                    all_3sigma_values_away['precursor_found'] = True
                    
            except:
                pass
           

            #### Add info to ion data
            try:
                for frag_types in fileinfo['lowend_peaks']:
                    for ions in fileinfo['lowend_peaks'][frag_types]:
                        if fileinfo['lowend_peaks'][frag_types][ions]['type'] == "fragment_ion":
                            
                            try:
                                if file not in ion_three_sigma_table:
                                    ion_three_sigma_table[file] = []
                                
                                if self.metadata['files'][file]['spectra_stats']['fragmentation_type'].startswith('HR'):
                                    ion_three_sigma_table[file].append({
                                        "ion": ions,
                                        "three_sigma_lower": fileinfo["lowend_peaks"][frag_types][ions]["peak"]["extended"]["three_sigma_ppm_lower"],
                                        "three_sigma_upper": fileinfo["lowend_peaks"][frag_types][ions]["peak"]["extended"]["three_sigma_ppm_upper"],
                                        "fit_mz": fileinfo["lowend_peaks"][frag_types][ions]['peak']["fit"]['mz'],
                                        "intensity": fileinfo["lowend_peaks"][frag_types][ions]['peak']["fit"]['intensity'],
                                    })

                                elif self.metadata['files'][file]['spectra_stats']['fragmentation_type'].startswith('LR'):
                                    ion_three_sigma_table[file].append({
                                        "ion": ions,
                                        "three_sigma_lower": fileinfo["lowend_peaks"][frag_types][ions]["peak"]["extended"]["three_sigma_mz_lower"],
                                        "three_sigma_upper": fileinfo["lowend_peaks"][frag_types][ions]["peak"]["extended"]["three_sigma_mz_upper"],
                                        "fit_mz": fileinfo["lowend_peaks"][frag_types][ions]['peak']["fit"]['mz'],
                                        "intensity": fileinfo["lowend_peaks"][frag_types][ions]['peak']["fit"]['intensity'],
                                    })              
                            except:
                                pass
            except:
                if file not in ion_three_sigma_table:
                    ion_three_sigma_table[file] = []
                ion_three_sigma_table[file].append("No peaks")
            #### Gather info for summary table 
            try:
                    
                info_dict = {
                    "file": file,
                    "acquisition datetime": fileinfo.get('start_timestamp', ''),
                    "relative acquisition time (days)": fileinfo.get('delta_start_days', '?'),
                    "labeling": labeling,
                    "file_instrument": file_instrument,
                    "acquisition type": fileinfo['spectra_stats'].get('acquisition_type', ''),
                    "High accuracy precursor": high_accuracy_precursors,
                }
                # Fragmentation type
                try:
                    frag_type = fileinfo['summary']['combined summary']['fragmentation type']
                except:
                    frag_type = "N/A"
                info_dict["fragmentation type"] = frag_type

                # Fragmentation tolerances
                try:
                    #if 'warning' not in frag_tol:
                    if self.fragmentation_units in {"mz", "ppm"}:
                        units = self.fragmentation_units
                    else:
                        frag_type = self.metadata['files'][file]['spectra_stats']['fragmentation_type']
                        if frag_type.startswith('HR'):
                            units = "ppm"
                        elif frag_type.startswith('LR'):
                            units = "mz"
                        else:
                            units = None
                    if units:
                        if units == "ppm":
                            lower_key, upper_key = ('fragment_tolerance_ppm_lower','fragment_tolerance_ppm_upper')
                            error = self.ppm_error
                            unit_label = "ppm"
                            rec_fmt = lambda v: math.ceil(v)
                            frag_tol = fileinfo['summary']['combined summary']['fragmentation tolerance']
                            low_tol = frag_tol[lower_key]
                            high_tol = frag_tol[upper_key]
                        else: 
                            lower_key, upper_key = ('fragment_tolerance_mz_lower', 'fragment_tolerance_mz_upper')
                            error = self.mz_error
                            unit_label = "m/z"
                            rec_fmt = lambda v: round(v, 4)
                            frag_tol = fileinfo['summary']['combined summary']['fragmentation tolerance']
                            low_tol = rec_fmt(frag_tol[lower_key])
                            high_tol = rec_fmt(frag_tol[upper_key])
                            
                        combined = math.sqrt(error**2 + max(abs(low_tol), abs(high_tol))**2)
                        rec_tol = rec_fmt(combined)
                        frag_units = unit_label
                        
                    else:
                        low_tol = high_tol = rec_tol = "multiple"
                        frag_units = None


                except (KeyError, TypeError):
                    if units == "ppm":
                        low_tol = high_tol = "N/A"
                        rec_tol = 20
                        frag_units = "ppm"
                    else:
                        low_tol = high_tol = "N/A"
                        rec_tol = 0.6
                        frag_units = "m/z" 

                info_dict["fragment tolerance lower_three_sigma"] = low_tol
                info_dict["fragment tolerance lower_three_sigma units"] = frag_units
                info_dict["fragment tolerance upper_three_sigma"] = high_tol
                info_dict["fragment tolerance upper_three_sigma units"] = frag_units
                info_dict["recommended fragment tolerance"] = rec_tol
                info_dict["recommended fragment tolerance units"] = frag_units

                
                # Dynamic exclusion time
                try:
                    dynamic_exclusion_time = round(fileinfo['summary']['precursor stats']['dynamic exclusion window']['fit_pulse_time']['pulse start'], 2)
                except:
                    dynamic_exclusion_time = "N/A"
                info_dict["dynamic exclusion time (s)"] = dynamic_exclusion_time

                # Precursor tolerances
                try:
                    precursor_low = round(fileinfo['summary']['precursor stats']['precursor tolerance']['fit_ppm']['lower_three_sigma (ppm)'], 2)
                    precursor_high = round(fileinfo['summary']['precursor stats']['precursor tolerance']['fit_ppm']['upper_three_sigma (ppm)'], 2)
                    recommended_precursor = fileinfo['summary']['precursor stats']['precursor tolerance']['fit_ppm']['recommended precursor tolerance (ppm)']
                    if recommended_precursor == 0:
                        recommended_precursor = "N/A"

                except:
                    precursor_low = precursor_high = recommended_precursor =  "N/A"
                info_dict["precursor tolerance three_sigma_lower (ppm)"] = precursor_low
                info_dict["precursor tolerance three_sigma_higher (ppm)"] = precursor_high
                info_dict['recommended precursor tolerance (ppm)'] = recommended_precursor

                # Isolation window
                try:
                    isolation_window = fileinfo['spectra_stats']['isolation_window_full_widths']
                    
                    if isinstance(isolation_window, dict):
                        if len(isolation_window) <= 0:
                            iso_str = "N/A"
                        elif len(isolation_window) > 3:
                            first_three = list(isolation_window.items())[:3]
                            iso_str = ', '.join(f"{{{k}: {v}}}" for k, v in first_three) + " ..."
                        else:
                            iso_str = ', '.join(f"{{{k}: {v}}}" for k, v in isolation_window.items())
                    else:
                        iso_str = "N/A"
                except:
                    iso_str = "N/A"
                info_dict["isolation window"] = iso_str

                # Water/phospho info
                try:
                    has_water = fileinfo['summary']['combined summary']['has water_loss']
                except:
                    has_water = False
                info_dict["has water_loss"] = has_water

                try:
                    has_phospho = fileinfo['summary']['combined summary']['has phospho_spectra']
                except:
                    has_phospho = False
                info_dict["has phospho_spectra"] = has_phospho
            
                if has_phospho:
                    phosphoenrichment['has_phospho'] += 1
                else:
                    phosphoenrichment['no_phospho'] += 1


                try:
                    total_phospho = fileinfo['summary']['combined summary']['total intensity of z=2 phosphoric_acid']
                    total_water = fileinfo['summary']['combined summary']['total intensity of z=2 water_loss']
                    ratio = round(fileinfo['summary']['combined summary']['total z=2 phosphoric_acid to z=2 water_loss intensity ratio'], 2)
                except:
                    total_phospho = total_water = ratio = "N/A"
                info_dict["total intensity of z=2 phospho_spectra"] = total_phospho
                info_dict["total intensity of z=2 water_loss_spectra"] = total_water
                info_dict["total z=2 phosphoric_acid to z=2 water_loss intensity ratio"] = ratio

                #Info about phospho shifts
                for keys in fileinfo['summary']:
                    if "HR" in keys:
                        try:
                            info_dict[keys + " absolute difference between delta m/z for z=2"] = round(fileinfo['summary'][keys]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'], 5)
                        except:
                            info_dict[keys + " absolute difference between delta m/z for z=2"] =  "N/A"
                    elif "LR" in keys:
                        try:
                            info_dict[keys + " absolute difference between delta m/z for z=2"] = round(fileinfo['summary'][keys]['absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss'], 2)
                        except:
                            info_dict[keys + " absolute difference between delta m/z for z=2"] =  "N/A"

                # Save
                info.append(info_dict)

            except:
               info.append({"file": file, "labeling": "N/A"})

        #### Write file summary table
        self.write_summary_table(info)

        #### Set the the main tolerance based on three_sigma values from all files
        self.set_main_tolerance(all_3sigma_values_away)

        ### Set phosphoenrichment
        criteria.setdefault('phosphoenrichment', {})
        if phosphoenrichment['has_phospho'] > 0 and phosphoenrichment['no_phospho'] == 0:
            criteria['phosphoenrichment']['has_phosphoenrichment'] = 'true'
            criteria['phosphoenrichment']['phospho_enriched_spectra'] = phosphoenrichment['has_phospho']
        elif phosphoenrichment['has_phospho'] == 0 and phosphoenrichment['no_phospho'] > 0:
            criteria['phosphoenrichment']['has_phosphoenrichment'] = 'false'
            criteria['phosphoenrichment']['non_phospho_enriched_spectra'] = phosphoenrichment['no_phospho']
        elif phosphoenrichment['has_phospho'] > 0 and phosphoenrichment['no_phospho'] > 0:
            criteria['phosphoenrichment']['has_phosphoenrichment'] = 'mixed'
            criteria['phosphoenrichment']['phospho_enriched_spectra'] = phosphoenrichment['has_phospho']
            criteria['phosphoenrichment']['non_phospho_enriched_spectra'] = phosphoenrichment['no_phospho']
        else:
            criteria['phosphoenrichment']['has_phosphoenrichment'] = 'unknown'
            criteria['phosphoenrichment']['phospho_enriched_spectra'] = phosphoenrichment['has_phospho']
            criteria['phosphoenrichment']['non_phospho_enriched_spectra'] = phosphoenrichment['no_phospho']

        #### write ion data into a table
        if self.write_ions == True:
            self.write_ion_table(ion_three_sigma_table)
            eprint("INFO: Ion Three-Sigma Table generated and stored")

    ####################################################################################################
    #### If standard deviations have been found, set them in the file
    def set_main_tolerance(self, three_sigma_dict):
        criteria = self.metadata['search_criteria']
        ppm_error = 5
        criteria.setdefault('tolerances', {})
        percentile = 90 #Picks this percentile for the three_sigma value
        lower = three_sigma_dict['Lowest']
        upper = three_sigma_dict['Highest']       
        criteria.setdefault('tolerances', {})
        
        
        if (three_sigma_dict['Status']):
            three_sigma_upper = numpy.percentile(upper, percentile) # Gives x percentile for a sorted list
            three_sigma_lower = numpy.percentile(lower, 100-percentile)

            if criteria["fragmentation_type"].startswith('HR'):
                self.metadata['search_criteria']['tolerances']['overall_upper_fragment_tolerance_ppm'] = three_sigma_upper
                self.metadata['search_criteria']['tolerances']['overall_lower_fragment_tolerance_ppm'] = three_sigma_lower
                self.metadata['search_criteria']['tolerances']['recommended overall fragment tolerance (ppm)'] = math.ceil(math.sqrt(ppm_error**2 + (max(abs(three_sigma_lower), abs(three_sigma_upper))**2)))

            elif criteria["fragmentation_type"].startswith('LR'):
                self.metadata['search_criteria']['tolerances']['overall_lower_m/z'] =three_sigma_lower
                self.metadata['search_criteria']['tolerances']['overall_upper_m/z'] = three_sigma_upper
                self.metadata['search_criteria']['tolerances']['recommended overall fragment tolerance (m/z)'] = math.ceil(math.sqrt(ppm_error**2 + (max(abs(three_sigma_lower), abs(three_sigma_upper))**2)))

            elif criteria['fragmentation_type'] == "multiple":
                self.metadata['search_criteria']['tolerances']['recommended overall fragment tolerance'] = "multiple"

        if three_sigma_dict['precursor_found']:
            precursor = three_sigma_dict['precursors']
            upper_precursor = numpy.percentile(precursor, percentile)
            self.metadata['search_criteria']['tolerances']['recommended overall precursor tolerance (ppm)'] = upper_precursor
        else:
            self.metadata['search_criteria']['tolerances']['recommended overall precursor tolerance (ppm)'] = None


    ###################################################################################################
    #### Generates a table of every file and their fit ions with upper and lower three_sigma values
    def write_ion_table(self, ion_data):
        rows = []
        for file, ions_list in ion_data.items():            
                for ion_info in ions_list:
                    if ion_info == "No peaks":
                        pass
                    else:
                        rows.append({
                            "file": file,
                            "ion type": ion_info["ion"],
                            "three_sigma_ppm_lower": ion_info["three_sigma_lower"],
                            "three_sigma_ppm_upper": ion_info["three_sigma_upper"],
                            "fit_mz": ion_info['fit_mz'],
                            "fit_intensity": ion_info['intensity']
                        })

        df = pd.DataFrame(rows)
        ion_file = self.metadata_filepath.replace(".json", ".ion_summary.tsv")
        df.to_csv(ion_file, sep="\t", index=False)


    ####################################################################################################
    #### Generates a table summarizing all the info in 
    def write_summary_table(self, info):
        headers = ["file", "acquisition datetime", "relative acquisition time (days)", "acquisition type", "file_instrument", "labeling",
                    "High accuracy precursor", "fragmentation type", 
                    "fragment tolerance lower_three_sigma", 'fragment tolerance lower_three_sigma units', "fragment tolerance upper_three_sigma", 'fragment tolerance upper_three_sigma units', 
                    "recommended fragment tolerance", "recommended fragment tolerance units", 
                    "dynamic exclusion time (s)",
                    "precursor tolerance three_sigma_lower (ppm)", "precursor tolerance three_sigma_higher (ppm)", 
                    "recommended precursor tolerance (ppm)",
                    "isolation window", 
                    "has water_loss", "has phospho_spectra",
                    "total intensity of z=2 phospho_spectra", "total intensity of z=2 water_loss_spectra",
                    "total z=2 phosphoric_acid to z=2 water_loss intensity ratio"]
        for row in info:
            for key in row:
                if key not in headers:
                    headers.append(key)

        summary_file = self.metadata_filepath.replace(".json", ".file_summary.tsv")
        with open(summary_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t", extrasaction='ignore')
            writer.writeheader()
            writer.writerows(info)
       

    ####################################################################################################
    #### Generate SDRF table data
    def generate_sdrf_table(self, include_provenance=None):

        self.sdrf_table_column_titles = []
        self.sdrf_table_rows = []
        keys_dict = self.sdrf_hints.get('keys', {})

        if len(keys_dict) == 0:
            if self.verbose > 0:
                eprint(f"INFO: Skip generating an SDRF file. Study metadata txt template is not available")
            return

        sdrf_hints_to_skip = [ 'dataset_identifier', 'experiment_identifier', 'sdrf_template' ]

        for key in self.sdrf_hints['key_order']:
            if key in sdrf_hints_to_skip:
                continue
            if key not in self.sdrf_hints['keys']:
                eprint("ERROR E362")
                return
            for value in self.sdrf_hints['keys'][key]:
                self.sdrf_table_column_titles.append(key)

        files = list(self.metadata['files'])
        files.sort()
        i_sample = 1

        for file in files:
            row = []
            for key in self.sdrf_hints['key_order']:
                if key in sdrf_hints_to_skip:
                    continue
                if key not in self.sdrf_hints['keys']:
                    eprint("ERROR E362")
                    return
                for value in self.sdrf_hints['keys'][key]:

                    if value != '' and include_provenance:
                        value += '[HUMAN]'

                    if key == 'source name':
                        value = f"sample {i_sample}"
                        if include_provenance:
                            value += '[AG-TBA]'

                    if key == 'comment[proteomics data acquisition method]':
                        if value == '':
                            if self.metadata['files'][file]['spectra_stats']['acquisition_type'] == 'DDA':
                                value = 'NT=data-dependent acquisition;AC=MS:1003221'
                            elif self.metadata['files'][file]['spectra_stats']['acquisition_type'] == 'DIA':
                                value = 'NT=data-independent acquisition;AC=MS:1003215'
                            if include_provenance:
                                value += '[DATA-I]'

                    if key == 'comment[data file]':
                        value = f"{file}"
                        if include_provenance:
                            value += '[DATA]'

                    if key == 'comment[precursor mass tolerance]':
                        if value == '':
                            if self.metadata['files'][file]['spectra_stats']['high_accuracy_precursors'] == 'true' and self.metadata['files'][file]['summary']['combined summary']['recommended precursor tolerance (ppm)'] != None:
                                value = str(self.metadata['files'][file]['summary']['combined summary']['recommended precursor tolerance (ppm)']) + " ppm"
                            else:
                                value = '3.1 Da'
                            if include_provenance:
                                value += '[DATA-I]'

                    if key == 'comment[fragment mass tolerance]':
                        if value == '':
                            if isinstance(self.metadata['files'][file]['summary']['combined summary']['fragmentation tolerance'], dict):
                                    if self.metadata['files'][file]['summary']['combined summary']['fragmentation tolerance']['recommended fragment tolerance units'] == 'mz':
                                        units = 'Da'
                                        rec_fmt = lambda v: round(v, 2)
                                    else:
                                        units = 'ppm'
                                        rec_fmt = lambda v: round(v, 4)
                                    value = f"{rec_fmt(self.metadata['files'][file]['summary']['combined summary']['fragmentation tolerance']['recommended fragment tolerance'])} {units}"
                            else:
                                if self.metadata['files'][file]['spectra_stats']['fragmentation_type'].startswith('LR'):
                                    value = '0.6 Da'
                                else: 
                                    #Uses this value for HR or multiple fragmentation type
                                    value = '20 ppm'

                    if key == 'comment[instrument]':
                        value = f"NT={self.metadata['files'][file]['instrument_model']['name']};AC={self.metadata['files'][file]['instrument_model']['accession']}"
                        if include_provenance:
                            value += '[DATA]'

                    if key == 'comment[tool metadata]':
                        value = 'RunAssessor v0.1'
                        if include_provenance:
                            value += '[TOOL]'

                    row.append(value)
            self.sdrf_table_rows.append(row)
            i_sample += 1



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
    argparser = argparse.ArgumentParser(description='Handler for study/experiment metadata files')
    argparser.add_argument('--metadata_filepath', action='store', help='Filepath of the metadata file (defaults to study_metadata.json)')
    argparser.add_argument('--verbose', action='count' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    params = argparser.parse_args()

    #### Create metadata handler object and read or create the metadata structure
    metadata = MetadataHandler(params.metadata_filepath, params.verbose)
    metadata.read_or_create()

    #### Try to infer search criteria based on the information we have available
    metadata.infer_search_criteria()

    #### Store the metadata structure, infer the SDRF file name and write the SDRF.tsv file
    metadata.store()
    filename = metadata.infer_filename()
    metadata.write_sdrf_file(filename)

if __name__ == "__main__": main()
