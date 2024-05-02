#!/usr/bin/env python3
import sys
import os
import argparse
import os.path
import re
import json
import numpy
import urllib.parse
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

#sys.path.append("G:\Repositories\GitHub\SpectralLibraries\lib")
#from universal_spectrum_identifier import UniversalSpectrumIdentifier

# Define column offsets for peak_list
PL_I_PEAK = 0
PL_MZ = 1
PL_INTENSITY = 2
PL_INTERPRETATION_STRING = 3
PL_AGGREGATION_INFO = 4
PL_INTERPRETATIONS = 5
PL_ATTRIBUTES = 6
# Define column offsets for peak_list attributes
PLA_CHARGE = 0
PLA_N_ISOTOPES = 1
PLA_IS_ISOTOPE = 2
PLA_DELTA_PPM = 3
PLA_PARENT_PEAK = 4
PLA_N_NEUTRAL_LOSSES = 5
PLA_IS_NEUTRAL_LOSS = 6
PLA_IS_PRECURSOR = 7
PLA_IS_REPORTER = 8
PLA_DIAGNOSTIC_CATEGORY = 9
PLA_IS_DELETED = 10
# Define column offsets for interpretations
INT_MZ = 0
INT_REFERENCE_PEAK = 1
INT_INTERPRETATION_STRING = 2
INT_DELTA_PPM = 3
INT_SCORE = 4
INT_DELTA_SCORE = 5
INT_COMMONNESS_SCORE = 6
INT_DIAGNOSTIC_CATEGORY = 7


####################################################################################################
#### Spectrum class
class Spectrum:

    #### Constructor
    def __init__(self, attribute_list=None, peak_list=None, analyte_list=None):
        """
        __init__ - Spectrum constructor

        Parameters
        ----------
        attribute_list: list
            A list of attribute [key, value (, group)] sets to initialize to.
        """
        if attribute_list is None:
            attribute_list = []
        if peak_list is None:
            peak_list = []
        if analyte_list is None:
            analyte_list = []
        self.attribute_list = attribute_list
        self.peak_list = peak_list
        self.analyte_list = analyte_list

        self.attributes = {}
        self.analytes = { '1': {} }
        self.peak_index = {}

        self.network = {}
        self.metrics = {}

        self.extended_data = {}

        self.reporter_ions = {}
        self.set_up_reporter_ions()


    ####################################################################################################
    #### Fetch a spectrum from PeptideAtlas given a USI
    def fetch_spectrum(self, usi_string):

        import requests
        import requests_cache
        requests_cache.install_cache('spectrum_cache')

        #usi = UniversalSpectrumIdentifier()
        #usi.parse(usi_string)
        #if usi.is_valid is False:
        #    print("ERROR: Provided USI is not valid")
        #    print(json.dumps(usi.__dict__,sort_keys=True,indent=2))
        #    exit()

        #### URL to fetch a spectrum from PeptideAtlas PROXI client
        url = f"https://peptideatlas.org/api/proxi/v0.1/spectra?resultType=full&usi={urllib.parse.quote_plus(usi_string)}"
        #url = f"https://www.ebi.ac.uk/pride/proxi/archive/v0.1/spectra?resultType=full&usi={usi_string}"

        #### Alternative URL for fetching a predicted spectrum
        #url = f"https://proteomecentral.proteomexchange.org/api/proxi/v0.1/spectra?resultType=full&accession=SEQ2MS&usi={usi_string}"

        if usi_string.startswith('http'):
            url = usi_string

        response_content = requests.get(url)
        status_code = response_content.status_code
        if status_code != 200:
            print("ERROR returned with status "+str(status_code))
            print(response_content)
            exit()

        #### Unpack the response spectrum into internal data
        if usi_string.startswith('http'):
            proxi_spectrum = [ { 'mzs': [], 'intensities': [], 'attributes': [] } ]
            lines = str(response_content.content).split(r'\n')
            for line in lines:
                print(line)
                if line[0] >= '0' and line[0] <= '9':
                    values = line.split(r'\t')
                    proxi_spectrum[0]['mzs'].append(values[0])
                    proxi_spectrum[0]['intensities'].append(values[1])
                else:
                    match = re.match(r'CHARGE=(\d+)',line)
                    if match:
                        attribute = { 'accession': 'MS:1000041', 'name': 'charge state', 'value': match.group(1) }
                        proxi_spectrum[0]['attributes'].append(attribute)
                    match = re.match(r'PEPMASS ([\d\.]+)',line)
                    if match:
                        attribute = { 'accession': 'MS:1000827', 'name': 'isolation window target m/z', 'value': match.group(1) }
                        proxi_spectrum[0]['attributes'].append(attribute)

        else:
            proxi_spectrum = response_content.json()


        self.peak_list = []
        n_peaks = len(proxi_spectrum[0]['mzs'])
        for i_peak in range(n_peaks):

            # mz should always be there
            mz = float(proxi_spectrum[0]['mzs'][i_peak])

            # Extract the intensity value for the peak if present
            intensity = 1.0
            if 'intensities' in proxi_spectrum[0]:
                intensity = float(proxi_spectrum[0]['intensities'][i_peak])

            # Extract the interpretation_string value for the peak if present
            interpretation_string = '?'
            if 'interpretations' in proxi_spectrum[0] and proxi_spectrum[0]['interpretations'] is not None:
                interpretation_string = proxi_spectrum[0]['interpretations'][i_peak]

            # Current PROXI spectrum doesn't have aggregation information, interpretations as a list or peak attributes,
            # so just set those to empty
            aggregation_info = ''
            interpretations = []
            attributes = [ 0, 0, 0, 0, -1, 0, 0, 0, 0, 'unexplained', False ]

            # Store the peak data as a data list in the peak_list
            self.peak_list.append( [ i_peak, mz, intensity, interpretation_string, aggregation_info, interpretations, attributes ] )

        # Extract the attribute list from the fetched JSON
        self.attribute_list = []
        if 'attributes' in proxi_spectrum[0]:
            self.attribute_list = proxi_spectrum[0]['attributes']
            for attribute in self.attribute_list:
                if attribute['accession'] == 'MS:1000827' or attribute['name'] == 'isolation window target m/z':
                    self.attributes['isolation window target m/z'] = float(attribute['value'])
                    self.analytes['1']['precursor_mz'] = float(attribute['value'])
                if attribute['accession'] == 'MS:1000041' or attribute['name'] == 'charge state':
                    self.analytes['1']['charge state'] = int(attribute['value'])

        # Add a few attributes by key
        self.attributes['usi'] = usi_string
        self.attributes['number of peaks'] = n_peaks
        self.attributes['n_identified_peptide_low_mass_ions'] = 0
        self.attributes['n_identified_reporter_ions'] = 0
        self.attributes['mass_accuracy'] = {
            'offset': 0.0,
            'siqr': 4.0,
            'is_optimized': False,
            'best_tolerance': 4.0,
            'middle_tolerance': 8.0,
            'outer_tolerance': 12.0,
#            'best_tolerance': 19.0,
#            'middle_tolerance': 19.1,
#            'outer_tolerance': 19.2,
            'max_tolerance': 20.0,
        }

        return self


    ####################################################################################################
    #### Fill a spectrum object with a basic set of data elements
    def fill(self, mzs=None, intensities=None, interpretations=None, precursor_mz=None, charge_state=None, usi_string=None):

        #### Check the input
        if mzs is None or intensities is None:
            eprint(f"ERROR: Spectrum.fill is missing mzs or intensities")
            return
        if len(mzs) != len(intensities) is None:
            eprint(f"ERROR: Spectrum.fill: Length of mzs does not match length of intensities")
            return

        #### Fill the peak list in the correct format
        self.peak_list = []
        n_peaks = len(mzs)
        for i_peak in range(n_peaks):

            # mz should always be there
            mz = float(mzs[i_peak])

            # Extract the intensity value for the peak
            intensity = float(intensities[i_peak])

            # Extract the interpretation_string value for the peak if present
            interpretation_string = '?'
            if interpretations is not None:
                try:
                    interpretation_string = interpretations[i_peak]
                except:
                    pass

            # Current PROXI spectrum doesn't have aggregation information, interpretations as a list or peak attributes,
            # so just set those to empty
            aggregation_info = ''
            peak_interpretations = []
            attributes = [ 0, 0, 0, 0, -1, 0, 0, 0, 0, 'unexplained', False ]

            # Store the peak data as a data list in the peak_list
            self.peak_list.append( [ i_peak, mz, intensity, interpretation_string, aggregation_info, peak_interpretations, attributes ] )

        # Set up the attribute list from the proided information
        self.attribute_list = []
        if precursor_mz is not None:
            self.attributes['isolation window target m/z'] = float(precursor_mz)
            self.analytes['1']['precursor_mz'] = float(precursor_mz)
        if charge_state is not None:
            self.analytes['1']['charge state'] = int(charge_state)

        # Add a few attributes by key
        self.attributes['usi'] = usi_string
        self.attributes['number of peaks'] = n_peaks
        self.attributes['n_identified_peptide_low_mass_ions'] = 0
        self.attributes['n_identified_reporter_ions'] = 0
        self.attributes['mass_accuracy'] = {
            'offset': 0.0,
            'siqr': 4.0,
            'is_optimized': False,
            'best_tolerance': 4.0,
            'middle_tolerance': 8.0,
            'outer_tolerance': 12.0,
            'max_tolerance': 20.0,
        }

        return self


    ####################################################################################################
    #### Compute a set of metrics from the spectrum such as SNR and quality score
    def compute_spectrum_metrics(self):

        spectrum = self
        n_peaks = spectrum.attributes['number of peaks']
        precursor_mz = None
        selection_mz = None
        selection_floor_mz = None
        selection_ceiling_mz = None

        selection_half_width = 1.5          #### FIXME hardcoded for the moment
        isobaric_label_mode = 'TMTpro'      #### FIXME hardcoded to TMTpro for the moment

        if 'isolation window target m/z' in spectrum.attributes:
            precursor_mz = spectrum.attributes['isolation window target m/z']
            selection_mz = precursor_mz
            selection_floor_mz = selection_mz - selection_half_width
            selection_ceiling_mz = selection_mz + selection_half_width

        charge = None
        if 'charge state' in spectrum.analytes['1']:
            charge = spectrum.analytes['1']['charge state']

        #### TMT spectra tend to have a lot of junk within p^1 - 180, so just throw all that out
        #### Probably should only do this in TMT mode, but hard-coded for now. FIXME
        charge1_precursor_losses_floor = None
        reporter_ion_floor = None
        reporter_ion_ceiling = None
        reporter_ion_intensities = None
        potential_reporter_peaks = {}
        if isobaric_label_mode.startswith('TMT') and precursor_mz is not None and charge is not None:
            charge1_precursor_losses_floor = precursor_mz * charge - ( charge - 1) - 180
            reporter_ion_floor = 126
            reporter_ion_ceiling = 135.5

        intensity_list = []
        mz_list = []

        #### Loop through the peaks and create lists
        for i_peak in range(n_peaks):

            #### Add to the list of all peaks
            mz = spectrum.peak_list[i_peak][PL_MZ]
            intensity = spectrum.peak_list[i_peak][PL_INTENSITY]
            mz_list.append(mz)
            intensity_list.append(intensity)

            #### Add to the special list of potential reporter peaks for later matching
            if reporter_ion_floor is not None and mz < reporter_ion_ceiling and mz > reporter_ion_floor:
                int_mz = int(mz)
                if int_mz not in potential_reporter_peaks:
                    potential_reporter_peaks[int_mz] = []
                potential_reporter_peaks[int_mz].append( { 'mz': mz, 'intensity': intensity })

        #### If we're processing reporter ions, scan through the list of desired reporter ions and record the intensities
        if reporter_ion_floor is not None:
            reporter_ion_tolerance = 10.0
            reporter_ion_intensities = [None] * len(self.reporter_ions[isobaric_label_mode])
            i_reporter_ion = 0
            for reporter_ion in self.reporter_ions[isobaric_label_mode]:
                int_mz = int(reporter_ion['mz'])
                if int_mz in potential_reporter_peaks:
                    for potential_reporter_peak in potential_reporter_peaks[int_mz]:
                        if abs(potential_reporter_peak['mz'] - reporter_ion['mz']) / int_mz * 1e6 < reporter_ion_tolerance:
                            reporter_ion_intensities[i_reporter_ion] = potential_reporter_peak['intensity']
                i_reporter_ion += 1


        #### Convert ms_list to numpy array and compute metrics
        mz_array = numpy.asarray(mz_list)
        sorted_array = numpy.sort(mz_array)
        spectrum.attributes['minimum mz'] = sorted_array[0]
        spectrum.attributes['maximum mz'] = sorted_array[-1]

        #### Convert intensities list to numpy array and compute metrics
        intensity_array = numpy.asarray(intensity_list)
        sorted_array = numpy.sort(intensity_array)
        if len(intensity_list) >= 10:
            noise = sorted_array[3]
            signal = sorted_array[-3]
        else:
            index = int(len(intensity_list)/5)
            noise = sorted_array[index]
            signal = sorted_array[-1*index]

        spectrum.attributes['minimum intensity'] = sorted_array[0]
        spectrum.attributes['maximum intensity'] = sorted_array[-1]
        if spectrum.attributes['maximum intensity'] == 0:
            spectrum.attributes['maximum intensity'] = 0.01

        #### Prevent division by 0
        if noise == 0:
            noise = 0.03

        # Assume that the smallest peaks are signal to noise of 3.0. An arbitrary guess
        # At least in HCD spectra, the smallest peaks are usually not true noise but can be real signal.
        # They're the smallest peaks above the noise that the peak-finding algorithm selected presumably
        # Which might be around S/N of 3 maybe?
        noise = noise / 3.0
        signal_to_noise = signal / noise

        spectrum.attributes['noise level'] = noise
        spectrum.attributes['signal level'] = signal
        spectrum.attributes['signal to noise level'] = signal_to_noise

        #### Now that we know the noise level, compile per-bin S/N stats
        stats_per_bin = {}
        bin_size = 100
        ibin = 0
        bin_floor = 0
        bin_ceiling = bin_floor + bin_size

        informative_peak_intensities_list = []

        for i_peak in range(n_peaks):

            #### All to the list of all peaks
            mz = spectrum.peak_list[i_peak][PL_MZ]
            intensity = spectrum.peak_list[i_peak][PL_INTENSITY]

            if selection_floor_mz is not None and selection_ceiling_mz is not None:
                if precursor_mz is not None and mz > selection_floor_mz and mz < selection_ceiling_mz:
                    continue
            if charge1_precursor_losses_floor is not None and mz > charge1_precursor_losses_floor:
                continue

            #### Add to the binned information stats
            if mz >= bin_ceiling:
                bin_floor += bin_size
                bin_ceiling += bin_size
                ibin += 1
            if str(ibin) not in stats_per_bin:
                stats_per_bin[str(ibin)] = { 'bin_floor': bin_floor, 'bin_ceiling': bin_ceiling, 'n_peaks': 0,
                                       'n_peaks_sn5': 0, 'n_peaks_sn10': 0, 'n_peaks_sn20': 0 }
            stats_per_bin[str(ibin)]['n_peaks'] += 1
            for sn_factor in [ 5, 10, 20 ]:
                if intensity > noise * sn_factor:
                    stats_per_bin[str(ibin)][f"n_peaks_sn{sn_factor}"] += 1

            #### Add the reporter ion information stats
            if reporter_ion_floor is not None and mz < reporter_ion_ceiling and mz > reporter_ion_floor:
                bin_name = 'reporter_ions'
                if bin_name not in stats_per_bin:
                    stats_per_bin[bin_name] = { 'bin_floor': reporter_ion_floor, 'bin_ceiling': reporter_ion_ceiling, 'n_peaks': 0,
                                        'n_peaks_sn10': 0, 'n_peaks_sn30': 0, 'n_peaks_sn100': 0, 'snrs': [], 'intensities': reporter_ion_intensities }
                stats_per_bin[bin_name]['n_peaks'] += 1
                stats_per_bin[bin_name]['snrs'].append(intensity / noise)
                for sn_factor in [ 10, 30, 100 ]:
                    if intensity > noise * sn_factor:
                        stats_per_bin[bin_name][f"n_peaks_sn{sn_factor}"] += 1

            if ibin >= 2:
                informative_peak_intensities_list.append(intensity)

        if 'reporter_ions' in stats_per_bin:
            n_peaks = stats_per_bin['reporter_ions']['n_peaks']
            median_position = int(n_peaks / 2)
            snrs_array = numpy.asarray(stats_per_bin[bin_name]['snrs'])
            sorted_snrs_array = numpy.sort(snrs_array)
            stats_per_bin['reporter_ions']['median_snr'] = sorted_snrs_array[median_position]

        spectrum.metrics['stats_per_bin'] = stats_per_bin

        informative_peak_intensities_array = numpy.asarray(informative_peak_intensities_list)
        sorted_informative_peak_intensities_array = numpy.sort(informative_peak_intensities_array)
        if len(informative_peak_intensities_array) == 0:
            weighted_snr = 0.0
        elif len(informative_peak_intensities_array) < 5:
            weighted_snr = sorted_informative_peak_intensities_array[0] / noise
        elif len(informative_peak_intensities_array) < 12:
            weighted_snr = sorted_informative_peak_intensities_array[3] / noise
        else:
            weighted_snr = sorted_informative_peak_intensities_array[-10] / noise

        #### Compute a spectrum quality score
        start_bin = 2
        mz_ceiling = spectrum.attributes['maximum mz']
        if charge1_precursor_losses_floor is not None:
            mz_ceiling = charge1_precursor_losses_floor
        end_bin= int(mz_ceiling / 100.0) - 1
        if end_bin < 6:
            end_bin = 6
        if end_bin > 20:
            end_bin = 20
        ibin = start_bin
        quality_score = 0.0
        while ibin <= end_bin:
            str_ibin = str(ibin)
            bin_score = 0.0
            if str_ibin in stats_per_bin:
                n_great_peaks = stats_per_bin[str_ibin]['n_peaks_sn20']
                n_good_peaks = stats_per_bin[str_ibin]['n_peaks_sn10'] - n_great_peaks
                n_okay_peaks = stats_per_bin[str_ibin]['n_peaks_sn5'] - n_great_peaks - n_good_peaks
                if n_great_peaks >= 4:
                    bin_score = 1.0
                elif n_great_peaks + n_good_peaks >= 4:
                    bin_score = 0.5 + n_great_peaks * 0.1
                elif n_great_peaks + n_good_peaks + n_okay_peaks >= 4:
                    bin_score = 0.2 + n_great_peaks * 0.1 + n_good_peaks * 0.05
                else:
                    bin_score = n_great_peaks * 0.1 + n_good_peaks * 0.05 + n_okay_peaks * 0.025
            quality_score += bin_score
            #print(f"** For ibin {ibin}, bin_score = {bin_score}")
            ibin += 1
        quality_score = quality_score / (ibin - 2) * 100
        #print(f"** Final score = {quality_score}")

        spectrum.metrics['quality_score'] = quality_score
        spectrum.metrics['weighted_snr'] = weighted_snr
        spectrum.metrics['stats_per_bin'] = stats_per_bin
        spectrum.metrics['isobaric_label_mode'] = isobaric_label_mode
        spectrum.metrics['precursor_mz'] = precursor_mz
        spectrum.metrics['selection_mz'] = selection_mz
        spectrum.metrics['selection_half_width'] = selection_half_width
        spectrum.metrics['selection_floor_mz'] = selection_floor_mz
        spectrum.metrics['selection_ceiling_mz'] = selection_ceiling_mz
        spectrum.metrics['charge1_precursor_losses_floor'] = charge1_precursor_losses_floor
        spectrum.metrics['reporter_ion_floor'] = reporter_ion_floor
        spectrum.metrics['reporter_ion_ceiling'] = reporter_ion_ceiling
        if 'reporter_ions' in stats_per_bin:
            spectrum.metrics['reporter_ion_intensities'] = stats_per_bin['reporter_ions']['intensities']
            #print(spectrum.metrics['reporter_ion_intensities'])



    ####################################################################################################
    #### Return a printable buffer string of the details of the peptidoform
    def show(self, show_all_annotations=False, verbose=0):

        buf = ''

        if verbose > 0:

            if 'peptidoform sequence' in self.analytes['1']:
                buf += f"peptidoform sequence={self.analytes['1']['peptidoform sequence']}\n"

            if 'charge state' in self.analytes['1']:
                buf += f"charge state={self.analytes['1']['charge state']}\n"

            for attribute_name, value in self.attributes.items():
                if attribute_name in [ 'psm_score', 'mass_accuracy' ]:
                    continue
                buf += f"  - {attribute_name}={value}\n"

            for attribute_name in [ 'psm_score', 'mass_accuracy' ]:
                if attribute_name in self.attributes:
                    buf += f"  - {attribute_name} = " + json.dumps(self.attributes[attribute_name], indent=2, sort_keys=True) + "\n"

            if len(self.metrics) > 0:
                buf += f"  - {self.metrics} = " + json.dumps(self.metrics, indent=2, sort_keys=True) + "\n"

        for peak in self.peak_list:
            i_peak = peak[PL_I_PEAK]
            mz = peak[PL_MZ]
            intensity = peak[PL_INTENSITY]
            interpretations_string = peak[PL_INTERPRETATION_STRING]
            diagnostic_category = peak[PL_ATTRIBUTES][PLA_DIAGNOSTIC_CATEGORY]

            is_deleted = peak[PL_ATTRIBUTES][PLA_IS_DELETED]
            if is_deleted:
                continue

            buf += '{:4d}'.format(i_peak) + '{:10.4f}'.format(mz) + '{:10.1f}'.format(intensity) + '  ' + interpretations_string + "\n"
            #buf += '{:4d}'.format(i_peak) + '{:10.4f}'.format(mz) + '{:10.1f}'.format(intensity) + '  ' + interpretations_string + '  ' + diagnostic_category + "\n"

            if show_all_annotations is True:
                attrs = peak[PL_ATTRIBUTES]
                #if attrs[PLA_CHARGE] > 0 or attrs[PLA_N_ISOTOPES] > 0 or attrs[PLA_IS_ISOTOPE] > 0 or attrs[PLA_PARENT_PEAK] > -1:
                if len(peak[PL_INTERPRETATIONS]) > 0:
                    buf += f"                             ** ch={attrs[PLA_CHARGE]}, n_iso={attrs[PLA_N_ISOTOPES]}, is_iso={attrs[PLA_IS_ISOTOPE]}, "
                    buf += f"parent={attrs[PLA_PARENT_PEAK]}, n_NL={attrs[PLA_N_NEUTRAL_LOSSES ]}, is_NL={attrs[PLA_IS_NEUTRAL_LOSS ]}, "
                    buf += f"is_P={attrs[PLA_IS_PRECURSOR]}, is_reporter={attrs[PLA_IS_REPORTER]}\n"
                    for interpretation in peak[PL_INTERPRETATIONS]:
                        diagnostic_category = interpretation[INT_DIAGNOSTIC_CATEGORY]
                        buf += "                               ++ " + '{:5.3f}'.format(interpretation[INT_SCORE])
                        buf += '{:8.1f}'.format(interpretation[INT_DELTA_PPM]) + "  "
                        buf += '{:8.2f}'.format(interpretation[INT_COMMONNESS_SCORE]) + "  " + interpretation[INT_INTERPRETATION_STRING] + "\n"
                        #buf += '{:8.2f}'.format(interpretation[INT_COMMONNESS_SCORE]) + "  " + interpretation[INT_INTERPRETATION_STRING] + '  ' + diagnostic_category + "\n"

        return buf


    ####################################################################################################
    #### Set up reporter ion information
    def set_up_reporter_ions(self):

        self.reporter_ions = {}

        # Define reporter ions to look for
        self.reporter_ions['TMTpro'] = [
            # Jimmy's numbers from https://proteomicsresource.washington.edu/protocols03/isotopic_labeling.php
            { 'name': 'TMT126',  'type': 'TMT', 'mz': 126.127726 },
            { 'name': 'TMT127N', 'type': 'TMT', 'mz': 127.124761 },
            { 'name': 'TMT127C', 'type': 'TMT', 'mz': 127.131081 },
            { 'name': 'TMT128N', 'type': 'TMT', 'mz': 128.128116 },
            { 'name': 'TMT128C', 'type': 'TMT', 'mz': 128.134436 },
            { 'name': 'TMT129N', 'type': 'TMT', 'mz': 129.131471 },
            { 'name': 'TMT129C', 'type': 'TMT', 'mz': 129.137790 },
            { 'name': 'TMT130N', 'type': 'TMT', 'mz': 130.134825 },
            { 'name': 'TMT130C', 'type': 'TMT', 'mz': 130.141145 },
            { 'name': 'TMT131N', 'type': 'TMT', 'mz': 131.138180 },
            { 'name': 'TMT131C', 'type': 'TMT', 'mz': 131.1445 },
            { 'name': 'TMT132N', 'type': 'TMT', 'mz': 132.141535 },
            { 'name': 'TMT132C', 'type': 'TMT', 'mz': 132.147855 },
            { 'name': 'TMT133N', 'type': 'TMT', 'mz': 133.14489 },
            { 'name': 'TMT133C', 'type': 'TMT', 'mz': 133.15121 },
            { 'name': 'TMT134N', 'type': 'TMT', 'mz': 134.148245 },
            { 'name': 'TMT134C', 'type': 'TMT', 'mz': 134.154565 },
            { 'name': 'TMT135N', 'type': 'TMT', 'mz': 135.151600 },
            ]



    ####################################################################################################
    #### Get the mzs, intensities, and interpretations
    def get_peaks(self):

        mzs = []
        intensities = []
        interpretations = []

        for peak in self.peak_list:
            mzs.append(peak[PL_MZ])
            intensities.append(peak[PL_INTENSITY])
            interpretations.append(peak[PL_INTERPRETATION_STRING])

        return mzs, intensities, interpretations


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Class representing a peptidoform')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('usi', type=str, nargs='*', help='Optional USIs to load')
    params = argparser.parse_args()

    # Set verbose mode
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    # If there are supplied peptidoform strings, use those, else set a default example
    if params.usi is not None and len(params.usi) > 0:
        usi = params.usi[0]
    else:
        usi = 'mzspec:PXD007058:SF_200217_pPeptideLibrary_pool1_HCDOT_rep1:scan:4336:GQEY[Phospho]LILEK/2'

    spectrum = Spectrum()
    spectrum.fetch_spectrum(usi)
    spectrum.compute_spectrum_metrics()
    print(spectrum.show())


#### For command line usage
if __name__ == "__main__": main()
