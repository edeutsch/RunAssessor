#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import re
import json

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

        url = f"http://www.peptideatlas.org/api/proxi/v0.1/spectra?resultType=full&usi={usi_string}"

        response_content = requests.get(url)
        status_code = response_content.status_code
        if status_code != 200:
            print("ERROR returned with status "+str(status_code))
            print(response_content)
            exit()

        #### Unpack the response spectrum into internal data
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
            if 'interpretations' in proxi_spectrum[0]:
                interpretation_string = proxi_spectrum[0]['interpretations'][i_peak]

            # Current PROXI spectrum doesn't have aggregation information, interpretations as a list or peak attributes,
            # so just set those to empty
            aggregation_info = ''
            interpretations = []
            attributes = [ 0, 0, 0, 0, -1, 0, 0, 0, 0, 'unexplained' ]

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
            'max_tolerance': 20.0,
        }

        return self


    ####################################################################################################
    #### Return a printable buffer string of the details of the peptidoform
    def show(self, show_all_annotations=False):

        buf = ''
        if 'peptidoform sequence' in self.analytes['1']:
            buf += f"peptidoform sequence={self.analytes['1']['peptidoform sequence']}\n"
        if 'charge state' in self.analytes['1']:
            buf += f"charge state={self.analytes['1']['charge state']}\n"
        for attribute,value in self.attributes.items():
            buf += f"{attribute}={value}\n"
        if 'psm_score' in self.attributes:
            print(json.dumps(self.attributes['psm_score'], indent=2, sort_keys=True))

        for peak in self.peak_list:
            i_peak = peak[PL_I_PEAK]
            mz = peak[PL_MZ]
            intensity = peak[PL_INTENSITY]
            interpretations_string = peak[PL_INTERPRETATION_STRING]
            diagnostic_category = peak[PL_ATTRIBUTES][PLA_DIAGNOSTIC_CATEGORY]
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
    print(spectrum.show())


#### For command line usage
if __name__ == "__main__": main()
