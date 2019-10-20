#!/usr/bin/env python3

#### Define eprint() as printing to stderr
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import timeit
import re
import json
import numpy
import pickle
from lxml import etree

#### Import technical modules and pyteomics
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy import exp
from pyteomics import mzml, auxiliary

#### Import the metadata handler
from metadata_handler import MetadataHandler


####################################################################################################
#### mzML Assessor class
class MzMLAssessor:


    ####################################################################################################
    #### Constructor
    def __init__(self, mzml_file, metadata=None, verbose=None):
        self.mzml_file = mzml_file

        #### Store the provided metadata or create a template to work in
        if metadata is None:
            self.metadata = { 'files': { mzml_file: {} } }
        else:
            self.metadata = metadata
            metadata['files'][mzml_file] = {}

        #### Create a place to store our composite spectra for analysis
        self.composite = {}

        #### Set verbosity
        if verbose is None: verbose = 0
        self.verbose = verbose


    ####################################################################################################
    #### Read spectra
    def read_spectra(self):

        #### Set up information
        t0 = timeit.default_timer()
        stats = { 'n_spectra': 0, 'n_ms1_spectra': 0, 'n_ms2_spectra': 0, 'n_HCD_spectra': 0, 'n_IT_spectra': 0,
            'high_accuracy_precursors': 'unknown', 'fragmentation_type': 'unknown' }
        self.metadata['files'][self.mzml_file]['spectra_stats'] = stats

        #### Show information
        if self.verbose >= 1:
            eprint(f"INFO: Assessing mzML file {self.mzml_file}")
            progress_intro = False

        #### Read spectra from the file
        with mzml.read(self.mzml_file) as reader:
            for spectrum in reader:

                #### Testing. Print the data structure of the first spectrum
                #if stats['n_spectra'] == 0:
                #    auxiliary.print_tree(spectrum)

                #### Set a default spectrum type
                spectrum_type = 'default'

                #### Look for a filter string and parse it
                if 'filter string' in spectrum['scanList']['scan'][0]:
                    filter_string = spectrum['scanList']['scan'][0]['filter string']
                    match = re.search(' ms ',filter_string)

                    #### If this is an MS1 scan, learn what we can from it
                    if match:
                        stats['n_ms1_spectra'] += 1
                        spectrum_type = 'MS1'
                        match = re.search('^(\S+)',filter_string)
                        if match:
                            if match.group(0) == 'FTMS':
                                if stats['high_accuracy_precursors'] == 'unknown':
                                    stats['high_accuracy_precursors'] = 'true'
                                elif stats['high_accuracy_precursors'] == 'true' or stats['high_accuracy_precursors'] == 'mixed':
                                    pass
                                elif stats['high_accuracy_precursors'] == 'false':
                                    stats['high_accuracy_precursors'] = 'mixed'
                                else:
                                    self.log_event('ERROR','InternalStateError',f"high_accuracy_precursors has a strange state '{stats['high_accuracy_precursors']}'")
                                    break
                            elif match.group(0) == 'ITMS':
                                if stats['high_accuracy_precursors'] == 'unknown':
                                    stats['high_accuracy_precursors'] = 'false'
                                elif stats['high_accuracy_precursors'] == 'false' or stats['high_accuracy_precursors'] == 'mixed':
                                    pass
                                elif stats['high_accuracy_precursors'] == 'true':
                                    stats['high_accuracy_precursors'] = 'mixed'
                                else:
                                    self.log_event('ERROR','InternalStateError',f"high_accuracy_precursors has a strange state '{stats['high_accuracy_precursors']}'")
                                    break
                            else:
                                self.log_event('ERROR','UnknownPrecursorScanType',f"Unknown precursor scan type '{match.group(0)}'")
                                break

                    #### Else it's an MSn (n>1) scan so learn what we can from that
                    else:
                        match = re.search(' ms2 ',filter_string)
                        if match:
                            stats['n_ms2_spectra'] += 1
                            match = re.search('@hcd',filter_string)
                            if match:
                                stats['n_HCD_spectra'] += 1
                                spectrum_type = 'HCD'
                            else:
                                match = re.search('@cid',filter_string)
                                if match:
                                    stats['n_IT_spectra'] += 1
                                    spectrum_type = 'IT'
                                else:
                                    self.log_event('ERROR','FragNotInFilterStr',f"Did not find @hcd or @cid in filter string '{filter_string}'")
                                    break
                        else:
                            self.log_event('ERROR','CantParseFilterStr',f"Unable to parse filter string '{filter_string}'")
                            break

                #### If the ms level is greater than 2, fail
                if spectrum['ms level'] > 2:
                    self.log_event('ERROR','MSnTooHigh',f"MS level is greater than we can handle at '{spectrum['ms level']}'")
                    break

                #### If the ms level is 2, then examine it for information
                if spectrum['ms level'] == 2:
                    precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                    self.add_spectrum(spectrum,spectrum_type,precursor_mz)

                    #### Record the fragmentation type for this run
                    #print(stats['fragmentation_type'])
                    if spectrum_type != stats['fragmentation_type']:
                        if stats['fragmentation_type'] == 'unknown':
                            stats['fragmentation_type'] = spectrum_type
                        elif stats['fragmentation_type'] == 'mixed!':
                            pass
                        else:
                            stats['fragmentation_type'] = 'mixed!'
                            self.log_event('ERROR','MixedFragTypes',f"There are multiple fragmentation types in this MS run. Split them.")


                #### Update counters and print progress
                stats['n_spectra'] += 1
                if self.verbose >= 1:
                    if stats['n_spectra']/1000 == int(stats['n_spectra']/1000):
                        if not progress_intro:
                            eprint("INFO: Processing spectra.. ", end='')
                            progress_intro = True
                        eprint(f"{stats['n_spectra']}.. ", end='', flush=True)

        #### Show progress information
        if self.verbose >= 1:
            eprint("")

        #### Write the composite spectrum to a file
        with open('zzcomposite_spectrum.pickle','wb') as outfile:
            pickle.dump(self.composite,outfile)

        #### Print final timing information
        t1 = timeit.default_timer()
        print(f"INFO: Read {stats['n_spectra']} spectra from {self.mzml_file}")
        print(f"INFO: Elapsed time: {t1-t0}")
        print(f"INFO: Processed {stats['n_spectra']/(t1-t0)} spectra per second")


    ####################################################################################################
    #### Add spectrum
    def add_spectrum(self,spectrum,spectrum_type,precursor_mz):
        if 'n_spectra_examined' not in self.metadata['files'][self.mzml_file]['spectra_stats']:
            self.metadata['files'][self.mzml_file]['spectra_stats']['n_spectra_examined'] = 0
        self.metadata['files'][self.mzml_file]['spectra_stats']['n_spectra_examined'] += 1

        destination = f"lowend_{spectrum_type}"
        destination2 = f"neutral_loss_{spectrum_type}"

        if destination not in self.composite:
            if spectrum_type == 'HCD':
                #print(f"INFO: Creating a composite spectrum {destination}")
                minimum = 100
                maximum = 400
                binsize = 0.0002
                array_size = int( (maximum - minimum ) / binsize ) + 1
                self.composite[destination] = { 'minimum': minimum, 'maximum': maximum, 'binsize': binsize }
                self.composite[destination]['intensities'] = numpy.zeros(array_size,dtype=numpy.float32)
                self.composite[destination]['n_peaks'] = numpy.zeros(array_size,dtype=numpy.int32)

                print(f"INFO: Creating a composite spectrum {destination2}")
                minimum = 0
                maximum = 165
                binsize = 0.001
                array_size = int( (maximum - minimum ) / binsize ) + 1
                self.composite[destination2] = { 'minimum': minimum, 'maximum': maximum, 'binsize': binsize }
                self.composite[destination2]['intensities'] = numpy.zeros(array_size,dtype=numpy.float32)
                self.composite[destination2]['n_peaks'] = numpy.zeros(array_size,dtype=numpy.int32)

            elif spectrum_type == 'IT':
                print(f"INFO: Creating a composite spectrum {destination}")
                minimum = 100
                maximum = 400
                binsize = 0.05
                array_size = int( (maximum - minimum ) / binsize ) + 1
                self.composite[destination] = { 'minimum': minimum, 'maximum': maximum, 'binsize': binsize }
                self.composite[destination]['intensities'] = numpy.zeros(array_size,dtype=numpy.float32)
                self.composite[destination]['n_peaks'] = numpy.zeros(array_size,dtype=numpy.int32)

                print(f"INFO: Creating a composite spectrum {destination2}")
                minimum = 0
                maximum = 120
                binsize = 0.2
                array_size = int( (maximum - minimum ) / binsize ) + 1
                self.composite[destination2] = { 'minimum': minimum, 'maximum': maximum, 'binsize': binsize }
                self.composite[destination2]['intensities'] = numpy.zeros(array_size,dtype=numpy.float32)
                self.composite[destination2]['n_peaks'] = numpy.zeros(array_size,dtype=numpy.int32)

            else:
                print("ERROR: Unrecognized spectrum type {spectrum_type}")
                return

        ipeak = -1
        mz_intensity = []
        ipeak_subarray = 0
        for mz in spectrum['m/z array']:
            ipeak += 1
            if mz >= self.composite[destination]['minimum'] and mz <= self.composite[destination]['maximum']:
                bin = int( ( mz - self.composite[destination]['minimum'] ) / self.composite[destination]['binsize'])
                #print(bin,self.regions[region]['n_peaks'],mz)
                self.composite[destination]['n_peaks'][bin] += 1
                self.composite[destination]['intensities'][bin] += spectrum['intensity array'][ipeak]

            intensity = spectrum['intensity array'][ipeak]
            if mz > 132 and abs(mz - precursor_mz) > 5:
                mz_intensity.append( [ipeak_subarray,mz,intensity] )
                ipeak_subarray += 1
            #print("* ",ipeak,mz,intensity)

            # if mz >= (precursor_mz + self.composite[destination2]['minimum']) and mz <= (precursor_mz + self.composite[destination2]['maximum']):
                # bin = int( ( mz - (precursor_mz + self.composite[destination2]['minimum']) ) / self.composite[destination2]['binsize'])
                # #print(bin,self.regions[region]['n_peaks'],mz)
                # self.composite[destination2]['n_peaks'][bin] += 1
                # self.composite[destination2]['intensities'][bin] += spectrum['intensity array'][ipeak]
            # if mz > self.composite[destination2]['maximum'] and mz > (precursor_mz + self.composite[destination2]['maximum']):
                # break

        if len(mz_intensity) < 8:
            return

        sorted_mz_intensity = sorted(mz_intensity, key=lambda pair: pair[2], reverse=True)
        length = 30
        if len(mz_intensity) < 30:
            length = len(mz_intensity)
        clipped_mz_intensity = sorted_mz_intensity[0:(length-1)]
        mz_intensity = sorted(clipped_mz_intensity, key=lambda pair: pair[1], reverse=False)
        i = 0
        while i < length-1:
            mz_intensity[i][0] = i
            i += 1
        sorted_mz_intensity = sorted(mz_intensity, key=lambda pair: pair[2], reverse=True)


        imajor_peak = 0
        while imajor_peak < 20 and imajor_peak < len(mz_intensity):
            ipeak,mz,major_intensity = sorted_mz_intensity[imajor_peak]
            #print(ipeak,mz,major_intensity)
            if major_intensity == 0.0:
                major_intensity = 1.0
            imz = mz
            while ipeak > 0 and mz - imz < self.composite[destination2]['maximum'] and imz > 132 and abs(imz - precursor_mz) > 5:
                bin = int( ( mz - imz ) / self.composite[destination2]['binsize'])
                #print("  - ",ipeak,imz,bin,intensity)
                self.composite[destination2]['n_peaks'][bin] += 1
                self.composite[destination2]['intensities'][bin] += intensity
                ipeak -= 1
                imz = mz_intensity[ipeak][1]
                intensity = mz_intensity[ipeak][2]
            imajor_peak += 1



    ####################################################################################################
    #### Read header
    def read_header(self):

        #### Read the header into a text buffer
        with open(self.mzml_file) as infile:
            buffer = ''
            counter = 0

            #### Read line by line
            for line in infile:

                #### Completely skip the <indexedmzML> tag if present
                match = re.search('<indexedmzML ',line)
                if match:
                    continue

                #### Look for first tag after the header and end when found
                match = re.search('<run ',line)
                if match:
                    break
                if counter > 0:
                    buffer += line
                counter += 1

            #### Finish the XML by closing the tags
            buffer += '  </mzML>\n'

            #### Get the root of the XML and the namespace
            xmlroot = etree.fromstring(buffer)
            namespace = xmlroot.nsmap
            namespace = namespace[None]

            #### List of instruments we know
            instrument_by_category = { 'pureHCD': [ 'MS:1001911|Q Exactive' ],
                'variable': [ 'MS:1001910|LTQ Orbitrap Elite', 'MS:1001742|LTQ Orbitrap Velos', 'MS:1000556|LTQ Orbitrap XL',
                    'MS:1000555|LTQ Orbitrap Discovery' ] }
            instrument_attributes = {}
            for instrument_category in instrument_by_category:
                for instrument_string in instrument_by_category[instrument_category]:
                    accession,name = instrument_string.split('|')
                    instrument_attributes[accession] = { 'category': instrument_category, 'accession': accession, 'name': name }

            #### Get all the CV params in the header and look for ones we know about
            cv_params = xmlroot.findall('.//{http://psi.hupo.org/ms/mzml}cvParam')
            found_instrument = 0
            for cv_param in cv_params:
                accession = cv_param.get('accession')
                if accession in instrument_attributes:
                    #print(f"Found instrument: {instrument_attributes[accession]['name']}")
                    #### Store the attributes about this instrument model
                    model_data = { 'accession': accession, 'name': instrument_attributes[accession]['name'], 'category': instrument_attributes[accession]['category'] }
                    self.metadata['files'][self.mzml_file]['instrument_model'] = model_data
                    found_instrument = 1

            if not found_instrument:
                print("ERROR: Did not recognize the instrument. Please teach me about this instrument.")

        return


    ####################################################################################################
    #### Assess composite spectra
    def assess_composite_spectra(self):
        with open('zzcomposite_spectrum.pickle','rb') as infile:
            spec = pickle.load(infile)

        composite_type = 'lowend_HCD'
        if composite_type not in spec:
            if 'lowend_IT' in spec:
                composite_type = 'lowend_IT'
            else:
                print("ERROR: Unrecognized fragmentation type. Need more training")

        minimum = spec[composite_type]['minimum']
        maximum = spec[composite_type]['maximum']
        binsize = spec[composite_type]['binsize']
        spec[composite_type]['mz'] = numpy.arange(spec[composite_type]['minimum'],spec[composite_type]['maximum'] + spec[composite_type]['binsize'],spec[composite_type]['binsize'])

        ROIs = {
            'TMT6_126': { 'type': 'TMT', 'mz': 126.127725, 'initial_window': 0.01 },
            'TMT6_127': { 'type': 'TMT', 'mz': 127.124760, 'initial_window': 0.01 },
            'TMT6_128': { 'type': 'TMT', 'mz': 128.134433, 'initial_window': 0.01 },
            'TMT6_129': { 'type': 'TMT', 'mz': 129.131468, 'initial_window': 0.01 },
            'TMT6_130': { 'type': 'TMT', 'mz': 130.141141, 'initial_window': 0.01 },
            'TMT6_131': { 'type': 'TMT', 'mz': 131.138176, 'initial_window': 0.01 },
            'TMT6_nterm': { 'type': 'TMT', 'mz': 230.1702, 'initial_window': 0.01 },
            'iTRAQ_114': { 'type': 'iTRAQ', 'mz': 114.11068, 'initial_window': 0.01 },
            'iTRAQ_115': { 'type': 'iTRAQ', 'mz': 115.107715, 'initial_window': 0.01 },
            'iTRAQ_116': { 'type': 'iTRAQ', 'mz': 116.111069, 'initial_window': 0.01 },
            'iTRAQ_117': { 'type': 'iTRAQ', 'mz': 117.114424, 'initial_window': 0.01 },
            'iTRAQ4_nterm': { 'type': 'iTRAQ4', 'mz': 145.109, 'initial_window': 0.01 },
            'iTRAQ8_118': { 'type': 'iTRAQ8', 'mz': 118.111459, 'initial_window': 0.01 },
            'iTRAQ8_119': { 'type': 'iTRAQ8', 'mz': 119.114814, 'initial_window': 0.01 },
            'iTRAQ8_121': { 'type': 'iTRAQ8', 'mz': 121.121524, 'initial_window': 0.01 },
            'iTRAQ8_113': { 'type': 'iTRAQ8', 'mz': 113.107325, 'initial_window': 0.01 },
        }

        #### What should be look at for ion trap data?
        if composite_type == 'lowend_IT':
            ROIs = {
                'TMT6_nterm': { 'type': 'TMT', 'mz': 230.1702, 'initial_window': 1.0 },
                'iTRAQ4_nterm': { 'type': 'iTRAQ4', 'mz': 145.109, 'initial_window': 1.0 },
            }


        for ROI in ROIs:
            #print(f"INFO: Looking for {ROI} at {ROIs[ROI]['mz']}")
            center = ROIs[ROI]['mz']
            range = ROIs[ROI]['initial_window']
            first_bin = int((center - range/2.0 - minimum) / binsize)
            last_bin = int((center + range/2.0 - minimum) / binsize)

            #### Look for the largest peak
            peak = self.find_peak(spec,ROIs,ROI,composite_type)
            ROIs[ROI]['peak'] = peak
            #print(peak)

        self.metadata['files'][self.mzml_file]['ROIs'] = ROIs


    ####################################################################################################
    #### Peak finding routine
    def find_peak(self,spec,ROIs,ROI,composite_type):
        minimum = spec[composite_type]['minimum']
        maximum = spec[composite_type]['maximum']
        binsize = spec[composite_type]['binsize']
        center = ROIs[ROI]['mz']
        range = ROIs[ROI]['initial_window']
        first_bin = int((center - range/2.0 - minimum) / binsize)
        last_bin = int((center + range/2.0 - minimum) / binsize)

        peak = { 'mode_bin': { 'ibin': 0, 'mz': 0, 'intensity': 0, 'n_spectra': 0 } }

        #### Find the tallest peak in the window
        ibin = first_bin
        while ibin < last_bin:
            if spec[composite_type]['intensities'][ibin] > peak['mode_bin']['intensity']:
                peak['mode_bin']['ibin'] = ibin
                peak['mode_bin']['mz'] = float(spec[composite_type]['mz'][ibin])
                peak['mode_bin']['intensity'] = float(spec[composite_type]['intensities'][ibin])
                peak['mode_bin']['n_spectra'] = int(spec[composite_type]['n_peaks'][ibin])
            ibin += 1

        #### If there were 0 peaks found, return a verdict now
        if peak['mode_bin']['ibin'] == 0:
            peak['assessment'] = { 'is_found': False, 'fraction': 0.0 }
            return(peak)

        #### Grow the extent to get most of the signal
        ibin = peak['mode_bin']['ibin']
        iextent = 1
        done = 0
        prev_n_spectra = peak['mode_bin']['n_spectra']
        peak['extended'] = { 'extent': 0, 'intensity': float(peak['mode_bin']['intensity']), 'n_spectra': int(peak['mode_bin']['n_spectra']) }
        #print(spec[composite_type]['n_peaks'][ibin-10:ibin+10])
        while not done:
            prev_intensity = peak['extended']['intensity']
            prev_n_spectra = peak['extended']['n_spectra']
            new_intensity = prev_intensity + spec[composite_type]['intensities'][ibin-iextent] + spec[composite_type]['intensities'][ibin+iextent]
            new_n_spectra = prev_n_spectra + spec[composite_type]['n_peaks'][ibin-iextent] + spec[composite_type]['n_peaks'][ibin+iextent]
            #print(new_intensity/prev_intensity, new_n_spectra/prev_n_spectra)
            peak['extended'] = { 'extent': iextent, 'intensity': float(new_intensity), 'n_spectra': int(new_n_spectra) }
            if new_n_spectra/prev_n_spectra < 1.02:
                #print('Extent reached by n_spectra')
                done = 1
            if iextent > 15:
                #print('Extent reached by max extent')
                done = 1
            iextent += 1

        if peak['extended']['n_spectra'] > 100:
            extent = peak['extended']['extent'] * 2
            x = spec[composite_type]['mz'][ibin-extent:ibin+extent]
            y = spec[composite_type]['intensities'][ibin-extent:ibin+extent]

            n = len(x)
            center = int(n/2)
            binsize = x[center]-x[center-1]
            popt,pcov = curve_fit(gaussian_function,x,y,p0=[y[center],x[center],binsize])

            peak['fit'] = { 'mz': popt[1], 'intensity': popt[0], 'sigma': popt[2], 'delta_mz': popt[1]-ROIs[ROI]['mz'], 'delta_ppm': (popt[1]-ROIs[ROI]['mz'])/ROIs[ROI]['mz']*1e6 }
            peak['assessment'] = { 'is_found': True, 'fraction': 0.0, 'comment': 'Peak found and fit' }

            if 0:
                plt.step(x+binsize/2.0,y,'b+:')
                plt.plot(x,gaussian_function(x,*popt),'ro:')
                title = "Fit results: mz = %.4f, sigma = %.4f" % (popt[1], popt[2])
                plt.title(title)
                plt.show()

        else:
            peak['assessment'] = { 'is_found': False, 'fraction': 0.0, 'comment': 'Too few peaks found in ROI' }

        return(peak)


    ####################################################################################################
    #### Assess Regions of Interest to make calls
    def assess_ROIs(self):
        results = { 'labeling': { 'scores': { 'TMT': 0, 'TMT6': 0, 'TMT10': 0, 'iTRAQ': 0, 'iTRAQ4': 0, 'iTRAQ8': 0 } } }

        #### Determine what the denominator of MS2 spectra should be
        n_ms2_spectra = self.metadata['files'][self.mzml_file]['spectra_stats']['n_ms2_spectra']
        if ( 'n_HCD_spectra' in self.metadata['files'][self.mzml_file]['spectra_stats'] and
                self.metadata['files'][self.mzml_file]['spectra_stats']['n_HCD_spectra'] > 0 and 
                self.metadata['files'][self.mzml_file]['spectra_stats']['n_HCD_spectra'] <
                self.metadata['files'][self.mzml_file]['spectra_stats']['n_ms2_spectra'] ):
            n_ms2_spectra = self.metadata['files'][self.mzml_file]['spectra_stats']['n_HCD_spectra']
 
        for ROI in self.metadata['files'][self.mzml_file]['ROIs']:
            ROIobj = self.metadata['files'][self.mzml_file]['ROIs'][ROI]
            if ROIobj['peak']['assessment']['is_found']:
                type = ROIobj['type']
                n_spectra = ROIobj['peak']['extended']['n_spectra']
                ROIobj['peak']['assessment']['fraction'] = n_spectra / n_ms2_spectra
                results['labeling']['scores'][type] += ROIobj['peak']['assessment']['fraction']

        #### Make the call for TMT or iTRAQ
        if results['labeling']['scores']['TMT'] > results['labeling']['scores']['iTRAQ']:
            if results['labeling']['scores']['TMT'] > 2:
                results['labeling']['call'] = 'TMT'
            elif results['labeling']['scores']['TMT'] < 1:
                results['labeling']['call'] = 'none'
            else:
                results['labeling']['call'] = 'ambiguous'
        else:
            if results['labeling']['scores']['iTRAQ'] > 1:
                if results['labeling']['scores']['iTRAQ4'] > results['labeling']['scores']['iTRAQ8']:
                    results['labeling']['call'] = 'iTRAQ4'
                else:
                    results['labeling']['call'] = 'iTRAQ8'
            elif results['labeling']['scores']['iTRAQ'] < 0.5:
                results['labeling']['call'] = 'none'
            else:
                results['labeling']['call'] = 'ambiguous'

        results['fragmentation_type'] = self.metadata['files'][self.mzml_file]['spectra_stats']['fragmentation_type']
        self.metadata['files'][self.mzml_file]['summary'] = results


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
#### Gaussian function used for curve fitting
def gaussian_function(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Creates an index for an MSP spectral library file')
    argparser.add_argument('--use_cache', action='store', help='Set to true to use the cached file instead of regenerating')
    argparser.add_argument('--refresh', action='count', default=0, help='If set, existing metadata for a file will be overwritten rather than skipping the file')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('files', type=str, nargs='+', help='Filenames of one or more mzML files to read')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None: verbose = 1

    #### Loop over all the files to ensure that they are really there before starting work
    for file in params.files:
        if not os.path.isfile(file):
            print(f"ERROR: File '{file}' not found or not a file")
            return

    #### Load the current metadata file in order to update it
    study = MetadataHandler(verbose=verbose)
    study.read_or_create()

    #### Loop over all the files processing them
    for file in params.files:

        #### If this file is already in the metadata store, purge it and generate new
        if file in study.metadata['files']:
            if params.refresh is not None and params.refresh > 0:
                study.metadata['files'][file] = None
            else:
                if verbose >= 1: eprint(f"INFO: Already have results for '{file}'. Skipping..")
                continue

        #### Assess the mzML file
        assessor = MzMLAssessor(file, metadata=study.metadata, verbose=verbose)
        assessor.read_header()
        if not params.use_cache:
            assessor.read_spectra()
        assessor.assess_composite_spectra()
        assessor.assess_ROIs()

    #### Infer parameters based on the latest data
    study.infer_search_criteria()

    #### Write out our state of mind
    study.store()


#### For command line usage
if __name__ == "__main__": main()
