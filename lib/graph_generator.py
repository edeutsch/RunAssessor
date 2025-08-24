import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from math import sqrt, pi
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader, PdfWriter
import sys
import os
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

class GraphGenerator:
    def __init__(self, metadata_filepath=None, verbose=None):
        default_metadata_filepath = 'study_metadata.json'
        if metadata_filepath == None:
            self.metadata_file = default_metadata_filepath
        elif not metadata_filepath.endswith('.json'):
            self.metadata_file = metadata_filepath + ".json"
        else:
            self.metadata_file = metadata_filepath

         #### Load JSON data
        try:
            with open(self.metadata_file, "r") as f:
                self.data = json.load(f)
        except:
            eprint("INFO: Unable to read JSON file")
            return
        
        if verbose != None and verbose != 0:
            self.verbose = 1
        else:
            self.verbose = 0
        
        self.files = {}

    def fmt(self, v):
        return f"{v:.3f}" if isinstance(v, (int, float)) else str(v)

    def buildGraphs(self):

        self.files = self.data.get("files", {})

        coverpage = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Delta_graphs_documentation.pdf")

        #### Output PDF file
        output_pdf_path = self.metadata_file.replace(".json", ".histograms_with_chi2.pdf")
        with PdfPages(output_pdf_path) as pdf:

            for filename, file_data in self.files.items():
                if file_data.get("spectra_stats", {}).get("acquisition_type", 'none') == "DIA":
                    continue
                pre_tol = {}
                hist_ppm = {}
                amplitude_ppm = None
                mu_ppm = None
                sigma_ppm = None
                y_offset_ppm = None
                chi_squared_ppm = None

                hist_time = {}
                time_initial_level = None
                pulse_start = None
                peak_level = None
                decay_constant = None
                final_level = None
                pulse_duration = None
                chi_squared_time = None

                time_bin_centers = None
                time_counts = None
                ppm_bin_centers = None
                ppm_counts = None
                precursor_stats = file_data.get("summary", {}).get('precursor stats', {})
                dynamic_exclusion_window = precursor_stats.get('dynamic exclusion window', {})

                try:
                    pre_tol = precursor_stats.get('precursor tolerance', {})
                    hist_ppm = pre_tol.get("histogram_ppm", {})
                    amplitude_ppm = pre_tol.get('fit_ppm', {}).get("intensity", [])
                    mu_ppm = pre_tol.get('fit_ppm', {}).get("delta_ppm peak", [])
                    sigma_ppm = pre_tol.get('fit_ppm', {}).get("sigma (ppm)", [])
                    y_offset_ppm = pre_tol.get('fit_ppm', {}).get("y_offset", [])
                    
                    chi_squared_ppm = pre_tol.get("fit_ppm", {}).get("chi_squared", "N/A")
                except:
                    pass

                try:
                    hist_time = dynamic_exclusion_window.get("histogram_time", {})
                    time_initial_level = dynamic_exclusion_window.get("fit_pulse_time", {}).get("inital level", [])
                    pulse_start = dynamic_exclusion_window.get("fit_pulse_time", {}).get("pulse start", [])
                    peak_level = dynamic_exclusion_window.get("fit_pulse_time", {}).get("peak level", [])
                    decay_constant = dynamic_exclusion_window.get("fit_pulse_time", {}).get("decay constant", [])
                    final_level = dynamic_exclusion_window.get("fit_pulse_time", {}).get("final level", [])
                    pulse_duration = dynamic_exclusion_window.get("fit_pulse_time", {}).get("pulse duration", [])
                    chi_squared_time = dynamic_exclusion_window.get("fit_pulse_time", {}).get("chi_squared", "N/A")
                except:
                    pass

                

                try:
                    time_bin_centers = np.array(hist_time.get("bin_centers", []))
                    time_counts = np.array(hist_time.get("counts", []))
                except:
                    pass
                
                try:
                    ppm_bin_centers = np.array(hist_ppm.get("bin_centers", []))
                    ppm_counts = np.array(hist_ppm.get("counts", []))
                except:
                    pass


                fig, axs = plt.subplots(2, 2, figsize=(12, 8), gridspec_kw={'height_ratios': [4, 1]})
                fig.suptitle(filename, fontsize=14)

                #### Time histogram
                ax_time = axs[0, 0]
                try:
                    width_time = (time_bin_centers[1] - time_bin_centers[0]) if len(time_bin_centers) > 1 else 1
                    ax_time.bar(time_bin_centers, time_counts, width=width_time, color='skyblue', edgecolor='black', label='Time Data')
                
                    try:
                        if peak_level <= 0:
                            pass
                        else:
                            t_fine = np.linspace(min(time_bin_centers), max(time_bin_centers), 500)
                            y = np.full_like(t_fine, time_initial_level, dtype=float)

                            pulse_end = pulse_start + pulse_duration
                            pulse_mask = (t_fine >= pulse_start) & (t_fine < pulse_end)
                            y[pulse_mask] = peak_level
                            decay_mask = t_fine >= pulse_end
                            safe_tau = max(decay_constant, 1e-8)
                            exponent = -(t_fine[decay_mask] - pulse_end) / safe_tau
                            exponent = np.clip(exponent, -700, 700)
                            y[decay_mask] = final_level + (peak_level - final_level) * np.exp(exponent)

                            ax_time.plot(t_fine, y, color='blue', linewidth=2, label='Time Decay Fit')
                    except:
                        pass
                except:
                    ax_time.bar([], [])

                ax_time.set_title("ΔTime Between Precursors")
                ax_time.set_xlabel("ΔTime (seconds)")
                ax_time.set_ylabel("Counts")
                ax_time.legend()

                #### PPM histogram
                ax_ppm = axs[0, 1]
                try: 
                    width_ppm = (ppm_bin_centers[1] - ppm_bin_centers[0]) if len(ppm_bin_centers) > 1 else 1
                    ax_ppm.bar(ppm_bin_centers, ppm_counts, width=width_ppm, color='salmon', edgecolor='black', label='PPM Data')

                
                    try:
                        
                        if amplitude_ppm <= 0:
                            pass
                        else:
                            x_ppm = np.linspace(min(ppm_bin_centers), max(ppm_bin_centers), 500)
                            scale = amplitude_ppm * sigma_ppm * sqrt(2 * pi)
                            y_ppm = scale * norm.pdf(x_ppm, mu_ppm, sigma_ppm) + y_offset_ppm
                            ax_ppm.plot(x_ppm, y_ppm, color='darkred', linewidth=2, label='PPM Gaussian Fit')
                    except:
                        pass
                except:
                    ax_ppm.bar([], [])
                    
                ax_ppm.set_title("Precursor m/z Variance in PPM")
                ax_ppm.set_xlabel("PPM Variance")
                ax_ppm.set_ylabel("Counts")
                ax_ppm.legend()

                #### Tables
                axs[1, 0].axis('off')
                axs[1, 1].axis('off')


                ## Table under Time histogram
                try:
                    if peak_level == None or peak_level == 0:
                        time_table_data = [["No fit found", ""]]
                    else:
                        time_table_data = [
                            ["Chi² Time", self.fmt(chi_squared_time)],
                            ["Initial Level", self.fmt(time_initial_level)],
                            ["Pulse Start", self.fmt(pulse_start)],
                            ['Pulse Duration', self.fmt(pulse_duration)],
                            ["Peak Level",self.fmt(peak_level)],
                            ["Decay Constant",self.fmt(decay_constant)],
                            ["Final Level", self.fmt(final_level)],
                        ]
                except:
                    time_table_data = [["No fit found", ""]]
                    
                table_time = axs[1, 0].table(
                        cellText=time_table_data,
                        colLabels=["Time Fit", "Value"],
                        loc='center',
                        cellLoc='left',
                        colWidths=[0.35, 0.3],
                        bbox=[0, -0.25, 1, 1] 
                    )
                table_time.scale(1.2, 3)

                ## Table under PPM histogram
                try:
                    if amplitude_ppm == None or amplitude_ppm  == 0:
                        ppm_table_data = [["No fit found", ""]]
                    else:
                        ppm_table_data = [
                            ["Chi² PPM", self.fmt(chi_squared_ppm)],
                            ["Amplitude", self.fmt(amplitude_ppm)],
                            ["Mu (Peak Center)", self.fmt(mu_ppm)],
                            ["Sigma", self.fmt(sigma_ppm)],
                            ["Y Offset", self.fmt(y_offset_ppm)],
                        ]
                except:
                    ppm_table_data = [["No fit found", ""]]

                table_ppm = axs[1, 1].table(
                    cellText=ppm_table_data,
                    colLabels=["PPM Fit", "Value"],
                    loc='center',
                    cellLoc='left',
                    colWidths=[0.35, 0.3],
                    bbox=[0, -0.25, 1, 1]  
                )
                table_ppm.scale(1.2, 1.6)


                pdf.savefig(fig)
                plt.close(fig)

            if self.verbose >= 1:
                print(f"Saved plot for {filename}")

        #Add explanatory Cover Page
        writer = PdfWriter()
        cover_reader = PdfReader(coverpage)
        for page in cover_reader.pages:
            writer.add_page(page)

        try:
            output_reader = PdfReader(output_pdf_path)
            for page in output_reader.pages:
                writer.add_page(page)

            # Write the final combined PDF
            with open(output_pdf_path, "wb") as f:
                writer.write(f)
        except:
            eprint("WARNING: No doumentation page added to chi_sqared graphs")

        if self.verbose >= 1:
            eprint("Saving PDF to:", os.path.abspath(output_pdf_path))
        return self.metadata_file


    
    def plot_precursor_loss_composite_spectra(self, assessor, pdf):
        destinations_list = list(assessor.composite.keys())
        file_name_root = assessor.mzml_file.split('.')[0]
        for destination in destinations_list:
            if destination.startswith("precursor_loss_"):

                intensities = assessor.composite[destination]['intensities']
                maximum = assessor.composite[destination]['maximum']
                minimum = assessor.composite[destination]['minimum']
                binsize = assessor.composite[destination]['binsize']
                mz = np.arange((maximum - minimum) / binsize + 1) * binsize + minimum

                water_z2 = 9.00528235
                phosphoric_acid_z2 = 48.98844785
                try:
                    water_z2_extent_bins = assessor.metadata['files'][assessor.mzml_file]['neutral_loss_peaks'][destination]['water_z2']['peak']['extended']['extent']
                except:
                    water_z2_extent_bins = 0
                try:
                    phosphoric_acid_z2_extent_bins = assessor.metadata['files'][assessor.mzml_file]['neutral_loss_peaks'][destination]['phosphoric_acid_z2']['peak']['extended']['extent']
                except:
                    phosphoric_acid_z2_extent_bins = 0

                water_z2_axis_min = water_z2 - 30 * binsize
                water_z2_axis_max = water_z2 + 30 * binsize
                water_z2_extent_min = water_z2 - water_z2_extent_bins * binsize
                water_z2_extent_max = water_z2 + water_z2_extent_bins * binsize
                water_z2_fitting_min = water_z2 - 2 * water_z2_extent_bins * binsize
                water_z2_fitting_max = water_z2 + 2 * water_z2_extent_bins * binsize

                phosphoric_acid_z2_axis_min = phosphoric_acid_z2 - 30 * binsize
                phosphoric_acid_z2_axis_max = phosphoric_acid_z2 + 30 * binsize
                phosphoric_acid_z2_extent_min = phosphoric_acid_z2 - phosphoric_acid_z2_extent_bins * binsize
                phosphoric_acid_z2_extent_max = phosphoric_acid_z2 + phosphoric_acid_z2_extent_bins * binsize
                phosphoric_acid_z2_fitting_min = phosphoric_acid_z2 - 2 * phosphoric_acid_z2_extent_bins * binsize
                phosphoric_acid_z2_fitting_max = phosphoric_acid_z2 + 2 * phosphoric_acid_z2_extent_bins * binsize

                sum_type = destination.replace("precursor_loss_", "")
                summary = self.files[assessor.mzml_file]['summary'][sum_type]

                # Colors based on detection
                water_color = 'palegreen' if summary.get('has water_loss') else 'navajowhite'
                phospho_color = 'palegreen' if summary.get('has phospho_spectra') else 'navajowhite'

                correct_bar_water = 'limegreen' if summary.get('has water_loss') else 'orange'
                correct_bar_phospho = 'limegreen' if summary.get('has phospho_spectra') else 'orange'

                ratio = summary.get('z=2 phosphoric_acid to z=2 water_loss intensity ratio', 'N/A')
                num_water = assessor.metadata['files'][assessor.mzml_file]['neutral_loss_peaks'][destination]['water_z2']['peak']['mode_bin']['n_spectra']
                num_phospho = assessor.metadata['files'][assessor.mzml_file]['neutral_loss_peaks'][destination]['phosphoric_acid_z2']['peak']['mode_bin']['n_spectra']
                abs_diff_delta_mz = summary.get('absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss', 'N/A')

                # Plot water z=2
                if water_z2_extent_bins * 2 < 30:
                    plt.axvspan(xmin=water_z2_fitting_min, xmax=water_z2_fitting_max, color='lightsteelblue', alpha=0.25, lw=0)
                plt.axvspan(xmin=water_z2_extent_min, xmax=water_z2_extent_max, color=water_color, alpha=0.75, lw=0)
                plt.axvline(x=water_z2, color=correct_bar_water, alpha=0.75)
                water_z2_range = (mz >= water_z2_axis_min) & (mz <= water_z2_axis_max)
                plt.plot(mz[water_z2_range], intensities[water_z2_range])
                
                ## Plot possible fit
                little_label = ''
                try:
                    if self.files[assessor.mzml_file]['neutral_loss_peaks'][destination]['water_z2']['peak']['assessment']['is_found']:
                        fit_param = self.files[assessor.mzml_file]['neutral_loss_peaks'][destination]['water_z2']['peak']['fit']
                        sigma_mz = fit_param['sigma_mz']
                        mu_mz = fit_param['mz']
                        y_offset = fit_param['y_offset']
                        mz_subset = mz[water_z2_range]
                        intensity_subset = intensities[water_z2_range]

                        gaussian = norm.pdf(mz_subset, mu_mz, sigma_mz)
                        scaled_gaussian = gaussian * np.max(intensity_subset) / np.max(gaussian) + y_offset/ np.max(gaussian)

                        plt.plot(mz_subset, scaled_gaussian, color='darkred', linewidth=2)
                    else:
                        if self.verbose >= 1:
                            eprint(f"No water loss z =2 peak found in {assessor.mzml_file}")
                        little_label = " (no fit found)"
                except:
                    pass
                plt.xlabel("m/z loss (water z=2)" + little_label + "\nMode of water z=2 spectra: " + str(num_water))
                plt.ylabel("Intensity")
                try:
                    plt.title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz:.4f}", fontsize=8.5)
                except:
                    plt.title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz}", fontsize=8.5)
                plt.tight_layout()
                pdf.savefig()
                plt.close()

                # Plot phosphoric acid z=2
                if phosphoric_acid_z2_extent_bins * 2 < 30:
                    plt.axvspan(xmin=phosphoric_acid_z2_fitting_min, xmax=phosphoric_acid_z2_fitting_max, color='lightsteelblue', alpha=0.25, lw=0)
                plt.axvspan(xmin=phosphoric_acid_z2_extent_min, xmax=phosphoric_acid_z2_extent_max, color=phospho_color, alpha=0.75, lw=0)
                plt.axvline(x=phosphoric_acid_z2, color=correct_bar_phospho, alpha=0.75)
                phosphoric_acid_z2_range = (mz >= phosphoric_acid_z2_axis_min) & (mz <= phosphoric_acid_z2_axis_max)
                plt.plot(mz[phosphoric_acid_z2_range], intensities[phosphoric_acid_z2_range])
                plt.xlabel("m/z loss (phosphoric acid z=2)")
                little_label = ""
                try:
                    if self.files[assessor.mzml_file]['neutral_loss_peaks'][destination]['phosphoric_acid_z2']['peak']['assessment']['is_found']:
                        fit_param = self.files[assessor.mzml_file]['neutral_loss_peaks'][destination]['phosphoric_acid_z2']['peak']['fit']
                        sigma_mz = fit_param['sigma_mz']
                        mu_mz = fit_param['mz']
                        y_offset = fit_param['y_offset']
                        mz_subset = mz[phosphoric_acid_z2_range]
                        intensity_subset = intensities[phosphoric_acid_z2_range]

                        gaussian = norm.pdf(mz_subset, mu_mz, sigma_mz)
                        scaled_gaussian = gaussian * np.max(intensity_subset) / np.max(gaussian) + y_offset

                        plt.plot(mz_subset, scaled_gaussian, color='darkred', linewidth=2)
                    else:
                        if self.verbose >= 1:
                            eprint(f"No phosphoric acid z =2 peak found in {assessor.mzml_file}")
                        little_label = " (no fit found)"
                except:
                    pass

                plt.xlabel("m/z loss (phosphoric acid z=2)" + little_label + "\nMode of phosphoric acid z=2 spectra: " + str(num_phospho))
                plt.ylabel("Intensity")
                try:
                    plt.title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz:.4f}\nPhosphoric acid to water intensity ratio: {ratio:.2f}", fontsize=8.5)
                except:
                    try:
                        plt.title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz}\nPhosphoric acid to water intensity ratio: {ratio:.2f}", fontsize=8.5)
                    except:
                        plt.title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz}\nPhosphoric acid to water intensity ratio: {ratio}", fontsize=8.5)
                plt.tight_layout()
                pdf.savefig()
                plt.close()

                # Full neutral loss spectrum
                plt.plot(mz, intensities)
                plt.xlabel("m/z loss")
                plt.ylabel("Intensity")
                plt.title(f"{file_name_root} ({destination})\nFull Neutral Loss Spectrum", fontsize=8.5)
                plt.tight_layout()
                pdf.savefig()
                plt.close()
        nl_pdf = self.metadata_file.replace(".json", ".NLplots.pdf")
        if self.verbose >= 1:
            eprint(f"All {assessor.mzml_file} neutral loss spectra saved to: {nl_pdf}")
    
    # Manuscript figure graphing: neutral loss example and low-end example (ADPhosTMT3R1.mzML.gz from PXD010807)
    def plot_NL_LE_figures(self, assessor):
        destinations_list = list(assessor.composite.keys())
        file_name_root = assessor.mzml_file.split('.')[0]
        for destination in destinations_list:

            intensities = assessor.composite[destination]['intensities']
            maximum = assessor.composite[destination]['maximum']
            minimum = assessor.composite[destination]['minimum']
            binsize = assessor.composite[destination]['binsize']
            mz = np.arange((maximum - minimum) / binsize + 1) * binsize + minimum

            if destination.startswith("precursor_loss_"):

                water_z2 = 9.00528235
                phosphoric_acid_z2 = 48.98844785

                try:
                    water_z2_extent_bins = assessor.metadata['files'][assessor.mzml_file]['neutral_loss_peaks'][destination]['water_z2']['peak']['extended']['extent']
                    phosphoric_acid_z2_extent_bins = assessor.metadata['files'][assessor.mzml_file]['neutral_loss_peaks'][destination]['phosphoric_acid_z2']['peak']['extended']['extent']
                except:
                    water_z2_extent_bins = 0
                    phosphoric_acid_z2_extent_bins = 0

                water_z2_axis_min = water_z2 - 30 * binsize
                water_z2_axis_max = water_z2 + 30 * binsize
                water_z2_extent_min = water_z2 - water_z2_extent_bins * binsize
                water_z2_extent_max = water_z2 + water_z2_extent_bins * binsize
                water_z2_fitting_min = water_z2 - 2 * water_z2_extent_bins * binsize
                water_z2_fitting_max = water_z2 + 2 * water_z2_extent_bins * binsize

                phosphoric_acid_z2_axis_min = phosphoric_acid_z2 - 30 * binsize
                phosphoric_acid_z2_axis_max = phosphoric_acid_z2 + 30 * binsize
                phosphoric_acid_z2_extent_min = phosphoric_acid_z2 - phosphoric_acid_z2_extent_bins * binsize
                phosphoric_acid_z2_extent_max = phosphoric_acid_z2 + phosphoric_acid_z2_extent_bins * binsize
                phosphoric_acid_z2_fitting_min = phosphoric_acid_z2 - 2 * phosphoric_acid_z2_extent_bins * binsize
                phosphoric_acid_z2_fitting_max = phosphoric_acid_z2 + 2 * phosphoric_acid_z2_extent_bins * binsize

                sum_type = destination.replace("precursor_loss_", "")
                summary = self.files[assessor.mzml_file]['summary'][sum_type]

                # Colors based on detection
                water_color = 'palegreen' if summary.get('has water_loss') else 'navajowhite'
                phospho_color = 'palegreen' if summary.get('has phospho_spectra') else 'navajowhite'

                correct_bar_water = 'limegreen' if summary.get('has water_loss') else 'orange'
                correct_bar_phospho = 'limegreen' if summary.get('has phospho_spectra') else 'orange'

                ratio = summary.get('z=2 phosphoric_acid to z=2 water_loss intensity ratio', 'N/A')
                num_water = assessor.metadata['files'][assessor.mzml_file]['neutral_loss_peaks'][destination]['water_z2']['peak']['mode_bin']['n_spectra']
                num_phospho = assessor.metadata['files'][assessor.mzml_file]['neutral_loss_peaks'][destination]['phosphoric_acid_z2']['peak']['mode_bin']['n_spectra']
                abs_diff_delta_mz = summary.get('absolute difference in delta m/z for z=2 phosphoric_acid_loss and z=2 water_loss', 'N/A')

                fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16, 6), layout='constrained')

                # Plot water z=2
                if water_z2_extent_bins * 2 < 30:
                    ax[0].axvspan(xmin=water_z2_fitting_min, xmax=water_z2_fitting_max, color='lightsteelblue', alpha=0.25, lw=0)
                ax[0].axvspan(xmin=water_z2_extent_min, xmax=water_z2_extent_max, color=water_color, alpha=0.75, lw=0)
                ax[0].axvline(x=water_z2, color=correct_bar_water, alpha=0.75)
                water_z2_range = (mz >= water_z2_axis_min) & (mz <= water_z2_axis_max)
                ax[0].plot(mz[water_z2_range], intensities[water_z2_range])

                ## Plot possible fit
                little_label = ''
                try:
                    if self.files[assessor.mzml_file]['neutral_loss_peaks'][destination]['water_z2']['peak']['assessment']['is_found']:
                        fit_param = self.files[assessor.mzml_file]['neutral_loss_peaks'][destination]['water_z2']['peak']['fit']
                        sigma_mz = fit_param['sigma_mz']
                        mu_mz = fit_param['mz']
                        y_offset = fit_param['y_offset']
                        mz_subset = mz[water_z2_range]
                        intensity_subset = intensities[water_z2_range]

                        gaussian = norm.pdf(mz_subset, mu_mz, sigma_mz)
                        scaled_gaussian = gaussian * np.max(intensity_subset) / np.max(gaussian) + y_offset/ np.max(gaussian)

                        ax[0].plot(mz_subset, scaled_gaussian, color='darkred', linewidth=2)
                    else:
                        if self.verbose >= 1:
                            eprint(f"No water loss z =2 peak found in {assessor.mzml_file}")
                        little_label = " (no fit found)"
                except:
                    pass

                ax[0].set_xlabel("m/z loss (water z=2)" + little_label + "\nMode of water z=2 spectra: " + str(num_water))
                ax[0].set_ylabel("Intensity")
                try:
                    ax[0].set_title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz:.4f}", fontsize=12)
                except:
                    ax[0].set_title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz}", fontsize=12)

                # Plot phosphoric acid z=2
                if phosphoric_acid_z2_extent_bins * 2 < 30:
                    ax[1].axvspan(xmin=phosphoric_acid_z2_fitting_min, xmax=phosphoric_acid_z2_fitting_max, color='lightsteelblue', alpha=0.25, lw=0)
                ax[1].axvspan(xmin=phosphoric_acid_z2_extent_min, xmax=phosphoric_acid_z2_extent_max, color=phospho_color, alpha=0.75, lw=0)
                ax[1].axvline(x=phosphoric_acid_z2, color=correct_bar_phospho, alpha=0.75)
                phosphoric_acid_z2_range = (mz >= phosphoric_acid_z2_axis_min) & (mz <= phosphoric_acid_z2_axis_max)
                ax[1].plot(mz[phosphoric_acid_z2_range], intensities[phosphoric_acid_z2_range])

                ## Plot possible fit
                little_label = ""
                try:
                    if self.files[assessor.mzml_file]['neutral_loss_peaks'][destination]['phosphoric_acid_z2']['peak']['assessment']['is_found']:
                        fit_param = self.files[assessor.mzml_file]['neutral_loss_peaks'][destination]['phosphoric_acid_z2']['peak']['fit']
                        sigma_mz = fit_param['sigma_mz']
                        mu_mz = fit_param['mz']
                        y_offset = fit_param['y_offset']
                        mz_subset = mz[phosphoric_acid_z2_range]
                        intensity_subset = intensities[phosphoric_acid_z2_range]

                        gaussian = norm.pdf(mz_subset, mu_mz, sigma_mz)
                        scaled_gaussian = gaussian * np.max(intensity_subset) / np.max(gaussian) + y_offset

                        ax[1].plot(mz_subset, scaled_gaussian, color='darkred', linewidth=2)
                    else:
                        if self.verbose >= 1:
                            eprint(f"No phosphoric acid z =2 peak found in {assessor.mzml_file}")
                        little_label = " (no fit found)"
                except:
                    pass

                ax[1].set_xlabel("m/z loss (phosphoric acid z=2)" + little_label + "\nMode of phosphoric acid z=2 spectra: " + str(num_phospho))
                ax[1].set_ylabel("Intensity")
                try:
                    ax[1].set_title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz:.4f}\nPhosphoric acid to water intensity ratio: {ratio:.2f}", fontsize=12)
                except:
                    try:
                        ax[1].set_title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz}\nPhosphoric acid to water intensity ratio: {ratio:.2f}", fontsize=12)
                    except:
                        ax[1].set_title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz}\nPhosphoric acid to water intensity ratio: {ratio}", fontsize=12)

                ax[0].text(-0.07, 1.1, "a", transform=ax[0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')
                ax[1].text(-0.06, 1.1, "b", transform=ax[1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

                plt.tight_layout()

                # Save as SVG
                plt.savefig(f"{file_name_root}_NL_multipanel.svg", format="svg", bbox_inches='tight')

            
            if destination.startswith("lowend_"):

                TMT10_129_C = 129.13779
                TMT6plex = 230.1702
                IH = 110.07127
                y_K = 147.1128
                
                try:
                    TMT10_129_C_extent_bins = assessor.metadata['files'][assessor.mzml_file]['lowend_peaks'][destination]['TMT10_129_C']['peak']['extended']['extent']
                    TMT6plex_extent_bins = assessor.metadata['files'][assessor.mzml_file]['lowend_peaks'][destination]['TMT6plex']['peak']['extended']['extent']
                    IH_extent_bins = assessor.metadata['files'][assessor.mzml_file]['lowend_peaks'][destination]['IH']['peak']['extended']['extent']
                    y_K_extent_bins = assessor.metadata['files'][assessor.mzml_file]['lowend_peaks'][destination]['y{K}']['peak']['extended']['extent']
                except:
                    TMT10_129_C_extent_bins = 0
                    TMT6plex_extent_bins = 0
                    IH_extent_bins = 0
                    y_K_extent_bins = 0

                TMT10_129_C_axis_min = TMT10_129_C - 10 * binsize
                TMT10_129_C_axis_max = TMT10_129_C + 10 * binsize
                TMT10_129_C_extent_min = TMT10_129_C - TMT10_129_C_extent_bins * binsize
                TMT10_129_C_extent_max = TMT10_129_C + TMT10_129_C_extent_bins * binsize
                TMT10_129_C_fitting_min = TMT10_129_C - 2 * TMT10_129_C_extent_bins * binsize
                TMT10_129_C_fitting_max = TMT10_129_C + 2 * TMT10_129_C_extent_bins * binsize

                TMT6plex_axis_min = TMT6plex - 10 * binsize
                TMT6plex_axis_max = TMT6plex + 10 * binsize
                TMT6plex_extent_min = TMT6plex - TMT6plex_extent_bins * binsize
                TMT6plex_extent_max = TMT6plex + TMT6plex_extent_bins * binsize
                TMT6plex_fitting_min = TMT6plex - 2 * TMT6plex_extent_bins * binsize
                TMT6plex_fitting_max = TMT6plex + 2 * TMT6plex_extent_bins * binsize

                IH_axis_min = IH - 10 * binsize
                IH_axis_max = IH + 10 * binsize
                IH_extent_min = IH - IH_extent_bins * binsize
                IH_extent_max = IH + IH_extent_bins * binsize
                IH_fitting_min = IH - 2 * IH_extent_bins * binsize
                IH_fitting_max = IH + 2 * IH_extent_bins * binsize

                y_K_axis_min = y_K - 10 * binsize
                y_K_axis_max = y_K + 10 * binsize
                y_K_extent_min = y_K - y_K_extent_bins * binsize
                y_K_extent_max = y_K + y_K_extent_bins * binsize
                y_K_fitting_min = y_K - 2 * y_K_extent_bins * binsize
                y_K_fitting_max = y_K + 2 * y_K_extent_bins * binsize

                sum_type = destination.replace("lowend_", "")
                summary = self.files[assessor.mzml_file]['summary'][sum_type]

                # Colors based on detection
                num_TMT10_129_C = assessor.metadata['files'][assessor.mzml_file]['lowend_peaks'][destination]['TMT10_129_C']['peak']['mode_bin']['n_spectra']
                num_TMT6plex = assessor.metadata['files'][assessor.mzml_file]['lowend_peaks'][destination]['TMT6plex']['peak']['mode_bin']['n_spectra']
                num_IH = assessor.metadata['files'][assessor.mzml_file]['lowend_peaks'][destination]['IH']['peak']['mode_bin']['n_spectra']
                num_y_K = assessor.metadata['files'][assessor.mzml_file]['lowend_peaks'][destination]['y{K}']['peak']['mode_bin']['n_spectra']

                TMT10_129_C_color = 'palegreen' if self.files[assessor.mzml_file]['lowend_peaks'][destination]['TMT10_129_C']['peak']['assessment']['is_found'] else 'navajowhite'
                TMT6plex_color = 'palegreen' if self.files[assessor.mzml_file]['lowend_peaks'][destination]['TMT6plex']['peak']['assessment']['is_found'] else 'navajowhite'
                IH_color = 'palegreen' if self.files[assessor.mzml_file]['lowend_peaks'][destination]['IH']['peak']['assessment']['is_found'] else 'navajowhite'
                y_K_color = 'palegreen' if self.files[assessor.mzml_file]['lowend_peaks'][destination]['y{K}']['peak']['assessment']['is_found'] else 'navajowhite'

                correct_bar_TMT10_129_C = 'limegreen' if self.files[assessor.mzml_file]['lowend_peaks'][destination]['TMT10_129_C']['peak']['assessment']['is_found'] else 'orange'
                correct_bar_TMT6plex = 'limegreen' if self.files[assessor.mzml_file]['lowend_peaks'][destination]['TMT6plex']['peak']['assessment']['is_found'] else 'orange'
                correct_bar_IH = 'limegreen' if self.files[assessor.mzml_file]['lowend_peaks'][destination]['IH']['peak']['assessment']['is_found'] else 'orange'
                correct_bar_y_K = 'limegreen' if self.files[assessor.mzml_file]['lowend_peaks'][destination]['y{K}']['peak']['assessment']['is_found'] else 'orange'

                fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(16, 10.9), layout='constrained')

                # Plot TMT10_129_C
                if TMT10_129_C_extent_bins * 2 < 10:
                    ax[0,0].axvspan(xmin=TMT10_129_C_fitting_min, xmax=TMT10_129_C_fitting_max, color='lightsteelblue', alpha=0.25, lw=0)
                if TMT10_129_C_extent_bins < 10:
                    ax[0,0].axvspan(xmin=TMT10_129_C_extent_min, xmax=TMT10_129_C_extent_max, color=TMT10_129_C_color, alpha=0.75, lw=0)
                ax[0,0].axvline(x=TMT10_129_C, color=correct_bar_TMT10_129_C, alpha=0.75)
                TMT10_129_C_range = (mz >= TMT10_129_C_axis_min) & (mz <= TMT10_129_C_axis_max)
                ax[0,0].plot(mz[TMT10_129_C_range], intensities[TMT10_129_C_range])

                try:
                    ax[0,0].set_title(f"{file_name_root} ({destination})")
                except:
                    ax[0,0].set_title(f"{file_name_root} ({destination})")

                ## Plot possible fit
                little_label = ''
                try:
                    if self.files[assessor.mzml_file]['lowend_peaks'][destination]['TMT10_129_C']['peak']['assessment']['is_found']:
                        fit_param = self.files[assessor.mzml_file]['lowend_peaks'][destination]['TMT10_129_C']['peak']['fit']
                        sigma_mz = fit_param['sigma_mz']
                        mu_mz = fit_param['mz']
                        y_offset = fit_param['y_offset']
                        mz_subset = mz[TMT10_129_C_range]
                        intensity_subset = intensities[TMT10_129_C_range]

                        gaussian = norm.pdf(mz_subset, mu_mz, sigma_mz)
                        scaled_gaussian = gaussian * np.max(intensity_subset) / np.max(gaussian) + y_offset/ np.max(gaussian)

                        ax[0,0].plot(mz_subset, scaled_gaussian, color='darkred', linewidth=2)
                    else:
                        little_label = " (no fit found)"
                except:
                    pass
               
                ax[0,0].set_xlabel("m/z (TMT10_129_C)" + little_label + "\nMode of TMT10_129_C spectra: " + str(num_TMT10_129_C))
                ax[0,0].set_ylabel("Intensity")

                # Plot TMT6plex
                if TMT6plex_extent_bins * 2 < 10:
                    ax[0,1].axvspan(xmin=TMT6plex_fitting_min, xmax=TMT6plex_fitting_max, color='lightsteelblue', alpha=0.25, lw=0)
                if TMT6plex_extent_bins < 10:
                    ax[0,1].axvspan(xmin=TMT6plex_extent_min, xmax=TMT6plex_extent_max, color=TMT6plex_color, alpha=0.75, lw=0)
                ax[0,1].axvline(x=TMT6plex, color=correct_bar_TMT6plex, alpha=0.75)
                TMT6plex_range = (mz >= TMT6plex_axis_min) & (mz <= TMT6plex_axis_max)
                ax[0,1].plot(mz[TMT6plex_range], intensities[TMT6plex_range])

                try:
                    ax[0,1].set_title(f"{file_name_root} ({destination})")
                except:
                    ax[0,1].set_title(f"{file_name_root} ({destination})")

                ## Plot possible fit
                little_label = ''
                try:
                    if self.files[assessor.mzml_file]['lowend_peaks'][destination]['TMT6plex']['peak']['assessment']['is_found']:
                        fit_param = self.files[assessor.mzml_file]['lowend_peaks'][destination]['TMT6plex']['peak']['fit']
                        sigma_mz = fit_param['sigma_mz']
                        mu_mz = fit_param['mz']
                        y_offset = fit_param['y_offset']
                        mz_subset = mz[TMT6plex_range]
                        intensity_subset = intensities[TMT6plex_range]

                        gaussian = norm.pdf(mz_subset, mu_mz, sigma_mz)
                        scaled_gaussian = gaussian * np.max(intensity_subset) / np.max(gaussian) + y_offset/ np.max(gaussian)

                        ax[0,1].plot(mz_subset, scaled_gaussian, color='darkred', linewidth=2)
                    else:
                        little_label = " (no fit found)"
                except:
                    pass

                ax[0,1].set_xlabel("m/z (TMT6plex)" + little_label + "\nMode of TMT6plex spectra: " + str(num_TMT6plex))
                ax[0,1].set_ylabel("Intensity")

                # Plot IH
                if IH_extent_bins * 2 < 10:
                    ax[1,0].axvspan(xmin=IH_fitting_min, xmax=IH_fitting_max, color='lightsteelblue', alpha=0.25, lw=0)
                if IH_extent_bins < 10:
                    ax[1,0].axvspan(xmin=IH_extent_min, xmax=IH_extent_max, color=IH_color, alpha=0.75, lw=0)
                ax[1,0].axvline(x=IH, color=correct_bar_IH, alpha=0.75)
                IH_range = (mz >= IH_axis_min) & (mz <= IH_axis_max)
                ax[1,0].plot(mz[IH_range], intensities[IH_range])

                try:
                    ax[1,0].set_title(f"{file_name_root} ({destination})")
                except:
                    ax[1,0].set_title(f"{file_name_root} ({destination})")

                ## Plot possible fit
                little_label = ''
                try:
                    if self.files[assessor.mzml_file]['lowend_peaks'][destination]['IH']['peak']['assessment']['is_found']:
                        fit_param = self.files[assessor.mzml_file]['lowend_peaks'][destination]['IH']['peak']['fit']
                        sigma_mz = fit_param['sigma_mz']
                        mu_mz = fit_param['mz']
                        y_offset = fit_param['y_offset']
                        mz_subset = mz[IH_range]
                        intensity_subset = intensities[IH_range]

                        gaussian = norm.pdf(mz_subset, mu_mz, sigma_mz)
                        scaled_gaussian = gaussian * np.max(intensity_subset) / np.max(gaussian) + y_offset/ np.max(gaussian)

                        ax[1,0].plot(mz_subset, scaled_gaussian, color='darkred', linewidth=2)
                    else:
                        little_label = " (no fit found)"
                except:
                    pass

                ax[1,0].set_xlabel("m/z (IH)" + little_label + "\nMode of IH spectra: " + str(num_IH))
                ax[1,0].set_ylabel("Intensity")

                # Plot y{K}
                if y_K_extent_bins * 2 < 10:
                    ax[1,1].axvspan(xmin=y_K_fitting_min, xmax=y_K_fitting_max, color='lightsteelblue', alpha=0.25, lw=0)
                if y_K_extent_bins < 10:
                    ax[1,1].axvspan(xmin=y_K_extent_min, xmax=y_K_extent_max, color=y_K_color, alpha=0.75, lw=0)
                ax[1,1].axvline(x=y_K, color=correct_bar_y_K, alpha=0.75)
                y_K_range = (mz >= y_K_axis_min) & (mz <= y_K_axis_max)
                ax[1,1].plot(mz[y_K_range], intensities[y_K_range])

                try:
                    ax[1,1].set_title(f"{file_name_root} ({destination})")
                except:
                    ax[1,1].set_title(f"{file_name_root} ({destination})")

                ## Plot possible fit
                little_label = ''
                try:
                    if self.files[assessor.mzml_file]['lowend_peaks'][destination]['y{K}']['peak']['assessment']['is_found']:
                        fit_param = self.files[assessor.mzml_file]['lowend_peaks'][destination]['y{K}']['peak']['fit']
                        sigma_mz = fit_param['sigma_mz']
                        mu_mz = fit_param['mz']
                        y_offset = fit_param['y_offset']
                        mz_subset = mz[y_K_range]
                        intensity_subset = intensities[y_K_range]

                        gaussian = norm.pdf(mz_subset, mu_mz, sigma_mz)
                        scaled_gaussian = gaussian * np.max(intensity_subset) / np.max(gaussian) + y_offset/ np.max(gaussian)

                        ax[1,1].plot(mz_subset, scaled_gaussian, color='darkred', linewidth=2)
                    else:
                        little_label = " (no fit found)"
                except:
                    pass

                ax[1,1].set_xlabel("m/z (y{K})" + little_label + "\nMode of y{K} spectra: " + str(num_y_K))
                ax[1,1].set_ylabel("Intensity")

                ax[0,0].ticklabel_format(style='plain', axis='x', useOffset=False)
                ax[0,1].ticklabel_format(style='plain', axis='x', useOffset=False)
                ax[1,0].ticklabel_format(style='plain', axis='x', useOffset=False)
                ax[1,1].ticklabel_format(style='plain', axis='x', useOffset=False)

                ax[0,0].text(-0.07, 1.04, "a", transform=ax[0,0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')
                ax[0,1].text(-0.07, 1.04, "b", transform=ax[0,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')
                ax[1,0].text(-0.07, 1.04, "c", transform=ax[1,0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')
                ax[1,1].text(-0.07, 1.04, "d", transform=ax[1,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

                plt.tight_layout()

                # Save as SVG
                plt.savefig(f"{file_name_root}_LE_multipanel.svg", format="svg", bbox_inches='tight')

                # Full low-end composite spectrum
                #plt.plot(mz, intensities)
                #plt.xlabel("m/z")
                #plt.ylabel("Intensity")
                #lt.title(f"{file_name_root} ({destination})\nFull Low-End Spectrum")
                #plt.tight_layout()
                #plt.savefig(f"{file_name_root}_LE_full.svg", format="svg", bbox_inches='tight')
                #plt.show()