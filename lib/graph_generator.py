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

    def buildGraphs(self):

        self.files = self.data.get("files", {})

        coverpage = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Delta_graphs_documentation.pdf")

        #### Output PDF file
        output_pdf_path = self.metadata_file.replace(".json", ".histograms_with_chi2.pdf")
        with PdfPages(output_pdf_path) as pdf:

            for filename, file_data in self.files.items():
                precursor_stats = file_data.get("summary", {}).get('precursor stats', {})
                dynamic_exclusion_window = precursor_stats.get('dynamic exclusion window', {})

                try:
                    pre_tol = precursor_stats.get('precursor tolerance', {})
                    hist_time = dynamic_exclusion_window.get("histogram_time", {})
                    hist_ppm = pre_tol.get("histogram_ppm", {})

                    amplitude_ppm = pre_tol.get('fit_ppm', {}).get("intensity", [])
                    mu_ppm = pre_tol.get('fit_ppm', {}).get("delta_ppm peak", [])
                    sigma_ppm = pre_tol.get('fit_ppm', {}).get("sigma (ppm)", [])
                    y_offset_ppm = pre_tol.get('fit_ppm', {}).get("y_offset", [])

                    time_initial_level = dynamic_exclusion_window.get("fit_pulse_time", {}).get("inital level", [])
                    pulse_start = dynamic_exclusion_window.get("fit_pulse_time", {}).get("pulse start", [])
                    pulse_duration = dynamic_exclusion_window.get("fit_pulse_time", {}).get("dynamic exclusion time offset", [])
                    peak_level = dynamic_exclusion_window.get("fit_pulse_time", {}).get("peak level", [])
                    decay_constant = dynamic_exclusion_window.get("fit_pulse_time", {}).get("decay constant", [])
                    final_level = dynamic_exclusion_window.get("fit_pulse_time", {}).get("final level", [])

                    chi_squared_time = hist_time.get("chi_squared", "N/A")
                    chi_squared_ppm = hist_ppm.get("chi_squared", "N/A")
                except:
                    continue

                if not hist_time or not hist_ppm:
                    eprint(f"Skipping {filename} due to missing histogram data")
                    continue

                time_bin_centers = np.array(hist_time.get("bin_centers", []))
                time_counts = np.array(hist_time.get("counts", []))

                ppm_bin_centers = np.array(hist_ppm.get("bin_centers", []))
                ppm_counts = np.array(hist_ppm.get("counts", []))

                if time_bin_centers.size == 0 or time_counts.size == 0 or ppm_bin_centers.size == 0 or ppm_counts.size == 0:
                    print(f"Incomplete histogram data in {filename}, skipping...")
                    continue

                fig, axs = plt.subplots(2, 2, figsize=(12, 8), gridspec_kw={'height_ratios': [4, 1]})
                fig.suptitle(filename, fontsize=14)

                #### Time histogram
                ax_time = axs[0, 0]
                width_time = (time_bin_centers[1] - time_bin_centers[0]) if len(time_bin_centers) > 1 else 1
                ax_time.bar(time_bin_centers, time_counts, width=width_time, color='skyblue', edgecolor='black', label='Time Data')

                try:
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

                ax_time.set_title("Time Between Precursors")
                ax_time.set_xlabel("Time (seconds)")
                ax_time.set_ylabel("Counts")
                ax_time.legend()

                #### PPM histogram 
                ax_ppm = axs[0, 1]
                width_ppm = (ppm_bin_centers[1] - ppm_bin_centers[0]) if len(ppm_bin_centers) > 1 else 1
                ax_ppm.bar(ppm_bin_centers, ppm_counts, width=width_ppm, color='salmon', edgecolor='black', label='PPM Data')

                try:
                    x_ppm = np.linspace(min(ppm_bin_centers), max(ppm_bin_centers), 500)
                    scale = amplitude_ppm * sigma_ppm * sqrt(2 * pi)
                    y_ppm = scale * norm.pdf(x_ppm, mu_ppm, sigma_ppm) + y_offset_ppm
                    ax_ppm.plot(x_ppm, y_ppm, color='darkred', linewidth=2, label='PPM Gaussian Fit')
                except:
                    pass

                ax_ppm.set_title("Precursor ppm Variance")
                ax_ppm.set_xlabel("PPM Varience")
                ax_ppm.set_ylabel("Counts")
                ax_ppm.legend()

                #### Tables
                axs[1, 0].axis('off')
                axs[1, 1].axis('off')

                def fmt(v):
                    return f"{v:.3f}" if isinstance(v, (int, float)) else str(v)

                ## Table under Time histogram
                time_table_data = [
                    ["Chi² Time", fmt(chi_squared_time)],
                    ["Initial Level", fmt(time_initial_level)],
                    ["Pulse Start", fmt(pulse_start)],
                    ["Pulse Duration", fmt(pulse_duration)],
                    ["Peak Level", fmt(peak_level)],
                    ["Decay Constant", fmt(decay_constant)],
                    ["Final Level", fmt(final_level)],
                ]
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
                ppm_table_data = [
                    ["Chi² PPM", fmt(chi_squared_ppm)],
                    ["Amplitude", fmt(amplitude_ppm)],
                    ["Mu (Peak Center)", fmt(mu_ppm)],
                    ["Sigma", fmt(sigma_ppm)],
                    ["Y Offset", fmt(y_offset_ppm)],
                ]
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

        #Add explanitory Cover Page
        writer = PdfWriter()
        cover_reader = PdfReader(coverpage)
        for page in cover_reader.pages:
            writer.add_page(page)

        # Add your generated PDF pages
        output_reader = PdfReader(output_pdf_path)
        for page in output_reader.pages:
            writer.add_page(page)

        # Write the final combined PDF
        with open(output_pdf_path, "wb") as f:
            writer.write(f)

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
                        eprint(f"No phosphoric acid z =2 peak found in {assessor.mzml_file}")
                        little_label = " (no fit found)"
                except:
                    pass

                plt.xlabel("m/z loss (phosphoric acid z=2)" + little_label + "\nMode of phosphoric acid z=2 spectra: " + str(num_phospho))
                plt.ylabel("Intensity")
                try:
                    plt.title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz:.4f}\nPhosphoric acid to water intensity ratio: {ratio:.2f}", fontsize=8.5)
                except: 
                    plt.title(f"{file_name_root} ({destination})\nAbsolute difference between delta m/z phosphoric acid and water: {abs_diff_delta_mz}\nPhosphoric acid to water intensity ratio: {ratio:.2f}", fontsize=8.5)
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
        nl_pdf= self.metadata_file.replace(".json", ".NLplots.pdf")
        if self.verbose >= 1:
            eprint(f"All {assessor.mzml_file} neutral loss spectra saved to: {nl_pdf}")
