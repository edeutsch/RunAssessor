import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from math import sqrt, pi
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

class GraphGenerator:
    def __init__(self, metadata_filepath=None, verbose=None):
        default_metadata_filepath = 'study_metadata.json'
        if metadata_filepath == None:
            self.metadata_file = default_metadata_filepath
        else:
            self.metadata_file = metadata_filepath


    def buildGraphs(self):

        #### Load JSON data
        try:
            with open(self.metadata_file, "r") as f:
                data = json.load(f)
        except:
            eprint("INFO: Unable to read JSON file")
            return

        files = data.get("files", {})

        #### Output PDF file
        output_pdf_path = "histograms_with_chi2.pdf"
        with PdfPages(output_pdf_path) as pdf:
            for filename, file_data in files.items():
                precursor_stats = file_data.get("summary", {}).get('precursor stats', {})
                dynamic_exclusion_window = precursor_stats.get('dynamic exclusion window', {})
                eprint()

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
                    print(f"Skipping {filename} due to missing histogram data")
                    continue

                time_bin_centers = np.array(hist_time.get("bin_centers", []))
                time_counts = np.array(hist_time.get("counts", []))

                ppm_bin_centers = np.array(hist_ppm.get("bin_centers", []))
                ppm_counts = np.array(hist_ppm.get("counts", []))

                if time_bin_centers.size == 0 or time_counts.size == 0 or ppm_bin_centers.size == 0 or ppm_counts.size == 0:
                    print(f"Incomplete histogram data in {filename}, skipping...")
                    continue

                fig, axs = plt.subplots(2, 2, figsize=(12, 8), gridspec_kw={'height_ratios': [4, 1]})
                fig.suptitle(f"File: {filename}", fontsize=14)

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

                ax_time.set_title("Time Histogram")
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

                ax_ppm.set_title("PPM Histogram")
                ax_ppm.set_xlabel("PPM Error")
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
                    ["Mu (Peak)", fmt(mu_ppm)],
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

                print(f"Saved plot for {filename}")

        print("Saving PDF to:", os.path.abspath(output_pdf_path))
