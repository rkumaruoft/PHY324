import numpy as np
from scipy.signal import find_peaks
from functions import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2
import csv

if __name__ == "__main__":
    x_uncert = 10
    y_uncert = 10

    time = []
    thetas = []
    thetas_uncert = []
    file = open('Data_files/data_IN_csv/200g_198_DATA_remove_after_102sec.csv', mode='r')
    # Create a CSV reader object
    csv_reader = csv.reader(file)
    # Skip the header row (if there is one)
    next(csv_reader, None)
    for row in csv_reader:
        time.append(float(row[0]))
        x, y = float(row[1]), float(row[2])
        thetas.append(calculate_theta_from_neg_y(x, y))
        thetas_uncert.append(propagate_uncertainty(x, y, 10, 10))

    """Convert data to theta"""
    plt.errorbar(time, thetas, yerr=thetas_uncert, marker=".", linestyle="None",
                 markersize=5, label="Raw Data", alpha=0.2, ecolor="yellow")
    plt.show()

    # Find the peaks in theta data to determine oscillations
    peaks, _ = find_peaks(thetas)
    valleys, _ = find_peaks(-np.array(thetas))
    # Peaks
    peaks = [0] + list(peaks)

    # Now section off the data into groups based on 5 peaks
    all_sections = []
    no_of_oscillations = 10
    for i in range(0, len(peaks), no_of_oscillations):
        # Define the start and end indices for each section (from peak i to peak i+5)
        start_peak = peaks[i]
        end_peak = peaks[i + no_of_oscillations - 1] if i + no_of_oscillations - 1 < len(peaks) else peaks[-1]
        # Extract the data for this section (from start_peak to end_peak)
        section_times = time[start_peak:end_peak + 1]  # Include the data from start to end peak
        section_thetas = thetas[start_peak:end_peak + 1]
        section_uncertainties = thetas_uncert[start_peak:end_peak + 1]
        # Store the data for each section
        all_sections.append((section_times, section_thetas, section_uncertainties))

    for i in range(0, len(valleys), no_of_oscillations):
        # Define the start a nd end indices for each section (from peak i to peak i+5)
        start_valley = valleys[i]
        end_valley = valleys[i + no_of_oscillations - 1] if i + no_of_oscillations - 1 < len(valleys) else valleys[-1]
        # Extract the data for this section (from start_peak to end_peak)
        section_times = time[start_valley:end_valley + 1]  # Include the data from start to end peak
        section_thetas = thetas[start_valley:end_valley + 1]
        section_uncertainties = thetas_uncert[start_valley:end_valley + 1]
        # Store the data for each section
        all_sections.append((section_times, section_thetas, section_uncertainties))

    periods = []
    taus = []
    periods_err = []
    taus_err = []
    start_theta = []

    for i in range(len(all_sections)):
        print(i)
        min_time = min(all_sections[i][0])
        this_time = np.array(all_sections[i][0]) - min_time
        this_thetas = np.array(all_sections[i][1])
        this_theta_uncert = np.array(all_sections[i][2])
        # A: amplitude (max value of the data)
        A_guess = np.max(this_thetas)
        omega_guess = 2 * np.pi / 1.2  # Angular frequency from period
        # phi: phase shift (initially assume 0)
        phi_guess = 0
        # tau: damping factor (initially guess 1)
        tau_guess = 100
        # B: base value (initial guess is the mean of the data)
        B_guess = np.mean(this_thetas)
        p0 = [A_guess, omega_guess, phi_guess, tau_guess, B_guess]
        popt, pcov = fit_damped_cosine(this_time, this_thetas, this_theta_uncert, p0)
        plt.errorbar(this_time, this_thetas, yerr=this_theta_uncert, marker=".", linestyle="None",
                     markersize=5, label="Raw Data", alpha=0.2, ecolor="yellow")
        plt.show()
        this_y = damped_cosine(this_time, *popt)

        periods.append(2 * np.pi / popt[1])
        periods_err.append((2 * np.pi / (popt[1] ** 2)) * np.sqrt(pcov[1][1]))
        taus.append(popt[3])
        taus_err.append(np.sqrt(pcov[3][3]))
        start_theta.append(this_thetas[0])

    # Specify the CSV filename
    csv_filename = "periods_and_taus_200_198.csv"

    # Write to CSV
    with open(csv_filename, mode='w', newline='') as file:
        writer = csv.writer(file)

        # Write header (optional)
        writer.writerow(["Start Angle", "Period", "Period Err", "Tau", "Tau Err"])

        # Write data
        for a, b, c, d, e in zip(start_theta, periods, periods_err, taus, taus_err):
            writer.writerow([a, b, c, d, e])
