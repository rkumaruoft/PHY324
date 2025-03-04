import math
import pickle
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2
import csv


def get_last_digit_error(arr):
    """
    Given an array of floating-point numbers, this function returns their
    last-digit uncertainty based on the least significant digit.

    Parameters:
    arr (list or np.array): Array of floating-point numbers.

    Returns:
    list of tuples: Each tuple contains (original_number, last_digit_error)
    """
    results = []

    for num in arr:
        # Convert to string and check number of decimal places
        num_str = f"{num:.15g}"  # Use general format to avoid scientific notation
        if '.' in num_str:
            decimal_places = len(num_str.split('.')[1])  # Count decimal places
            uncertainty = 10 ** (-decimal_places)  # Last digit uncertainty
        else:
            uncertainty = 1  # If no decimal, uncertainty is 1

        results.append(uncertainty)

    return results


def Gauss(x, A, mean, width, base):
    return A * np.exp(-(x - mean) ** 2 / (2 * width ** 2)) + base


def theta_to_lambda2(theta_angle, A, B):
    n_2 = ((((2 / math.sqrt(3)) * math.sin(np.radians(theta_angle))) + (1 / 2)) ** 2) + 3 / 4
    return A / ((np.sqrt(n_2)) - B)


angle = []
intensity = []
wavelength = []
with open('data/data_IN_csv/4V_1.csv', mode='r') as file:
    # Create a CSV reader object
    csv_reader = csv.reader(file)
    # Skip the header row (if there is one)
    next(csv_reader, None)
    # Iterate over each row in the CSV file
    for row in csv_reader:
        angle.append(float(row[0]))
        intensity.append(float(row[1]))

# Compute the error as the place value of the last significant digit
intensity_err = get_last_digit_error(intensity)

plt.xlabel("Sensor Position (Degrees)")
plt.ylabel("Intensity (Volts)")
plt.scatter(angle, intensity, label="Raw data", marker='.')
plt.legend()
plt.savefig("Graphs/raw_data.png", dpi=200)
plt.show()

small_angle_sec = []
small_intensity_sec = []
for i in range(len(angle)):
    if 76 <= angle[i] <= 82:
        small_angle_sec.append(angle[i])
        small_intensity_sec.append(intensity[i])

plt.scatter(small_angle_sec, small_intensity_sec, s=2, label="Small Peak data")
# Convert data to numpy arrays for curve fitting
x_data = np.array(small_angle_sec)
y_data = np.array(small_intensity_sec)

# Initial guesses for Gaussian parameters: amplitude, mean, and standard deviation
initial_guess = [max(y_data), np.mean(x_data), 1, 0.03]

# Fit the Gaussian model to the data
popt, pcov = curve_fit(Gauss, x_data, y_data, p0=initial_guess)

small_angle_sec = np.array(small_angle_sec)
y_fit = Gauss(small_angle_sec, *popt)

# Plot the data and the Gaussian fit
plt.plot(small_angle_sec, y_fit, '-', label=f'Gaussian Fit', linewidth=2, color="red")
plt.xlabel("Sensor Position (degrees)")
plt.ylabel("Intensity (Volts)")
plt.legend()
small_peak_at = popt[1]
small_peak_err = np.sqrt(pcov[1][1])
print("Small Peak mean = ", small_peak_at, "error = ", small_peak_err)
plt.savefig("Graphs/4V Small_Peak.png", dpi=200)
plt.show()

"""Residual for small peak fit"""
residuals_small_peak = y_data - y_fit
plt.scatter(small_angle_sec, residuals_small_peak, label='Residuals', s=2)
plt.xlabel("Sensor Position (degrees)")
plt.ylabel("Intensity Residuals (Volts)")
plt.axhline(y=0)
plt.legend()
plt.savefig("Graphs/small_peak_residual.png", dpi=200)
plt.show()

theta_init = small_peak_at
angle = np.array(angle) - theta_init

angle = -angle

plt.scatter(angle, intensity, label="Raw Data Calibrated", s=2)
plt.xlabel("Emergence angle (Theta) (degrees)")
plt.ylabel("Intensity (Volts)")
plt.legend()
plt.savefig("Graphs/4V raw_calibrated.png", dpi=200)
plt.show()

peak_angle_sec = []
peak_intensity_sec = []
for i in range(len(angle)):
    if 52 <= angle[i] <= 62:
        peak_angle_sec.append(angle[i])
        peak_intensity_sec.append(intensity[i])

plt.scatter(peak_angle_sec, peak_intensity_sec, s=2, label="Calibrated Data")

from error_propogation import *

popt, pcov = curve_fit(Gauss, peak_angle_sec, peak_intensity_sec, p0=(0.43, 57, 1, 0.05))

fit_data = Gauss(np.array(peak_angle_sec), *popt)
plt.plot(peak_angle_sec, fit_data, label="Gaussian Fit", c='k')
plt.xlabel("Emergence angle (Theta) (Degrees)")
plt.ylabel("Intensity (Volts)")
spectrum_peak_at = popt[1]
spectrum_peak_err = np.sqrt(pcov[1][1])
peak_wavelength = np.sqrt(theta_to_lambda2(spectrum_peak_at, 13900, 1.689))
error_in_peak = error_in_lambda(compute_lambda2_error(spectrum_peak_at, spectrum_peak_err, 13900, 1.689),
                                peak_wavelength)
print("Spectrum Peak mean = ", spectrum_peak_at, "err = ", spectrum_peak_err)
print("Peak wavelength = ", peak_wavelength, " nm", "err = ", error_in_peak)
plt.legend()
plt.savefig("Graphs/4V spectrum_peak.png", dpi=200)
plt.show()

"""Residual for big peak fit"""
residuals_big_peak = peak_intensity_sec - fit_data
plt.scatter(peak_angle_sec, residuals_big_peak, label='Residuals', s=2)
plt.xlabel("Emergence angle (Theta) (degrees)")
plt.ylabel("Intensity Residuals (Volts)")
plt.axhline(y=0)
plt.legend()
plt.savefig("Graphs/spectrum_residual.png", dpi=200)
plt.show()
