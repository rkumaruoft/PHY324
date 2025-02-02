import math
import pickle
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2
import csv


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

small_angle_sec = []
small_intensity_sec = []
for i in range(len(angle)):
    if 76 <= angle[i] <= 82:
        small_angle_sec.append(angle[i])
        small_intensity_sec.append(intensity[i])

plt.scatter(small_angle_sec, small_intensity_sec, marker=2)
# Convert data to numpy arrays for curve fitting
x_data = np.array(small_angle_sec)
y_data = np.array(small_intensity_sec)

# Initial guesses for Gaussian parameters: amplitude, mean, and standard deviation
initial_guess = [max(y_data), np.mean(x_data), 1, 0.03]

# Fit the Gaussian model to the data
popt, pcov = curve_fit(Gauss, x_data, y_data, p0=initial_guess)

y_fit = Gauss(small_angle_sec, *popt)

# Plot the data and the Gaussian fit
plt.plot(small_angle_sec, y_fit, '-', label=f'Gaussian Fit', linewidth=2)
plt.title("Small Peak")
plt.xlabel("Angle (degrees)")
plt.ylabel("Intensity")
plt.legend()
small_peak_at = popt[1]
small_peak_err = np.sqrt(pcov[1][1])
print("Small Peak mean = ", small_peak_at, "error = ", small_peak_err)
plt.show()

theta_init = small_peak_at
angle = np.array(angle) - theta_init

angle = -angle

plt.title("Angle vs Intensity")
plt.plot(angle, intensity)
plt.show()

peak_angle_sec = []
peak_intensity_sec = []
for i in range(len(angle)):
    if 52 <= angle[i] <= 62:
        peak_angle_sec.append(angle[i])
        peak_intensity_sec.append(intensity[i])


plt.scatter(peak_angle_sec, peak_intensity_sec, marker=1)

popt, pcov = curve_fit(Gauss, peak_angle_sec, peak_intensity_sec, p0=(0.43, 57, 1, 0.05))
plt.plot(peak_angle_sec, Gauss(np.array(peak_angle_sec), *popt))
plt.title("spectrum peak")
plt.xlabel("Angle (degrees)")
plt.ylabel("Intensity")
spectrum_peak_at = popt[1]
spectrum_peak_err = np.sqrt(pcov[1][1])
print("Spectrum Peak mean = ", spectrum_peak_at, "err = ", spectrum_peak_err)
# plt.axvline(x=spectrum_peak_at)
print("Peak wavelength = ", np.sqrt(theta_to_lambda2(spectrum_peak_at, 13900, 1.689)), "nm")
plt.show()
