import csv

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


# Function to convert x, y to theta in degrees
def to_theta_deg(x, y):
    return np.degrees(np.arctan2(x, -y))


def damped_cosine(t, theta_0, tau, T, phi_0):
    phi_0_rad = np.radians(phi_0)
    this_theta = theta_0 * np.exp(-t / tau) * np.cos((2 * np.pi / T) * t + phi_0_rad)
    return this_theta

def decay_function(t, theta_0, tau, b):
    return (theta_0 * np.exp(-t / tau)) + b


# Fast CSV reading using Pandas
data = pd.read_csv('t_vs_xy_mine.csv', skiprows=1, header=None).values  # Skipping header row

# Extract columns efficiently
time, x, y = data[:, 0], data[:, 1], data[:, 2]

# Compute theta
theta = to_theta_deg(x, y)

# Find maxima (peaks)
peaks, peak_properties = find_peaks(theta, height=0)  # Adjust 'height' if needed

# Find minima (troughs)
minima, min_properties = find_peaks(-theta, height=0)  # Finding peaks in inverted theta

# Extract peak and trough times and values
peak_times, peak_values = time[peaks], theta[peaks]
minima_times, minima_values = time[minima], theta[minima]

# Plot data with peaks and minima
plt.figure(figsize=(8, 5))
plt.plot(peak_times, peak_values, "r.", label="Peaks")  # Mark peaks in red
plt.plot(minima_times, minima_values, "g.", label="Minima")  # Mark minima in green
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time with Peaks and Minima')
plt.legend()
plt.show()

start_thetas = []
taus = []
tau_errors = []
t_max = []

print("============================")
valid_indices = time >= 0.36

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
t_max.append(np.max(time_filtered) - np.min(time_filtered))
peaks, _ = find_peaks(theta_filtered, height=0)
this_peak_times, this_peak_thetas = time_filtered[peaks], theta_filtered[peaks]
plt.figure(figsize=(8, 5))
plt.plot(this_peak_times, this_peak_thetas, label=r'$\theta$ vs. Time (Filtered)', color='b', linestyle="none",
         marker=".", markersize=2)
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(this_peak_thetas[0]))
"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = this_peak_thetas[0]  # First value in the sectioned data
b = 0
popt, pcov = curve_fit(decay_function, this_peak_times, this_peak_thetas, p0=[theta_0_guess, tau_guess, b])
fit_data = decay_function(this_peak_times, *popt)
plt.plot(this_peak_times, fit_data, linestyle="dashed", markersize=2,
         color='g', label="decay fit")
print("Start_theta =", this_peak_thetas[0])
print("tau = ", popt[1])
start_thetas.append(this_peak_thetas[0])
taus.append(popt[1])
tau_errors.append(np.sqrt(pcov[1][1]))
plt.grid(True)
plt.legend()
plt.show()
print("============================")

print("============================")
valid_indices = time >= 10.3
# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
t_max.append(np.max(time_filtered) - np.min(time_filtered))
peaks, _ = find_peaks(-theta_filtered, height=0)
this_peak_times, this_peak_thetas = time_filtered[peaks], theta_filtered[peaks]
plt.figure(figsize=(8, 5))
plt.plot(this_peak_times, this_peak_thetas, label=r'$\theta$ vs. Time (Filtered)', color='b', linestyle="none",
         marker=".", markersize=2)
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(this_peak_thetas[0]))
"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = this_peak_thetas[0]  # First value in the sectioned data
b = 0
popt, pcov = curve_fit(decay_function, this_peak_times, this_peak_thetas, p0=[theta_0_guess, tau_guess, b])
fit_data = decay_function(this_peak_times, *popt)
plt.plot(this_peak_times, fit_data, linestyle="dashed", markersize=2,
         color='g', label="decay fit")
print("Start_theta =", this_peak_thetas[0])
print("tau = ", popt[1])
start_thetas.append(this_peak_thetas[0])
taus.append(popt[1])
tau_errors.append(np.sqrt(pcov[1][1]))
plt.grid(True)
plt.legend()
plt.show()
print("============================")

print("============================")
valid_indices = time >= 33.36
# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
t_max.append(np.max(time_filtered) - np.min(time_filtered))
peaks, _ = find_peaks(theta_filtered, height=0)
this_peak_times, this_peak_thetas = time_filtered[peaks], theta_filtered[peaks]
plt.figure(figsize=(8, 5))
plt.plot(this_peak_times, this_peak_thetas, label=r'$\theta$ vs. Time (Filtered)', color='b', linestyle="none",
         marker=".", markersize=2)
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(this_peak_thetas[0]))
"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = this_peak_thetas[0]  # First value in the sectioned data
b = 0
popt, pcov = curve_fit(decay_function, this_peak_times, this_peak_thetas, p0=[theta_0_guess, tau_guess, b])
fit_data = decay_function(this_peak_times, *popt)
plt.plot(this_peak_times, fit_data, linestyle="dashed", markersize=2,
         color='g', label="decay fit")
print("Start_theta =", this_peak_thetas[0])
print("tau = ", popt[1])
start_thetas.append(this_peak_thetas[0])
taus.append(popt[1])
tau_errors.append(np.sqrt(pcov[1][1]))
plt.grid(True)
plt.legend()
plt.show()
print("============================")
print("============================")
valid_indices = time >= 50.31
# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
t_max.append(np.max(time_filtered) - np.min(time_filtered))
peaks, _ = find_peaks(-theta_filtered, height=0)
this_peak_times, this_peak_thetas = time_filtered[peaks], theta_filtered[peaks]
plt.figure(figsize=(8, 5))
plt.plot(this_peak_times, this_peak_thetas, label=r'$\theta$ vs. Time (Filtered)', color='b', linestyle="none",
         marker=".", markersize=2)
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(this_peak_thetas[0]))
"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = this_peak_thetas[0]  # First value in the sectioned data
b = 0
popt, pcov = curve_fit(decay_function, this_peak_times, this_peak_thetas, p0=[theta_0_guess, tau_guess, b])
fit_data = decay_function(this_peak_times, *popt)
plt.plot(this_peak_times, fit_data, linestyle="dashed", markersize=2,
         color='g', label="decay fit")
print("Start_theta =", this_peak_thetas[0])
print("tau = ", popt[1])
start_thetas.append(this_peak_thetas[0])
taus.append(popt[1])
tau_errors.append(np.sqrt(pcov[1][1]))
plt.grid(True)
plt.legend()
plt.show()
print("============================")
print("============================")
valid_indices = time >= 108.95
# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
t_max.append(np.max(time_filtered) - np.min(time_filtered))
peaks, _ = find_peaks(theta_filtered, height=0)
this_peak_times, this_peak_thetas = time_filtered[peaks], theta_filtered[peaks]
plt.figure(figsize=(8, 5))
plt.plot(this_peak_times, this_peak_thetas, label=r'$\theta$ vs. Time (Filtered)', color='b', linestyle="none",
         marker=".", markersize=2)
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(this_peak_thetas[0]))
"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = this_peak_thetas[0]  # First value in the sectioned data
b = 0
popt, pcov = curve_fit(decay_function, this_peak_times, this_peak_thetas, p0=[theta_0_guess, tau_guess, b])
fit_data = decay_function(this_peak_times, *popt)
plt.plot(this_peak_times, fit_data, linestyle="dashed", markersize=2,
         color='g', label="decay fit")
print("Start_theta =", this_peak_thetas[0])
print("tau = ", popt[1])
start_thetas.append(this_peak_thetas[0])
taus.append(popt[1])
tau_errors.append(np.sqrt(pcov[1][1]))
plt.grid(True)
plt.legend()
plt.show()
print("============================")
print("============================")
valid_indices = time >= 186.95
# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
t_max.append(np.max(time_filtered) - np.min(time_filtered))
peaks, _ = find_peaks(-theta_filtered, height=0)
this_peak_times, this_peak_thetas = time_filtered[peaks], theta_filtered[peaks]
plt.figure(figsize=(8, 5))
plt.plot(this_peak_times, this_peak_thetas, label=r'$\theta$ vs. Time (Filtered)', color='b', linestyle="none",
         marker=".", markersize=2)
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(this_peak_thetas[0]))
"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = this_peak_thetas[0]  # First value in the sectioned data
b = 0
popt, pcov = curve_fit(decay_function, this_peak_times, this_peak_thetas, p0=[theta_0_guess, tau_guess, b])
fit_data = decay_function(this_peak_times, *popt)
plt.plot(this_peak_times, fit_data, linestyle="dashed", markersize=2,
         color='g', label="decay fit")
print("Start_theta =", this_peak_thetas[0])
print("tau = ", popt[1])
start_thetas.append(this_peak_thetas[0])
taus.append(popt[1])
tau_errors.append(np.sqrt(pcov[1][1]))
plt.grid(True)
plt.legend()
plt.show()
print("============================")

print("============================")
valid_indices = time >= 263.01
# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
t_max.append(np.max(time_filtered) - np.min(time_filtered))
peaks, _ = find_peaks(theta_filtered, height=0)
this_peak_times, this_peak_thetas = time_filtered[peaks], theta_filtered[peaks]
plt.figure(figsize=(8, 5))
plt.plot(this_peak_times, this_peak_thetas, label=r'$\theta$ vs. Time (Filtered)', color='b', linestyle="none",
         marker=".", markersize=2)
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(this_peak_thetas[0]))
"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = this_peak_thetas[0]  # First value in the sectioned data
b = 0
popt, pcov = curve_fit(decay_function, this_peak_times, this_peak_thetas, p0=[theta_0_guess, tau_guess, b])
fit_data = decay_function(this_peak_times, *popt)
plt.plot(this_peak_times, fit_data, linestyle="dashed", markersize=2,
         color='g', label="decay fit")
print("Start_theta =", this_peak_thetas[0])
print("tau = ", popt[1])
start_thetas.append(this_peak_thetas[0])
taus.append(popt[1])
tau_errors.append(np.sqrt(pcov[1][1]))
plt.grid(True)
plt.legend()
plt.show()
print("============================")
print("============================")
valid_indices = time >= 341.71
# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
t_max.append(np.max(time_filtered) - np.min(time_filtered))
peaks, _ = find_peaks(-theta_filtered, height=0)
this_peak_times, this_peak_thetas = time_filtered[peaks], theta_filtered[peaks]
plt.figure(figsize=(8, 5))
plt.plot(this_peak_times, this_peak_thetas, label=r'$\theta$ vs. Time (Filtered)', color='b', linestyle="none",
         marker=".", markersize=2)
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(this_peak_thetas[0]))
"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = this_peak_thetas[0]  # First value in the sectioned data
b = 0
popt, pcov = curve_fit(decay_function, this_peak_times, this_peak_thetas, p0=[theta_0_guess, tau_guess, b])
fit_data = decay_function(this_peak_times, *popt)
plt.plot(this_peak_times, fit_data, linestyle="dashed", markersize=2,
         color='g', label="decay fit")
print("Start_theta =", this_peak_thetas[0])
print("tau = ", popt[1])
start_thetas.append(this_peak_thetas[0])
taus.append(popt[1])
tau_errors.append(np.sqrt(pcov[1][1]))
plt.grid(True)
plt.legend()
plt.show()
print("============================")
print("============================")
valid_indices = time >= 387.74
# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
t_max.append(np.max(time_filtered) - np.min(time_filtered))
peaks, _ = find_peaks(theta_filtered, height=0)
this_peak_times, this_peak_thetas = time_filtered[peaks], theta_filtered[peaks]
plt.figure(figsize=(8, 5))
plt.plot(this_peak_times, this_peak_thetas, label=r'$\theta$ vs. Time (Filtered)', color='b', linestyle="none",
         marker=".", markersize=2)
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(this_peak_thetas[0]))
"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = this_peak_thetas[0]  # First value in the sectioned data
b = 0
popt, pcov = curve_fit(decay_function, this_peak_times, this_peak_thetas, p0=[theta_0_guess, tau_guess, b],
                       maxfev=10000)
fit_data = decay_function(this_peak_times, *popt)
plt.plot(this_peak_times, fit_data, linestyle="dashed", markersize=2,
         color='g', label="decay fit")
print("Start_theta =", this_peak_thetas[0])
print("tau = ", popt[1])
start_thetas.append(this_peak_thetas[0])
taus.append(popt[1])
tau_errors.append(np.sqrt(pcov[1][1]))
plt.grid(True)
plt.legend()
plt.show()
print("============================")

print("============================")
valid_indices = time >= 392.71
# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
t_max.append(np.max(time_filtered) - np.min(time_filtered))
peaks, _ = find_peaks(-theta_filtered, height=0)
this_peak_times, this_peak_thetas = time_filtered[peaks], theta_filtered[peaks]
plt.figure(figsize=(8, 5))
plt.plot(this_peak_times, this_peak_thetas, label=r'$\theta$ vs. Time (Filtered)', color='b', linestyle="none",
         marker=".", markersize=2)
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(this_peak_thetas[0]))
"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = this_peak_thetas[0]  # First value in the sectioned data
b = 0
popt, pcov = curve_fit(decay_function, this_peak_times, this_peak_thetas, p0=[theta_0_guess, tau_guess, b],
                       maxfev=10000)
fit_data = decay_function(this_peak_times, *popt)
plt.plot(this_peak_times, fit_data, linestyle="dashed", markersize=2,
         color='g', label="decay fit")
print("Start_theta =", this_peak_thetas[0])
print("tau = ", popt[1])
start_thetas.append(this_peak_thetas[0])
taus.append(popt[1])
tau_errors.append(np.sqrt(pcov[1][1]))
plt.grid(True)
plt.legend()
plt.show()
print("============================")

# Specify the CSV filename
csv_filename = "for_taus.csv"

# Write to CSV
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.writer(file)

    # Write header (optional)
    writer.writerow(["Start Angle", "Tau", "Tau Err", "T max"])

    # Write data
    for a, b, c, d in zip(start_thetas, taus, tau_errors, t_max):
        writer.writerow([a, b, c, d])

