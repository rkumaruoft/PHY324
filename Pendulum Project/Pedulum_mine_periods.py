import csv

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit


# Function to convert x, y to theta in degrees
def to_theta_deg(x, y):
    return np.degrees(np.arctan2(x, -y))


def damped_cosine(t, theta_0, tau, T, phi_0):
    phi_0_rad = np.radians(phi_0)
    this_theta = theta_0 * np.exp(-t / tau) * np.cos((2 * np.pi / T) * t + phi_0_rad)
    return this_theta


# Fast CSV reading using Pandas
data = pd.read_csv('t_vs_xy_mine.csv', skiprows=1, header=None).values  # Skipping header row

# Extract columns efficiently
time, x, y = data[:, 0], data[:, 1], data[:, 2]

# Compute theta
theta = to_theta_deg(x, y)
theta_err = 1
start_theta = []
periods = []
periods_err = []

print("==========================")
valid_indices = time >= 0.36

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
plt.figure(figsize=(8, 5))
plt.plot(time_filtered, theta_filtered, label=r'$\theta$ vs. Time (Filtered)', color='b')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(theta_filtered[0]))
plt.grid(True)

"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = 60  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=time_filtered, ydata=theta_filtered,
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(time_filtered, *popt)
plt.plot(time_filtered, fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_filtered[0]) #63
print("Start_theta from curve fit = ", popt[0])
print("Period = ", popt[2])
print("tau = ", popt[1])
start_theta.append(theta_filtered[0])
periods.append(popt[2])
periods_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.show()
print("==========================")

print("==========================")
valid_indices = time >= 10.3

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
plt.figure(figsize=(8, 5))
plt.plot(time_filtered, theta_filtered, label=r'$\theta$ vs. Time (Filtered)', color='b')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(theta_filtered[0]))
plt.grid(True)

"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = 60  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=time_filtered, ydata=theta_filtered,
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(time_filtered, *popt)
plt.plot(time_filtered, fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_filtered[0]) # -58
print("Start_theta from curve fit = ", popt[0])
print("Period = ", popt[2])
print("tau = ", popt[1])
start_theta.append(theta_filtered[0])
periods.append(popt[2])
periods_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.show()
print("==========================")


print("==========================")
valid_indices = time >= 33.36

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
plt.figure(figsize=(8, 5))
plt.plot(time_filtered, theta_filtered, label=r'$\theta$ vs. Time (Filtered)', color='b')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(theta_filtered[0]))
plt.grid(True)

"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = 60  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=time_filtered, ydata=theta_filtered,
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(time_filtered, *popt)
plt.plot(time_filtered, fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_filtered[0]) # 44
print("Start_theta from curve fit = ", popt[0])
print("Period = ", popt[2])
print("tau = ", popt[1])
start_theta.append(theta_filtered[0])
periods.append(popt[2])
periods_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.show()
print("==========================")


print("==========================")
valid_indices = time >= 50.31

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
plt.figure(figsize=(8, 5))
plt.plot(time_filtered, theta_filtered, label=r'$\theta$ vs. Time (Filtered)', color='b')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(theta_filtered[0]))
plt.grid(True)

"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = 60  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=time_filtered, ydata=theta_filtered,
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(time_filtered, *popt)
plt.plot(time_filtered, fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_filtered[0]) # -41
print("Start_theta from curve fit = ", popt[0])
print("Period = ", popt[2])
print("tau = ", popt[1])
start_theta.append(theta_filtered[0])
periods.append(popt[2])
periods_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.show()
print("==========================")

print("==========================")
valid_indices = time >= 108.95

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
plt.figure(figsize=(8, 5))
plt.plot(time_filtered, theta_filtered, label=r'$\theta$ vs. Time (Filtered)', color='b')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(theta_filtered[0]))
plt.grid(True)

"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = 60  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=time_filtered, ydata=theta_filtered,
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(time_filtered, *popt)
plt.plot(time_filtered, fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_filtered[0]) #28
print("Start_theta from curve fit = ", popt[0])
print("Period = ", popt[2])
print("tau = ", popt[1])
start_theta.append(theta_filtered[0])
periods.append(popt[2])
periods_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.show()
print("==========================")

print("=======Start trusting taus from here===================")
valid_indices = time >= 186.95

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
plt.figure(figsize=(8, 5))
plt.plot(time_filtered, theta_filtered, label=r'$\theta$ vs. Time (Filtered)', color='b')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(theta_filtered[0]))
plt.grid(True)

"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = 60  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=time_filtered, ydata=theta_filtered,
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(time_filtered, *popt)
plt.plot(time_filtered, fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_filtered[0]) #-20
print("Start_theta from curve fit = ", popt[0])
print("Period = ", popt[2])
print("tau = ", popt[1], "err = ", np.sqrt(pcov[1][1]))
print("t_max = ", max(time_filtered) - min(time_filtered))
start_theta.append(theta_filtered[0])
periods.append(popt[2])
periods_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.show()
print("==========================")

print("==========================")
valid_indices = time >= 263.01

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
plt.figure(figsize=(8, 5))
plt.plot(time_filtered, theta_filtered, label=r'$\theta$ vs. Time (Filtered)', color='b')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(theta_filtered[0]))
plt.grid(True)

"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = 60  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=time_filtered, ydata=theta_filtered,
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(time_filtered, *popt)
plt.plot(time_filtered, fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_filtered[0]) #14
print("Start_theta from curve fit = ", popt[0])
print("Period = ", popt[2])
print("tau = ", popt[1], "err = ", np.sqrt(pcov[1][1]))
print("t_max = ", max(time_filtered) - min(time_filtered))
start_theta.append(theta_filtered[0])
periods.append(popt[2])
periods_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.show()
print("==========================")

print("==========================")
valid_indices = time >= 341.71

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
plt.figure(figsize=(8, 5))
plt.plot(time_filtered, theta_filtered, label=r'$\theta$ vs. Time (Filtered)', color='b')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(theta_filtered[0]))
plt.grid(True)

"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = 60  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=time_filtered, ydata=theta_filtered,
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(time_filtered, *popt)
plt.plot(time_filtered, fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_filtered[0]) #-10
print("Start_theta from curve fit = ", popt[0])
print("Period = ", popt[2])
print("tau = ", popt[1], "err = ", np.sqrt(pcov[1][1]))
print("t_max = ", max(time_filtered) - min(time_filtered))
start_theta.append(theta_filtered[0])
periods.append(popt[2])
periods_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.show()
print("==========================")

print("==========================")
valid_indices = time >= 387.74

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
plt.figure(figsize=(8, 5))
plt.plot(time_filtered, theta_filtered, label=r'$\theta$ vs. Time (Filtered)', color='b')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(theta_filtered[0]))
plt.grid(True)

"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = 60  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=time_filtered, ydata=theta_filtered,
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(time_filtered, *popt)
plt.plot(time_filtered, fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_filtered[0]) #8
print("Start_theta from curve fit = ", popt[0])
print("Period = ", popt[2])
print("tau = ", popt[1], "err = ", np.sqrt(pcov[1][1]))
print("t_max = ", max(time_filtered) - min(time_filtered))
start_theta.append(theta_filtered[0])
periods.append(popt[2])
periods_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.show()
print("==========================")

print("==========================")
valid_indices = time >= 392.71

# Filter data
time_filtered = time[valid_indices]
theta_filtered = theta[valid_indices]
plt.figure(figsize=(8, 5))
plt.plot(time_filtered, theta_filtered, label=r'$\theta$ vs. Time (Filtered)', color='b')
plt.xlabel('Time (s)')
plt.ylabel(r'$\theta$ (degrees)')
plt.title(r'$\theta$ vs. Time Starting at' + str(theta_filtered[0]))
plt.grid(True)

"""Curve fit this section"""
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = 60  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=time_filtered, ydata=theta_filtered,
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(time_filtered, *popt)
plt.plot(time_filtered, fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_filtered[0]) #-8
print("Start_theta from curve fit = ", popt[0])
print("Period = ", popt[2])
print("tau = ", popt[1], "err = ", np.sqrt(pcov[1][1]))
print("t_max = ", max(time_filtered) - min(time_filtered))
start_theta.append(theta_filtered[0])
periods.append(popt[2])
periods_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.show()
print("==========================")
periods_err = np.array([x + 0.017 for x in periods_err])
# Specify the CSV filename
csv_filename = "for_periods.csv"

# Write to CSV
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.writer(file)

    # Write header (optional)
    writer.writerow(["Start Angle", "Period", "Period Err"])

    # Write data
    for a, b, c in zip(start_theta, periods, periods_err):
        writer.writerow([a, b, c])

