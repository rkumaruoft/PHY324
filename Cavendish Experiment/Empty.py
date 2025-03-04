import math
import pickle
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2
from datetime import datetime
import csv


def damped_cosine(t, position_0, tau, omega, phase, base):
    return position_0 * np.exp(-t / tau) * np.cos((omega * t) + phase) + base


if __name__ == "__main__":
    time = []
    position = []
    position_uncert = []
    file = open('data/data_IN_csv/Empty_Feb7th.csv', mode='r')
    # Create a CSV reader object
    csv_reader = csv.reader(file)
    # Skip the header row (if there is one)
    next(csv_reader, None)
    for row in csv_reader:
        if float(row[1]) > 0.01:
            time.append(datetime.strptime(row[0], '%H:%M:%S.%f'))
            position.append(float(row[1]))
            position_uncert.append(0.001)
    fig, ax = plt.subplots(figsize=(8, 4))  # Create figure and axis
    # ax.scatter(time, position, marker=".", s=5, label="Raw Data")
    ax.errorbar(time, position, yerr=position_uncert, marker=".", markersize=5, label="Raw Data", alpha=0.5)

    # Format the x-axis (must use `ax.xaxis`, not `plt.xaxis`)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))  # Custom format
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=5))  # Major ticks every 5 minutes
    ax.xaxis.set_minor_locator(mdates.SecondLocator(interval=30))  # Minor ticks every 30 seconds

    """Fitting a damped cosine"""
    t_numeric = np.array([(t - time[0]).total_seconds() for t in time])


    # Estimate reasonable initial parameters from data
    amplitude_guess = max(position) - min(position)  # Approximate amplitude
    tau_guess = max(t_numeric) / 5  # Rough estimate of decay time
    omega_guess = 2 * np.pi / (max(t_numeric) / 10)  # Estimate period from observed data
    phase_guess = 0  # Assume initial phase 0
    base_guess = np.mean(position)  # Baseline position

    # Updated initial guesses
    p0 = [amplitude_guess, tau_guess, omega_guess, phase_guess, base_guess]

    popt, pcov = curve_fit(damped_cosine, xdata=t_numeric, ydata=position, sigma=position_uncert, p0=p0, maxfev=100000)
    position_fit = damped_cosine(t_numeric, *popt)
    ax.plot(time, position_fit, color="red", label="Fit Line")
    plt.xlabel("Time")
    plt.ylabel("Position")
    plt.xticks(rotation=45)  # Rotate labels to avoid overlap
    plt.grid(True)

    """From Curve fit"""
    print("Position Amplitude (p0) = ", popt[0], "err=", np.sqrt(pcov[0][0]))
    print("Decay Factor (tau) = ", popt[1], "err=", np.sqrt(pcov[1][1]))
    print("Angular Frequency (omega) = ", popt[2], "err=", np.sqrt(pcov[2][2]))
    print("Wave Phase (phi) = ", popt[3], "err=", np.sqrt(pcov[3][3]))
    print("Equilibrium position for Empty (base)= ", popt[4], "err=", np.sqrt(pcov[4][4]))

    print("Period (Empty) = ", 2*np.pi/popt[2])

    plt.show()
