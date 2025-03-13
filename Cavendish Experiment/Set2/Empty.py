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
    file = open('../data/data_IN_csv/Device2Empty.csv', mode='r')
    # Create a CSV reader object
    csv_reader = csv.reader(file)
    # Skip the header row (if there is one)
    next(csv_reader, None)
    for row in csv_reader:
        if float(row[1]) < 0.08:
            time.append(datetime.strptime(row[0], '%H:%M:%S.%f'))
            position.append(float(row[1]))
            position_uncert.append(0.00539 / 2)

    t_numeric = np.array([(t - time[0]).total_seconds() for t in time])

    fig, ax = plt.subplots(figsize=(8, 4))  # Create figure and axis
    # ax.scatter(time, position, marker=".", s=5, label="Raw Data")
    ax.errorbar(t_numeric, position, yerr=position_uncert, marker=".", linestyle="None",
                color='red', ecolor="lightblue", markersize=5, label="Raw Data", alpha=0.2)

    """Fitting a damped cosine"""


    # Estimate reasonable initial parameters from data
    amplitude_guess = max(position) - min(position)  # Approximate amplitude
    tau_guess = max(t_numeric) / 5  # Rough estimate of decay time
    omega_guess = 2 * np.pi / (max(t_numeric) / 10)  # Estimate period from observed data
    phase_guess = 0  # Assume initial phase 0
    base_guess = np.mean(position)  # Baseline position

    # Updated initial guesses
    p0 = [amplitude_guess, tau_guess, omega_guess, phase_guess, base_guess]

    popt, pcov = curve_fit(damped_cosine, xdata=t_numeric, ydata=position, sigma=position_uncert, p0=p0, maxfev=100000)
    position_fit = np.array(damped_cosine(t_numeric, *popt))
    ax.plot(t_numeric, position_fit, color="black", linestyle="--", label="Fit Line")


    """From Curve fit"""

    print("Position Amplitude (p0) = ", popt[0], "err=", np.sqrt(pcov[0][0]))
    print("Decay Factor (tau) = ", popt[1], "err=", np.sqrt(pcov[1][1]))
    print("Angular Frequency (omega) = ", popt[2], "err=", np.sqrt(pcov[2][2]))
    print("Wave Phase (phi) = ", popt[3], "err=", np.sqrt(pcov[3][3]))
    print("Equilibrium position for Empty (base)= ", popt[4], "err=", np.sqrt((np.sqrt(pcov[4][4]) ** 2) +
                                                                              ((0.00539/2) ** 2)))

    T_d = 2 * np.pi / popt[2]
    sigma_T_d = (2 * np.pi / popt[2] ** 2) * np.sqrt(pcov[2][2])

    T_err = np.sqrt(sigma_T_d ** 2 + (1/np.sqrt(12) ** 2))
    print(f"Period (Empty) = {T_d} ± {T_err}")


    # Compute damping ratio
    zeta = 1 / (popt[2] * popt[1])

    # Compute corrected undamped period
    T_0 = T_d / np.sqrt(1 - zeta ** 2)

    # Propagate uncertainty for T_0
    sigma_T_0 = sigma_T_d / np.sqrt(1 - zeta ** 2)  # Approximation for small zeta

    print(f"Undamped Period (T0) = {T_0} ± {sigma_T_0}")

    plt.xlabel("Time (seconds)")
    plt.ylabel("Position (meters)")
    plt.grid(True)
    plt.legend()
    plt.savefig("Graphs/empty_fit.png", dpi=200)
    plt.show()

    """Residuals"""
    # Compute residuals
    residuals = position - position_fit

    # Compute reduced chi-squared
    degrees_of_freedom = len(position) - len(popt)
    chi_squared = np.sum((residuals / position_uncert) ** 2)
    reduced_chi_squared = chi_squared / degrees_of_freedom
    
    # Compute fit probability (p-value)
    p_value = 1 - chi2.cdf(chi_squared, degrees_of_freedom)

    # Plot residuals
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.errorbar(t_numeric, residuals, yerr=position_uncert, marker=".", markersize=5, color='red', ecolor="lightblue",
                alpha=0.2, label="Residuals")
    ax.axhline(0, color="red", linestyle="--", linewidth=1)  # Zero residual line

    print(reduced_chi_squared)
    print(p_value)
    # Labels and grid
    plt.xlabel("Time (seconds)")
    plt.ylabel("Residuals (meters)")
    plt.legend()
    plt.grid(True)
    plt.savefig("Graphs/empty_residuals.png", dpi=200)
    plt.show()
