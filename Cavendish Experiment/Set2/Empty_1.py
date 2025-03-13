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


def damped_cosine(t, position_0, tau, T, phase, base):
    omega = 2 * np.pi / T
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
    T_guess = (max(t_numeric) / 10)  # Estimate period from observed data
    phase_guess = 0  # Assume initial phase 0
    base_guess = np.mean(position)  # Baseline position

    # Updated initial guesses
    p0 = [amplitude_guess, tau_guess, T_guess, phase_guess, base_guess]

    popt, pcov = curve_fit(damped_cosine, xdata=t_numeric, ydata=position, sigma=position_uncert, p0=p0, maxfev=100000)
    position_fit = np.array(damped_cosine(t_numeric, *popt))
    ax.plot(t_numeric, position_fit, color="black", linestyle="--", label="Fit Line")

    """From Curve fit"""

    print("Position Amplitude (p0) = ", popt[0], "err=", np.sqrt(pcov[0][0]))
    print("Decay Factor (tau) = ", popt[1], "err=", np.sqrt(pcov[1][1]))
    print("Period (T) = ", popt[2], "err=", np.sqrt(pcov[2][2]))
    print("Wave Phase (phi) = ", popt[3], "err=", np.sqrt(pcov[3][3]))
    print("Equilibrium position for Empty (base)= ", popt[4], "err=", np.sqrt((np.sqrt(pcov[4][4]) ** 2) +
                                                                              ((0.00539 / 2) ** 2)))
    plt.xlabel("Time (seconds)")
    plt.ylabel("Laser Position (meters)")
    plt.grid(True)
    plt.legend()
    plt.savefig("Graphs/empty_fit.png", dpi=200)
    plt.show()



    """Convert to theta vs position"""

    print("For theta vs position")

    L = 4.47291  # meters
    L_err = 0.001

    equilibrium = popt[4]
    equilibrium_err = np.sqrt(np.sqrt(pcov[4][4]) ** 2 + (0.00539/2) ** 2)

    thetas = []
    thetas_uncert = []
    for i in range(len(position)):  # Fixed iteration over indices
        z = (position[i] - equilibrium) / (2 * L)
        thetas.append(np.arctan(z))
        # Compute uncertainty
        uncert_theta = (1 / (1 + z ** 2)) * np.sqrt(
            (position_uncert[i] / (2 * L)) ** 2 +
            (equilibrium_err / (2 * L)) ** 2 +
            ((position[i] - equilibrium) * L_err / (2 * L ** 2)) ** 2
        )
        thetas_uncert.append(np.arctan(uncert_theta))
    amplitude_guess = max(thetas) - min(thetas)  # Approximate amplitude
    tau_guess = max(t_numeric) / 5  # Rough estimate of decay time
    T_guess = (max(t_numeric) / 10)  # Estimate period from observed data
    phase_guess = 0  # Assume initial phase 0
    base_guess = np.mean(thetas)  # Baseline position

    # Updated initial guesses
    p0 = [amplitude_guess, tau_guess, T_guess, phase_guess, base_guess]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.errorbar(t_numeric, thetas, yerr=thetas_uncert, marker=".", markersize=5, color='red', ecolor="lightblue",
                alpha=0.2, linestyle="None", label="Mirror Angle")
    popt, pcov = curve_fit(damped_cosine, xdata=t_numeric, p0=p0, ydata=thetas, maxfev=100000000)
    thetas_fit = np.array(damped_cosine(t_numeric, *popt))
    ax.plot(t_numeric, thetas_fit, color="red", label="Fit Line")

    """From Curve fit"""
    print("Theta init (p0) = ", popt[0], "err=", np.sqrt(pcov[0][0]))
    print("Decay Factor (tau) = ", popt[1], "err=", np.sqrt(pcov[1][1]))
    print("Period (T) = ", popt[2], "err=", np.sqrt(pcov[2][2]))
    print("Wave Phase (phi) = ", popt[3], "err=", np.sqrt(pcov[3][3]))
    print("Equilibrium position for Empty (base)= ", popt[4], "err=", np.sqrt((np.sqrt(pcov[4][4]) ** 2)))

    plt.xlabel("Time (seconds)")
    plt.ylabel("Mirror Angle (radians)")
    plt.grid(True)
    plt.legend()
    plt.savefig("Graphs/empty_angle_vs_time.png", dpi=200)
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
                alpha=0.2, linestyle="None", label="Residuals")
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
