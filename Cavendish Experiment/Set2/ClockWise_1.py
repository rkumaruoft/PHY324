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


def propagate_uncertainty_G(M1, M2, x1, x2, T, theta, d,
                            sigma_M1, sigma_M2, sigma_x1, sigma_x2, sigma_T, sigma_theta):
    """
    Computes the propagated uncertainty in G based on uncertainties in the measured quantities.

    Parameters:
    M1, M2: Masses of the large spheres
    x1, x2: Distances between the large and small masses
    T: Oscillation period
    theta: Angular displacement (in radians)
    d: Distance from the center of the rod to the small masses
    sigma_M1, sigma_M2: Uncertainties in M1 and M2
    sigma_x1, sigma_x2: Uncertainties in x1 and x2
    sigma_T: Uncertainty in T
    sigma_theta: Uncertainty in theta

    Returns:
    sigma_G: Propagated uncertainty in G
    """

    # Compute the denominator term (which appears in multiple places)
    denominator = T ** 3 * np.cos(theta) ** 2 * (M1 * x2 ** 2 + M2 * x1 ** 2) ** 2

    # First term: Contribution from uncertainty in theta
    term_theta = (T ** 2 * sigma_theta ** 2 * x1 ** 2 * x2 ** 2 * (M1 * x2 ** 2 + M2 * x1 ** 2) ** 2 * (
            theta * np.sin(theta) + np.cos(theta)) ** 2)

    # Second term: Contribution from uncertainties in M1, M2, x1, x2
    term_mass_distance = (T ** 2 * theta ** 2 * np.cos(theta) ** 2 *
                          (4 * M1 ** 2 * sigma_x1 ** 2 * x2 ** 6 +
                           4 * M2 ** 2 * sigma_x2 ** 2 * x1 ** 6 +
                           sigma_M1 ** 2 * x1 ** 2 * x2 ** 6 +
                           sigma_M2 ** 2 * x1 ** 6 * x2 ** 2))

    # Third term: Contribution from uncertainty in T
    term_T = (4 * sigma_T ** 2 * theta ** 2 * x1 ** 2 * x2 ** 2 * (M1 * x2 ** 2 + M2 * x1 ** 2) ** 2 * np.cos(
        theta) ** 2)

    # Compute propagated uncertainty in G
    sigma_G = (8 * np.pi ** 2 / denominator) * np.sqrt(
        d ** 2 * x1 ** 2 * x2 ** 2 * (term_theta + term_mass_distance + term_T))

    return sigma_G


def damped_cosine(t, position_0, tau, T, phase, base):
    omega = 2 * np.pi / T
    return position_0 * np.exp(-t / tau) * np.cos((omega * t) + phase) + base


if __name__ == "__main__":
    time = []
    position = []
    position_uncert = []
    file = open('../data/data_IN_csv/Device2Clockwise.csv', mode='r')
    # Create a CSV reader object
    csv_reader = csv.reader(file)
    # Skip the header row (if there is one)
    next(csv_reader, None)
    time_cutoff = datetime.strptime("12:27:00.000", '%H:%M:%S.%f')
    for row in csv_reader:
        this_time = datetime.strptime(row[0], '%H:%M:%S.%f')
        if float(row[1]) > 0.01:
            time.append(this_time)
            position.append(float(row[1]))
            position_uncert.append(0.00539 / 2)

    t_numeric = np.array([(t - time[0]).total_seconds() for t in time])
    fig, ax = plt.subplots(figsize=(8, 4))
    # ax.scatter(time, position, marker=".", s=5, label="Raw Data")
    ax.errorbar(t_numeric, position, yerr=position_uncert, marker=".", markersize=5, color='yellow', ecolor="lightblue",
                alpha=0.1, linestyle="None", label="Raw Data")

    # """Fitting a damped cosine"""

    # Estimate reasonable initial parameters from data
    amplitude_guess = max(position) - min(position)  # Approximate amplitude
    tau_guess = max(t_numeric) / 5  # Rough estimate of decay time
    T_guess = (max(t_numeric) / 100)  # Estimate period from observed data
    phase_guess = 0  # Assume initial phase 0
    base_guess = np.mean(position)  # Baseline position

    # Updated initial guesses
    p0 = [amplitude_guess, tau_guess, T_guess, phase_guess, base_guess]

    popt, pcov = curve_fit(damped_cosine, xdata=t_numeric, ydata=position, sigma=position_uncert, p0=p0, maxfev=100000)
    position_fit = np.array(damped_cosine(t_numeric, *popt))
    plt.plot(t_numeric, position_fit, color="red", label="Fit Line")

    """From Curve fit"""
    print("Position Amplitude (p0) = ", popt[0], "err=", np.sqrt(pcov[0][0]))
    print("Decay Factor (tau) = ", popt[1], "err=", np.sqrt(pcov[1][1]))
    print("Period (T) = ", popt[2], "err=", np.sqrt(pcov[2][2]) ** 2 + ())
    print("Wave Phase (phi) = ", popt[3], "err=", np.sqrt(pcov[3][3]))
    print("Equilibrium position for Clock (base)= ", popt[4], "err=", np.sqrt((np.sqrt(pcov[4][4]) ** 2) +
                                                                              ((0.00539 / 2) ** 2)))

    plt.xlabel("Time (seconds)")
    plt.ylabel("Laser Position (meters)")
    plt.grid(True)
    plt.legend()
    plt.savefig("Graphs/clock_fit.png", dpi=200)
    plt.show()

    """Convert to theta vs position"""

    print("For theta vs position")

    L = 4.47291  # meters
    L_err = 0.001

    equilibrium = popt[4]
    equilibrium_err = np.sqrt(np.sqrt(pcov[4][4]) ** 2 + (0.00539 / 2) ** 2)

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
        thetas_uncert.append(uncert_theta)
    amplitude_guess = max(thetas) - min(thetas)  # Approximate amplitude
    tau_guess = max(t_numeric) / 5  # Rough estimate of decay time
    T_guess = (max(t_numeric) / 10)  # Estimate period from observed data
    phase_guess = 0  # Assume initial phase 0
    base_guess = np.mean(thetas)  # Baseline position

    # Updated initial guesses
    p0 = [amplitude_guess, tau_guess, T_guess, phase_guess, base_guess]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.errorbar(t_numeric, thetas, yerr=thetas_uncert, marker=".", markersize=5, color='yellow', ecolor="lightblue",
                alpha=0.1, linestyle="None", label="Raw Data")
    popt, pcov = curve_fit(damped_cosine, xdata=t_numeric, p0=p0, ydata=thetas, maxfev=10000)
    thetas_fit = np.array(damped_cosine(t_numeric, *popt))
    ax.plot(t_numeric, thetas_fit, color="red", label="Fit Line")

    """From Curve fit"""
    print("Theta init (p0) = ", popt[0], "err=", np.sqrt(pcov[0][0]))
    print("Decay Factor (tau) = ", popt[1], "err=", np.sqrt(pcov[1][1]))
    print("Period (T) = ", popt[2], "err=", np.sqrt(pcov[2][2]))
    print("Wave Phase (phi) = ", popt[3], "err=", np.sqrt(pcov[3][3]))
    print("Equilibrium angle for Clock (base)= ", popt[4], "err=", np.sqrt((np.sqrt(pcov[4][4]) ** 2)))

    plt.xlabel("Time (seconds)")
    plt.ylabel("Mirror Angle (radians)")
    plt.grid(True)
    plt.legend()
    plt.savefig("Graphs/clock_angle_vs_time.png", dpi=200)
    plt.show()

    """Residuals"""
    # Compute residuals
    residuals = thetas - thetas_fit

    # Compute reduced chi-squared
    degrees_of_freedom = len(thetas) - len(popt)
    chi_squared = np.sum((residuals / thetas_uncert) ** 2)
    reduced_chi_squared = chi_squared / degrees_of_freedom

    # Compute fit probability (p-value)
    p_value = 1 - chi2.cdf(chi_squared, degrees_of_freedom)

    # Plot residuals
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.errorbar(t_numeric, residuals, yerr=thetas_uncert, marker=".", markersize=5, color='red', ecolor="lightblue",
                alpha=0.2, linestyle="None", label="Residuals")
    ax.axhline(0, color="red", linestyle="--", linewidth=1)  # Zero residual line

    print(reduced_chi_squared)
    print(p_value)
    # Labels and grid
    plt.xlabel("Time (seconds)")
    plt.ylabel("Residuals (radians)")
    plt.legend()
    plt.grid(True)
    plt.savefig("Graphs/clock_residuals.png", dpi=200)
    plt.show()

    """Calculate G"""
    print("Calculate G")

    # Equilibrium position for Empty(base)=  0.060560769666643435 err= 0.002695026884311751
    # Equilibrium angle for Empty(base)=  4.900580537830072e-10 err= 1.3456725614433343e-06
    # Equilibrium position for Clock(base)=  0.07664600062927877 err= 0.0026950326758820387

    theta_equil_empty = 4.900580537830072e-10
    theta_equil_empty_err = 0.002695026884311751
    pos_equil_empty = 0.060560769666643435
    pos_equil_empty_err = 1.3456725614433343e-06

    theta_equil_clock = popt[4]
    theta_equil_clock_err = np.sqrt((np.sqrt(pcov[4][4]) ** 2))
    pos_equil_clock = 0.07664600062368222
    pos_equil_clock_err = 0.0026950326758820387

    """Ball Specs"""
    w = 0.02938

    d = 0.05

    m1, m1_err = 1477.4e-3, 0.1e-3
    m2, m2_err = 1472.9e-3, 0.1e-3
    r1, r1_err = 63.72e-3, 0.01e-3
    r2, r2_err = 63.89e-3, 0.01e-3

    L = 4.47291  # meters
    L_err = 0.001

    # Period from fit
    T = 305.3897876812016  # seconds
    T_err = np.sqrt(0.09223622916845486 ** 2 + (1 / np.sqrt(12)) ** 2)  # seconds

    print("T", T, T_err)
    # **Method 1: Using theta equilibrium difference**
    theta_method1 = abs(theta_equil_clock - theta_equil_empty)
    theta_err_method1 = theta_equil_clock_err

    # **Method 2: Using position equilibrium difference**
    s = abs(pos_equil_clock - pos_equil_empty)
    s_err = max(pos_equil_clock_err, pos_equil_empty_err)
    print("S", s, s_err)

    theta_method2 = np.arctan(s / (2 * L)) / 2
    theta_err_method2 = theta_method2 * np.sqrt(((s_err / s) ** 2) + ((L_err / L) ** 2))
    print("Theta", theta_method2, theta_err_method2)

    # Compute x1 and x2 for both methods
    x1_method2 = r1 + (w / 2) - (d * np.sin(theta_method2))
    x1_err_method2 = np.sqrt((r1_err ** 2) + ((d * np.cos(theta_method2) * theta_err_method2) ** 2))

    x2_method2 = r2 + (w / 2) - (d * np.sin(theta_method2))
    x2_err_method2 = np.sqrt((r2_err ** 2) + ((d * np.cos(theta_method2) * theta_err_method2) ** 2))

    # Compute G using Method 2
    G_method2 = (8 * (np.pi ** 2) * d * theta_method2) / (
            np.cos(theta_method2) * (T ** 2) * ((m1 / (x1_method2 ** 2)) + (m2 / (x2_method2 ** 2)))
    )

    G_method2_err = propagate_uncertainty_G(m1, m2, x1_method2, x2_method2, T, theta_method2, d,
                                            m1_err, m2_err, x1_err_method2, x2_err_method2, T_err, theta_err_method2)

    # Display results for G using Method 2
    print(G_method2, G_method2_err)
