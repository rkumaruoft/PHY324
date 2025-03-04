# Reload data since previous execution state was lost

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chi2


# Define the absolute linear fit function
def fit_for_tau(x, a, b):
    return -a * np.abs(x) + b


def calculate_tau_error(tau_fit, tau_fit_err, t_max, delta_t):
    """
    Calculate the total uncertainty in tau (damping time constant)
    given the curve-fit error and uncertainty in time.

    Parameters:
    - tau_fit: Fitted value of tau
    - tau_fit_err: Error in tau from curve fitting
    - t_max: Total duration of measurement (max time - min time)
    - delta_t: Known uncertainty in time measurement

    Returns:
    - tau_err_total: Total uncertainty in tau
    """
    # Compute total tau error using error propagation formula
    tau_err_total = np.sqrt(tau_fit_err ** 2 + (tau_fit / t_max * delta_t) ** 2)
    return tau_err_total


# Load data from CSV
file_path = "for_taus.csv"
data = pd.read_csv(file_path, skiprows=1, header=None).values  # Skipping header row

# Extract columns
init_theta, taus, t_err, t_max = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
tau_err = []
for i in range(len(taus)):
    tau_err.append(calculate_tau_error(taus[i], t_err[i], t_max[i], 0.017))

# Perform absolute linear fit with uncertainties
abs_popt, abs_pcov = curve_fit(fit_for_tau, init_theta, taus, sigma=tau_err, absolute_sigma=True, maxfev=10000)

# Extract fit parameters and uncertainties
a_abs, b_abs = abs_popt
a_abs_err, b_abs_err = np.sqrt(np.diag(abs_pcov))

# Generate fit data
theta_fit = np.linspace(min(init_theta), max(init_theta), 100)
tau_fit_abs = fit_for_tau(theta_fit, *abs_popt)

# Compute reduced chi-squared for absolute linear fit
residuals_abs = taus - fit_for_tau(init_theta, *abs_popt)
chi_squared_abs = np.sum((residuals_abs / tau_err) ** 2)
chi_squared_reduced_abs = chi_squared_abs / (len(init_theta) - len(abs_popt))

# Compute fit probability
dof_abs = len(init_theta) - len(abs_popt)  # Degrees of freedom
fit_probability_abs = 1 - chi2.cdf(chi_squared_abs, dof_abs)  # p-value

# Plot the absolute linear fit
plt.figure(figsize=(8, 5))
plt.errorbar(init_theta, taus, yerr=tau_err, linestyle="none", marker=".", color='b', label="Data")
plt.plot(theta_fit, tau_fit_abs, 'r-', label="Absolute Linear Fit")

# Labels
plt.xlabel('Start Angle (Degrees)')
plt.ylabel('Decay Factor (Seconds)')

# Simplified legend
plt.legend()

# Display the plot
plt.grid(True)
plt.savefig("TausVsTheta.png", dpi=200)
plt.show()

# Print fit statistics for absolute linear fit
fit_statistics_abs = {
    "a (Slope)": f"{a_abs:.8f} ± {a_abs_err:.8f}",
    "b (Intercept)": f"{b_abs:.8f} ± {b_abs_err:.8f}",
    "Reduced Chi-Squared (χ²_r)": f"{chi_squared_reduced_abs:.4f}",
    "Fit Probability (p-value)": f"{fit_probability_abs:.6f}",
}

print(fit_statistics_abs)
