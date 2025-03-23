# Fast CSV reading using Pandas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2


def quadratic_fit(x, a, b, c):
    return a * x ** 3 + b * x**2 + c


def absolute_fit(x, a, b):
    return a * np.abs(x) + b


data = pd.read_csv('periods_and_taus_50.csv', skiprows=1, header=None).values  # Skipping header row

# Extract columns efficiently
init_theta, period, p_err, tau, tau_err = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]
p_err = [np.sqrt(x ** 2 + 0.033322 ** 2) for x in p_err]
# mask = tau >= 0
# tau = tau[mask]
# tau_err = tau_err[mask]

init_theta = np.array(init_theta)
period = np.array(period)
p_err = np.array(p_err)

sorted_ind = np.argsort(init_theta)
init_theta = init_theta[sorted_ind]
period = period[sorted_ind]
p_err = p_err[sorted_ind]
tau = tau[sorted_ind]
tau_err = tau_err[sorted_ind]

# Perform quadratic fit with uncertainties
quad_popt, quad_pcov = curve_fit(quadratic_fit, init_theta, period, sigma=p_err, absolute_sigma=True, maxfev=10000)
abs_popt, abs_pcov = curve_fit(absolute_fit, init_theta, period, sigma=p_err, absolute_sigma=True, maxfev=10000)

# Extract fit parameters and uncertainties
a_quad, b_quad, c_quad = quad_popt
a_quad_err, b_quad_err, c_quad_err = np.sqrt(np.diag(quad_pcov))

a_abs, b_abs = abs_popt
a_abs_err, b_abs_err = np.sqrt(np.diag(abs_pcov))

# Generate fit data
period_fit_quad = quadratic_fit(init_theta, *quad_popt)

# Compute reduced chi-squared for quadratic fit
residuals_quad = period - quadratic_fit(init_theta, *quad_popt)
chi_squared_quad = np.sum((residuals_quad / p_err) ** 2)
chi_squared_reduced_quad = chi_squared_quad / (len(init_theta) - len(quad_popt))


period_fit_abs = absolute_fit(init_theta, *abs_popt)
# Compute reduced chi-squared for absolute fit
residuals_abs = period - absolute_fit(init_theta, *abs_popt)
chi_squared_abs = np.sum((residuals_abs / p_err) ** 2)
chi_squared_reduced_abs = chi_squared_abs / (len(init_theta) - len(abs_popt))
# Compute fit probability (p-value) from reduced chi-squared
dof_abs = len(init_theta) - len(abs_popt)  # Degrees of freedom
fit_probability_abs = 1 - chi2.cdf(chi_squared_abs, dof_abs)  # p-value

# Compute fit probability (p-value) from reduced chi-squared
dof = len(init_theta) - len(quad_popt)  # Degrees of freedom
fit_probability = 1 - chi2.cdf(chi_squared_quad, dof)  # p-value

plt.errorbar(init_theta, period, yerr=p_err, linestyle="none", marker="o", color='b', label="Data")
plt.plot(init_theta, period_fit_quad, 'r-', label="Quadratic Fit")
# plt.plot(init_theta, period_fit_abs, 'y-', label="Absolute Fit")

# Labels
plt.xlabel('Start Angle (Degrees)')
plt.ylabel('Period (Seconds)')

# Simplified legend
plt.legend()
print(quad_popt)
# Display the plot
plt.grid(True)
plt.savefig("PeriodsVsTheta.png", dpi=200)
print(f"a (Quadratic Term) = {a_quad:.8f} ± {a_quad_err:.8f}")
# print(f"b (Linear Term) = {b_quad:.8f} ± {b_quad_err:.8f}")
print(f"c (Intercept) = {c_quad:.8f} ± {c_quad_err:.8f}")
print(f"Reduced Chi-Squared (χ²_r) = {chi_squared_reduced_quad:.4f}")
print(f"Fit Probability (p-value) = {fit_probability:.6f}")
# print("===================================================")
# print(f"a (Absolute Term) = {a_abs:.8f} ± {a_abs_err:.8f}")
# # print(f"b (Linear Term) = {b_quad:.8f} ± {b_quad_err:.8f}")
# print(f"c (Intercept) = {b_abs:.8f} ± {b_abs_err:.8f}")
# print(f"Reduced Chi-Squared (χ²_r) = {chi_squared_reduced_abs:.4f}")
# print(f"Fit Probability (p-value) = {fit_probability_abs:.6f}")
plt.show()

"""Quad residuals"""

# Plot the residuals
plt.figure(figsize=(10, 6))
plt.errorbar(init_theta, residuals_quad, yerr=p_err, label='Residuals', color='r', marker=".", linestyle="None")
plt.axhline(y=0, color='k', linestyle='--')  # Line at 0 for reference
plt.title('Residuals for Quadratic Fit')
plt.xlabel('Start Angle (Degrees)')
plt.ylabel('Residuals (Period - Fit)')
plt.grid(True)
plt.legend()
plt.show()
