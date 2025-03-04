import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chi2


def quadratic_fit(x, a, c):
    return a * x ** 2 + c


# Fast CSV reading using Pandas
data = pd.read_csv('for_periods.csv', skiprows=1, header=None).values  # Skipping header row

# Extract columns efficiently
init_theta, period, p_err = data[:, 0], data[:, 1], data[:, 2]
p_err = p_err
# Perform quadratic fit with uncertainties
quad_popt, quad_pcov = curve_fit(quadratic_fit, init_theta, period, sigma=p_err, absolute_sigma=True, maxfev=10000)

# Extract fit parameters and uncertainties
a_quad, c_quad = quad_popt
a_quad_err, c_quad_err = np.sqrt(np.diag(quad_pcov))

# Generate fit data
theta_fit = np.linspace(min(init_theta), max(init_theta), 100)
period_fit_quad = quadratic_fit(theta_fit, *quad_popt)

# Compute reduced chi-squared for quadratic fit
residuals_quad = period - quadratic_fit(init_theta, *quad_popt)
chi_squared_quad = np.sum((residuals_quad / p_err) ** 2)
chi_squared_reduced_quad = chi_squared_quad / (len(init_theta) - len(quad_popt))

# Compute fit probability (p-value) from reduced chi-squared
dof = len(init_theta) - len(quad_popt)  # Degrees of freedom
fit_probability = 1 - chi2.cdf(chi_squared_quad, dof)  # p-value

plt.errorbar(init_theta, period, yerr=p_err, linestyle="none", marker="o", color='b', label="Data")
plt.plot(theta_fit, period_fit_quad, 'r-', label="Quadratic Fit")

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
plt.show()
