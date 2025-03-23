# Fast CSV reading using Pandas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2


def continuous_piecewise_poly_rising(x_data, a, b, c):
    # Apply Region 1 for x <= 0 and Region 2 for x > 0
    return np.piecewise(x_data, [x_data <= 0, x_data > 0],
                        [lambda x: a * x**2 -b * x + c,
                         lambda x: a * x**2 + b * x + c])


def reduced_chi_squared(residuals, errors, n_params):
    return np.sum((residuals / errors) ** 2) / (len(residuals) - n_params)


# Load the data
data = pd.read_csv('periods_and_taus_200.csv', skiprows=1, header=None).values  # Skipping header row

# Extract columns efficiently
init_theta, period, p_err, tau, tau_err = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]
p_err = [np.sqrt(x ** 2 + 0.033322 ** 2) for x in p_err]

init_theta = np.array(init_theta)
period = np.array(period)
p_err = np.array(p_err)

sorted_ind = np.argsort(init_theta)
init_theta = init_theta[sorted_ind]
period = period[sorted_ind]
p_err = p_err[sorted_ind]
tau = tau[sorted_ind]
tau_err = tau_err[sorted_ind]

# Apply the mask to limit tau to values below 175
mask = (0 <= tau) & (tau < 175)
tau = tau[mask]
tau_err = tau_err[mask]
init_theta = init_theta[mask]

# Use curve fitting to find the optimal parameters
params, covariance = curve_fit(continuous_piecewise_poly_rising, init_theta, tau,
                                sigma=tau_err, absolute_sigma=True,
                                p0=[0.0426, 2.255, 95.1024])


# Generate the fitted curve based on the optimized parameters
x_for_fit = np.arange(np.min(init_theta), np.max(init_theta), 0.1)
y_fit = continuous_piecewise_poly_rising(x_for_fit, *params)

# Plot the data and the fitted curve
plt.errorbar(init_theta, tau, yerr=tau_err, linestyle="none", marker="o", color='b', label="Data")
plt.plot(x_for_fit, y_fit, label='Piecewise Polynomial Fit', color='r')

plt.xlabel('Initial Angle (degrees)')
plt.ylabel('Tau (seconds)')
plt.grid(True)
plt.legend()
plt.show()

# Calculate the residuals for the fit
residuals = tau - continuous_piecewise_poly_rising(init_theta, *params)

# Plot the residuals
plt.errorbar(init_theta, residuals, yerr=tau_err, linestyle="none", marker="o", color='g', label="Residuals")
plt.axhline(0, color='black', linestyle='--', label="Zero Line")


# Calculate the residuals
residuals = tau - continuous_piecewise_poly_rising(init_theta, *params)

# Number of parameters for the piecewise polynomial fit (3 parameters: a, b, c)
n_params = 3

# Calculate the reduced chi-squared
chi_squared_red = reduced_chi_squared(residuals, tau_err, n_params)

print(chi_squared_red)

plt.xlabel('Initial Angle (degrees)')
plt.ylabel('Residuals (τ - Fit) (seconds)')
plt.grid(True)
plt.legend()
plt.show()

