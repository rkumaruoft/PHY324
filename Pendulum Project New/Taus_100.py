# Fast CSV reading using Pandas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2


def abs_sin(x, a, b, c):
    return a * np.abs(np.sin(x + b)) + c


data = pd.read_csv('periods_and_taus_100.csv', skiprows=1, header=None).values  # Skipping header row

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


mask = tau >= 0
tau = tau[mask]
tau_err = tau_err[mask]
init_theta = init_theta[mask]

plt.errorbar(init_theta, tau, yerr=tau_err, linestyle="none", marker="o", color='b', label="Data")
plt.show()

# Perform quadratic fit with uncertainties
# sin_popt, sin_pcov = curve_fit(abs_sin, init_theta, tau, sigma=tau_err, absolute_sigma=True, maxfev=10000)


