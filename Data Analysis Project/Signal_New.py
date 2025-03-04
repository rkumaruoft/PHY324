# Re-import necessary libraries after execution state reset
import pickle
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2

# Set font for plots
font = {'family': 'DejaVu Sans', 'weight': 'normal', 'size': 12}
rc('font', **font)

# Define Gaussian function
def gaussian(x, A, mean, sigma, base):
    return A * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2)) + base

# Define exponential decay function
def decay_function(t, tau, A, B):
    return (A * np.exp(-t / tau)) + B

# Load signal data
file_path = "signal.pkl"
with open(file_path, "rb") as file:
    signal_data = pickle.load(file)

# Noise range and pulse fitting
noise_range = (-2.22420833288144, 2.2261688263127684)
pulse_fit = np.zeros(1000)

# Define pulse shape function
def pulse_shape(t_rise, t_fall):
    xx = np.linspace(0, 4095, 4096)
    yy = -(np.exp(-(xx - 1000) / t_rise) - np.exp(-(xx - 1000) / t_fall))
    yy[:1000] = 0
    yy /= np.max(yy)
    return yy

# Function for fitting pulses
def fit_pulse(x, A):
    _pulse_template = pulse_shape(20, 80)
    xx = np.linspace(0, 4095, 4096)
    return A * np.interp(x, xx, _pulse_template)

# Process signal data
for ievt in range(1000):
    current_data = signal_data['evt_%i' % ievt]
    baseline_avg = np.mean(current_data[0:1000])
    data_cleaned = np.array([x - baseline_avg for x in current_data])
    popt, _ = curve_fit(fit_pulse, range(len(data_cleaned)), data_cleaned)
    pulse_fit[ievt] = popt[0]

# Convert to keV
pulse_fit *= 1000  # Convert from V to mV
c_factor = 39.03182106539658
pulse_fit *= c_factor  # Convert to keV

# Apply noise threshold
pulse_fit = pulse_fit[(pulse_fit < noise_range[0]) | (pulse_fit > noise_range[1])]

# Histogram parameters
num_bins1 = 60
bin_range1 = (min(pulse_fit), max(pulse_fit))

# Generate histogram
n1, bin_edges1, _ = plt.hist(pulse_fit, bins=num_bins1, range=bin_range1, color='k', histtype='step', label='Data')

# Compute bin centers
bin_centers1 = 0.5 * (bin_edges1[1:] + bin_edges1[:-1])

# Compute uncertainties
sig1 = np.sqrt(n1)
sig1 = np.where(sig1 == 0, 1, sig1)  # Replace zeros with 1 for error bars

# Fit exponential decay function to the new histogram
initial_guess = (10, max(n1), min(n1))
popt, pcov = curve_fit(decay_function, bin_centers1, n1, p0=initial_guess)

# Compute residuals
residuals = n1 - decay_function(bin_centers1, *popt)

# Compute chi-squared statistic
chi_squared = np.sum((residuals / sig1) ** 2)
dof = len(n1) - len(popt)
reduced_chi_squared = chi_squared / dof

# Generate fit curve
fit_x = np.linspace(min(bin_centers1), max(bin_centers1), 300)
fit_y = decay_function(fit_x, *popt)
plt.plot(fit_x, fit_y, 'r--', label="Decay Curve")
# Plot error bars
plt.errorbar(bin_centers1, n1, yerr=sig1, fmt='none', c='k')

# Labels and legend
plt.xlabel('Particle Energy (keV)')
plt.ylabel('Number of Events')
plt.xlim(bin_range1)
plt.legend()
plt.savefig("new_plots/Signal_full_hist", dpi=200)
plt.show()

# Fit Gaussian to residuals
initial_guess_gauss = (8, np.mean(residuals), np.std(residuals), min(residuals))
popt_gauss, pcov_gauss = curve_fit(gaussian, bin_centers1, residuals, p0=initial_guess_gauss)

# Generate Gaussian fit curve
fit_x_gauss = np.linspace(min(bin_centers1), max(bin_centers1), 300)
fit_y_gauss = gaussian(fit_x_gauss, *popt_gauss)

# Plot the residuals with Gaussian fit
plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.errorbar(bin_centers1, residuals, yerr=sig1, fmt='.k', label='Residuals')
plt.plot(fit_x_gauss, fit_y_gauss, 'r--', label="Gaussian Fit Curve")
plt.xlabel('Particle Energy (keV)')
plt.ylabel('Residuals')
plt.xlim(bin_range1)
plt.legend()
plt.savefig("new_plots/Signal_full_hist_residuals", dpi=200)
plt.show()

# Compute residuals for Gaussian fit
gaussian_residuals = residuals - gaussian(bin_centers1, *popt_gauss)

# Plot residuals of Gaussian fit
plt.figure(figsize=(8, 5))
plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.errorbar(bin_centers1, gaussian_residuals, yerr=sig1, fmt='.k', label='Residuals')
plt.xlabel('Particle Energy (keV)')
plt.ylabel('Residuals')
plt.xlim(bin_range1)
plt.legend()
plt.savefig("new_plots/Signal_gauss_residuals", dpi=200)
plt.show()
