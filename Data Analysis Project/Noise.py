import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2


# Define Gaussian function for fitting
def myGauss(x, A, mean, width, base):
    return A * np.exp(-(x - mean) ** 2 / (2 * width ** 2)) + base


# Load the noise data
with open("noise.pkl", "rb") as file:
    noise_data = pickle.load(file)

# Gather all raw noise samples into a single array
all_noise_samples = []
for ievt in range(1000):
    current_data = noise_data[f'evt_{ievt}']
    all_noise_samples.extend(current_data)

# Convert to mV (raw values)
all_noise_samples = np.array(all_noise_samples) * 1000  # Convert from V to mV

# Apply calibration factor to convert mV to keV
calibration_factor = 39.03182106539658  # keV/mV
all_noise_samples_keV = all_noise_samples * calibration_factor

# Create histogram of raw noise fluctuations
num_bins = 100
bin_range = (min(all_noise_samples_keV), max(all_noise_samples_keV))
print(bin_range)
n, bin_edges, _ = plt.hist(all_noise_samples_keV, bins=num_bins, range=bin_range, color='k', histtype='step',
                           label='Data')
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# Add error bars
sig = np.sqrt(n)
sig = np.where(sig == 0, 1, sig)  # Replace 0 uncertainties with 1
plt.errorbar(bin_centers, n, yerr=sig, fmt='none', c='k')

# Fit Gaussian to the histogram
popt, pcov = curve_fit(myGauss, bin_centers, n, sigma=sig,
                       p0=[max(n), np.mean(all_noise_samples_keV), np.std(all_noise_samples_keV), min(n)],
                       absolute_sigma=True)

# Plot Gaussian fit
x_bestfit = np.linspace(bin_range[0], bin_range[1], 1000)
y_bestfit = myGauss(x_bestfit, *popt)
plt.plot(x_bestfit, y_bestfit, label='Gaussian Fit', color='r')

# Plot details
plt.xlabel('Raw Noise Energy (keV)')
plt.ylabel('Number of Events')
plt.xlim(bin_range)
plt.legend()
# plt.title('Histogram of Raw Noise Fluctuations with Gaussian Fit')

# Print results
mean = popt[1]
sigma = popt[2]
noise_range = (mean - (2 * sigma), mean + (2 * sigma))
print("Histogram of Raw Noise Fluctuations (in keV):")
print("Amp = ", popt[0], "err= ", np.sqrt(pcov[0][0]))
print("Mean (mu):", mean, "keV", "err = ", np.sqrt(pcov[1][1]))
print("Standard Deviation (sigma):", sigma, "keV", np.sqrt(pcov[2][2]))
print("Noise Range (mu ± 3σ):", noise_range)

# Calculate chi-squared
n1_fit = myGauss(bin_centers, *popt)  # Gaussian fit values
chisquared = np.sum(((n - n1_fit) / sig) ** 2)  # Chi-squared
dof = num_bins - len(popt)  # Degrees of freedom
reduced_chisquared = chisquared / dof  # Reduced chi-squared

# Print results
print("Reduced Chi-Squared (χ²/DOF):", reduced_chisquared)
chi_prob = 1 - chi2.cdf(chisquared, dof)
print("chi_prob= ", chi_prob)

# Save and show the plot
plt.savefig("plots/histogram_raw_noise_fluctuations_keV.png", dpi=200)
plt.show()

"""Residuals"""
# Compute residuals for the Gaussian fit after calibration
residuals = n - myGauss(bin_centers, *popt)

# Plot the residuals with black markers and error bars using '.' marker format
plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.errorbar(bin_centers, residuals, yerr=sig, fmt='.k', label='Residuals')
plt.xlabel('Particle Energy (keV)')
plt.ylabel('Residuals')
plt.legend()
plt.savefig("new_plots/noise_residuals", dpi=200)
plt.show()
