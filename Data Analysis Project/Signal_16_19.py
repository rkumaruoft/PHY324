import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Load the signal data from the provided file
with open("signal.pkl", "rb") as file:
    signal_data = pickle.load(file)

# Define the Gaussian function
def myGauss(x, A, mean, width, base):
    return A * np.exp(-(x - mean) ** 2 / (2 * width ** 2)) + base

# Define the pulse shape model
def pulse_shape(t_rise, t_fall):
    xx = np.linspace(0, 4095, 4096)
    yy = -(np.exp(-(xx - 1000) / t_rise) - np.exp(-(xx - 1000) / t_fall))
    yy[:1000] = 0
    yy /= np.max(yy)
    return yy

def fit_pulse(x, A):
    _pulse_template = pulse_shape(20, 80)
    xx = np.linspace(0, 4095, 4096)
    return A * np.interp(x, xx, _pulse_template)

# Filter noise and process data
noise_range = (np.float64(-1.482478806349072), np.float64(1.4844392997804006))
pulse_fit = np.zeros(1000)

for ievt in range(1000):
    current_data = signal_data[f'evt_{ievt}']
    baseline_avg = np.mean(current_data[0:1000])
    data_cleaned = [x - baseline_avg for x in current_data]
    popt, pcov = curve_fit(fit_pulse, range(len(data_cleaned)), data_cleaned)
    pulse_fit[ievt] = popt[0]

# Convert units to keV
pulse_fit *= 1000  # convert from V to mV
c_factor = 39.03182106539658
pulse_fit *= c_factor  # convert to keV

# Filter data for valid ranges
pulse_fit = pulse_fit[(pulse_fit < noise_range[0]) | (pulse_fit > noise_range[1])]

# Section the histogram between 15 and 20 keV
sectioned_signal_15_20 = pulse_fit[(pulse_fit >= 15) & (pulse_fit <= 20)]

# Create a histogram for this range
num_bins_reduced = 20  # Reduced number of bins
bin_range_15_20 = (min(sectioned_signal_15_20), max(sectioned_signal_15_20))
n_15_20_reduced, bin_edges_15_20_reduced, _ = plt.hist(
    sectioned_signal_15_20, bins=num_bins_reduced, range=bin_range_15_20,
    color='b', histtype='step', label='Original Histogram (15-20 keV)'
)

# Calculate bin centers and uncertainties
bin_centers_15_20_reduced = 0.5 * (bin_edges_15_20_reduced[1:] + bin_edges_15_20_reduced[:-1])
sig_15_20_reduced = np.sqrt(n_15_20_reduced)
sig_15_20_reduced = np.where(sig_15_20_reduced == 0, 1, sig_15_20_reduced)

# Extend the histogram range by adding zeros before 15 keV and after 20 keV
extended_bin_centers = np.concatenate((
    [bin_centers_15_20_reduced[0] - 0.5],
    bin_centers_15_20_reduced,
    [bin_centers_15_20_reduced[-1] + 0.5]
))

extended_counts = np.concatenate((
    [0],  # Add zero counts before
    n_15_20_reduced,
    [0]  # Add zero counts after
))

extended_errors = np.concatenate((
    [1],  # Error on the zero counts (set to 1)
    sig_15_20_reduced,
    [1]  # Error on the zero counts (set to 1)
))

# Fit a single Gaussian to the extended data
popt_gauss_extended, pcov_gauss_extended = curve_fit(
    myGauss, extended_bin_centers, extended_counts,
    sigma=extended_errors,
    p0=[max(extended_counts), 17.5, 0.5, 0],  # Initial guesses for A, mean, width, base
    absolute_sigma=True
)

# Calculate the fitted Gaussian curve
x_bestfit_gauss_extended = np.linspace(min(extended_bin_centers), max(extended_bin_centers), 1000)
y_bestfit_gauss_extended = myGauss(x_bestfit_gauss_extended, *popt_gauss_extended)

# Plot the extended histogram with the Gaussian fit
plt.hist(sectioned_signal_15_20, bins=num_bins_reduced, range=bin_range_15_20,
         color='b', histtype='step')
plt.errorbar(extended_bin_centers, extended_counts, yerr=extended_errors, fmt='none', c='k',)
plt.plot(x_bestfit_gauss_extended, y_bestfit_gauss_extended, label='Gaussian Fit', color='r')

# Extract Gaussian parameters
A_gauss_ext, mean_gauss_ext, sigma_gauss_ext, base_gauss_ext = popt_gauss_extended
chisquared_gauss_ext = np.sum(
    ((extended_counts - myGauss(extended_bin_centers, *popt_gauss_extended)) / extended_errors) ** 2
)
dof_gauss_ext = len(extended_counts) - len(popt_gauss_extended)
reduced_chisquared_gauss_ext = chisquared_gauss_ext / dof_gauss_ext

# Annotate the plot with fit parameters
plt.text(min(extended_bin_centers) + 0.5, max(extended_counts) * 0.8, r'A = {:.2f}'.format(A_gauss_ext))
plt.text(min(extended_bin_centers) + 0.5, max(extended_counts) * 0.7, r'$\mu$ = {:.2f} keV'.format(mean_gauss_ext))
plt.text(min(extended_bin_centers) + 0.5, max(extended_counts) * 0.6, r'$\sigma$ = {:.2f} keV'.format(sigma_gauss_ext))
plt.text(min(extended_bin_centers) + 0.5, max(extended_counts) * 0.5, r'$\chi^2$/DOF = {:.2f}'.format(reduced_chisquared_gauss_ext))

plt.xlabel('Particle Energy (keV)')
plt.ylabel('Number of Events')
plt.title('Gaussian Fit to Extended Sectioned Data (15-20 keV)')
plt.legend()
plt.show()