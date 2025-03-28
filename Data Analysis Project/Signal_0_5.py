import pickle
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2
from scipy.stats import iqr

font = {'family': 'DejaVu Sans',
        'weight': 'normal',
        'size': 10}
rc('font', **font)


# This changes the fonts for all graphs to make them bigger.


def myGauss(x, A, mean, width, base):
    return A * np.exp(-(x - mean) ** 2 / (2 * width ** 2)) + base


# This is my fitting function, a Guassian with a uniform background.

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


# fit_pulse can be used by curve_fit to fit a pulse to the pulse_shape

with open("signal.pkl", "rb") as file:
    signal_data = pickle.load(file)


# for itrace in range(500,600):
#     plt.plot(signal_data['evt_%i' % itrace], alpha=0.3)
# plt.xlabel('Sample Index')
# plt.ylabel('Readout (V)')
# plt.title('Signal data (10 sets)')
# plt.legend(loc=1)
# plt.show()


"""
This shows the first 10 data sets on top of each other.
Always a good idea to look at some of your data before analysing it!
It also plots our pulse template which has been scaled to be slightly
larger than any of the actual pulses to make it visible.
"""
noise_range =  (np.float64(-1.482478806349072), np.float64(1.4844392997804006))
pulse_fit = np.zeros(1000)


for ievt in range(1000):
    current_data = signal_data['evt_%i' % ievt]
    baseline_avg = np.mean(current_data[0:1000])
    data_cleaned = [x - baseline_avg for x in current_data]
    popt, pcov = curve_fit(fit_pulse, range(len(data_cleaned)), data_cleaned)
    pulse_fit[ievt] = popt[0]


pulse_fit *= 1000  # convert from V to mV
c_factor = 39.03182106539658
pulse_fit *= c_factor # convert to keV

pulse_fit = pulse_fit[(pulse_fit < noise_range[0]) | (pulse_fit > noise_range[1])]

num_bins1 = 60
bin_range1 = (min(pulse_fit), max(pulse_fit))
"""
These two values were picked by trial and error. You'll
likely want different values for each estimator.
"""

n1, bin_edges1, _ = plt.hist(pulse_fit, bins=num_bins1, range=bin_range1, color='k', histtype='step', label='Data')
# This plots the histogram AND saves the counts and bin_edges for later use


plt.xlabel('Particle Energy (keV)')
plt.ylabel('Number of Events')
plt.xlim(bin_range1)

bin_centers1 = 0.5 * (bin_edges1[1:] + bin_edges1[:-1])
# If the legend covers some data, increase the plt.xlim value, maybe (0,0.5)

sig1 = np.sqrt(n1)
sig1 = np.where(sig1 == 0, 1, sig1)
# The uncertainty on 0 count is 1, not 0. Replace all 0s with 1s.

plt.errorbar(bin_centers1, n1, yerr=sig1, fmt='none', c='k')
# This adds errorbars to the histograms, where each uncertainty is sqrt(y)
plt.show()



"""
Particle b/w 1 to 5
"""
sectioned_signal = pulse_fit[(pulse_fit >= noise_range[1]) & (pulse_fit <= 5.5)]

num_bins_sectioned = 80
bin_range_sectioned = (1.3, max(sectioned_signal))  # Directly specify the range

n_sectioned, bin_edges_sectioned, _ = plt.hist(
    sectioned_signal, bins=num_bins_sectioned, range=bin_range_sectioned, color='b', histtype='step', label='Sectioned Data'
)

plt.xlabel('Particle Energy (keV)')
plt.ylabel('Number of Events')
plt.xlim(bin_range_sectioned)

# Calculate bin centers and uncertainties
bin_centers_sectioned = 0.5 * (bin_edges_sectioned[1:] + bin_edges_sectioned[:-1])
sig_sectioned = np.sqrt(n_sectioned)
sig_sectioned = np.where(sig_sectioned == 0, 1, sig_sectioned)

# Add error bars
plt.errorbar(bin_centers_sectioned, n_sectioned, yerr=sig_sectioned, fmt='none', c='k')
# plt.title('Sectioned Signal Data (0 keV to 5 keV)')

popt, pcov = curve_fit(myGauss, bin_centers_sectioned, n_sectioned, sigma=sig_sectioned,
                       p0=[12, 2.0, 1, min(n_sectioned)], absolute_sigma=True)

# Calculate the fitted Gaussian curve
x_bestfit = np.linspace(bin_range_sectioned[0], bin_range_sectioned[1], 1000)
y_bestfit = myGauss(x_bestfit, *popt)

plt.plot(x_bestfit, y_bestfit, label="Fit Line", color='r')

# Annotate the Gaussian parameters
mean = popt[1]
sigma = popt[2]
chisquared = np.sum(((n_sectioned - myGauss(bin_centers_sectioned, *popt)) / sig_sectioned) ** 2)
dof = len(bin_centers_sectioned) - len(popt)  # Degrees of freedom
reduced_chisquared = chisquared / dof
print("chi_red= ", reduced_chisquared)
print("Amp = ", popt[0], "err= ", np.sqrt(pcov[0][0]))
print("Mean (mu):", mean, "keV" , "err = ", np.sqrt(pcov[1][1]))
print("Standard Deviation (sigma):", sigma, "keV", np.sqrt(pcov[2][2]))

dof = num_bins_sectioned - len(popt)
chi_prob = 1 - chi2.cdf(chisquared, dof)
print("chi_prob= ", chi_prob)
plt.legend()
plt.show()
