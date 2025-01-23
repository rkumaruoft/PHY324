import pickle
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2

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


for itrace in range(10):
    plt.plot(signal_data['evt_%i' % itrace], alpha=0.3)
plt.xlabel('Sample Index')
plt.ylabel('Readout (V)')
plt.title('Signal data (10 sets)')
plt.legend(loc=1)
plt.show()


"""
This shows the first 10 data sets on top of each other.
Always a good idea to look at some of your data before analysing it!
It also plots our pulse template which has been scaled to be slightly
larger than any of the actual pulses to make it visible.
"""

amp2 = np.zeros(1000)


for ievt in range(1000):
    current_data = signal_data['evt_%i' % ievt]
    baseline_avg = np.mean(current_data[0:1000])

amp2 *= 1000  # convert from V to mV
c_factor = 35.62321912581909
amp2 *= c_factor
num_bins1 = 80
print(min(amp2), max(amp2))
bin_range1 = (min(amp2), max(amp2))
"""
These two values were picked by trial and error. You'll
likely want different values for each estimator.
"""

n1, bin_edges1, _ = plt.hist(amp2, bins=num_bins1, range=bin_range1, color='k', histtype='step', label='Data')
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
plt.savefig("Signal_Amp2")
plt.show()
