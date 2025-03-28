import pickle
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2

font = {'family': 'DejaVu Sans',
        'weight': 'normal',
        'size': 12}
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

with open("calibration.pkl", "rb") as file:
    calibration_data = pickle.load(file)

# pulse_template = pulse_shape(20, 80)
# plt.plot(pulse_template / 2000, label='Pulse Template', color='r')
# for itrace in range(10):
#     plt.plot(calibration_data['evt_%i' % itrace], alpha=0.3)
# plt.xlabel('Sample Index')
# plt.ylabel('Readout (V)')
# plt.title('Calibration data (10 sets)')
# plt.legend(loc=1)
# plt.show()

area1 = np.zeros(1000)
# These are the 6 energy estimators as empty arrays of the correct size.

for ievt in range(1000):
    current_data = calibration_data['evt_%i' % ievt]
    area1[ievt] = np.sum(current_data)

area1 *= 1000  # convert from V to mV

num_bins1 = 60
bin_range1 = (min(area1), max(area1))
"""
These two values were picked by trial and error. You'll 
likely want different values for each estimator.
"""

n1, bin_edges1, _ = plt.hist(area1, bins=num_bins1, range=bin_range1, color='k', histtype='step', label='Data')
# This plots the histogram AND saves the counts and bin_edges for later use

plt.xlabel('Pulse Amplitude(mV)')
plt.ylabel('Number of events')
plt.xlim(bin_range1)
# If the legend covers some data, increase the plt.xlim value, maybe (0,0.5)

bin_centers1 = 0.5 * (bin_edges1[1:] + bin_edges1[:-1])
"""
This gives us the x-data which are the centres of each bin.
This is visually better for plotting errorbars.
More important, it's the correct thing to do for fitting the
Gaussian to our histogram.
It also fixes the shape -- len(n1) < len(bin_edges1) so we
cannot use 
plt.plot(n1, bin_edges1)
as it will give us a shape error.
"""

sig1 = np.sqrt(n1)
sig1 = np.where(sig1 == 0, 1, sig1)
# The uncertainty on 0 count is 1, not 0. Replace all 0s with 1s.

plt.errorbar(bin_centers1, n1, yerr=sig1, fmt='none', c='k')
# This adds errorbars to the histograms, where each uncertainty is sqrt(y)

popt1, pcov1 = curve_fit(myGauss, bin_centers1, n1,
                         sigma=sig1, p0=(60, 50, 10, 5), absolute_sigma=True)
n1_fit = myGauss(bin_centers1, *popt1)
"""
n1_fit is our best fit line using our data points.
Note that if you have few enough bins, this best fit
line will have visible bends which look bad, so you
should not plot n1_fit directly. See below.
"""

chisquared1 = np.sum(((n1 - n1_fit) / sig1) ** 2)
dof1 = num_bins1 - len(popt1)
# Number of degrees of freedom is the number of data points less the number of fitted parameters

x_bestfit1 = np.linspace(bin_edges1[0], bin_edges1[-1], 1000)
y_bestfit1 = myGauss(x_bestfit1, *popt1)
# Best fit line smoothed with 1000 datapoints. Don't use best fit lines with 5 or 10 data points!

fontsize = 18
plt.plot(x_bestfit1, y_bestfit1, label='Fit')
print("Before calibration")
print(bin_range1)
print("Amp = ", popt1[0], "err = ", np.sqrt(pcov1[0][0]))
print("mean = ", popt1[1], "err = ", np.sqrt(pcov1[1][1]))
print("sigma= ", popt1[2], "err = ", np.sqrt(pcov1[2][2]))
print("reduces_chi=", chisquared1/dof1)
plt.legend(loc='upper right')
fontsize = 14
plt.text(-95, 55, r'$\mu$ = %3.2f keV' % (popt1[1]), fontsize=fontsize)
plt.text(-95, 50, r'$\chi^2$/DOF=', fontsize=fontsize)
plt.text(-95, 45, r'%3.2f/%i' % (chisquared1, dof1), fontsize=fontsize)
plt.text(-95, 40, r'$\sigma$ = %3.2f keV' % (popt1[2]), fontsize=fontsize)
plt.text(-95, 35, r'$\chi^2$ prob.= %1.1f' % (1 - chi2.cdf(chisquared1, dof1)), fontsize=fontsize)
plt.savefig("Plots/area1_pre.png", dpi=200)
plt.show()

"""
Area1 calibration
"""
c_factor = 10 / popt1[1]  # in keV/mV
c_factor_error = 10 * (np.sqrt(pcov1[1][1])) / (popt1[1]**2)
print("c_factor = ",c_factor, "keV/mV error = ", c_factor_error)
area1 *= c_factor
num_bins1 = 60
bin_range1 = ((min(area1) * c_factor) - 20 , (max(area1) * c_factor) + 30)
n1, bin_edges1, _ = plt.hist(area1, bins=num_bins1, range=bin_range1, color='k', histtype='step', label='Data')
# This plots the histogram AND saves the counts and bin_edges for later use

plt.xlabel('Particle Energy(keV)')
plt.ylabel('Number of events')
plt.xlim(bin_range1)
# If the legend covers some data, increase the plt.xlim value, maybe (0,0.5)

bin_centers1 = 0.5 * (bin_edges1[1:] + bin_edges1[:-1])
sig1 = np.sqrt(n1)
sig1 = np.where(sig1 == 0, 1, sig1)
# The uncertainty on 0 count is 1, not 0. Replace all 0s with 1s.

plt.errorbar(bin_centers1, n1, yerr=sig1, fmt='none', c='k')
# This adds errorbars to the histograms, where each uncertainty is sqrt(y)

popt1, pcov1 = curve_fit(myGauss, bin_centers1, n1,
                         sigma=sig1, p0=(50, 10, 2, 5), absolute_sigma=True)
n1_fit = myGauss(bin_centers1, *popt1)
"""
n1_fit is our best fit line using our data points.
Note that if you have few enough bins, this best fit
line will have visible bends which look bad, so you
should not plot n1_fit directly. See below.
"""

chisquared1 = np.sum(((n1 - n1_fit) / sig1) ** 2)
dof1 = num_bins1 - len(popt1)
# Number of degrees of freedom is the number of data points less the number of fitted parameters

x_bestfit1 = np.linspace(bin_edges1[0], bin_edges1[-1], 1000)
y_bestfit1 = myGauss(x_bestfit1, *popt1)
# Best fit line smoothed with 1000 datapoints. Don't use best fit lines with 5 or 10 data points!
plt.plot(x_bestfit1, y_bestfit1, label='Fit')
print("After calibration")
print(bin_range1)
print("Amp = ", popt1[0], "err = ", np.sqrt(pcov1[0][0]))
print("mean = ", popt1[1], "err = ", np.sqrt(pcov1[1][1]))
print("sigma= ", popt1[2], "err = ", np.sqrt(pcov1[2][2]))
print("reduces_chi=", chisquared1/dof1)
plt.legend(loc='upper right')
fontsize = 14
plt.text(-22, 55, r'$\mu$ = %3.2f keV' % (popt1[1]), fontsize=fontsize)
plt.text(-22, 50, r'$\chi^2$/DOF=', fontsize=fontsize)
plt.text(-22, 45, r'%3.2f/%i' % (chisquared1, dof1), fontsize=fontsize)
plt.text(-22, 40, r'$\sigma$ = %3.2f keV' % (popt1[2]), fontsize=fontsize)
plt.text(-22, 35, r'$\chi^2$ prob.= %1.1f' % (1 - chi2.cdf(chisquared1, dof1)), fontsize=fontsize)
plt.savefig("Plots/area1_post.png", dpi=200)
plt.show()

"""Residuals"""
# Compute residuals for the Gaussian fit after calibration
residuals = n1 - myGauss(bin_centers1, *popt1)
plt.xlim(bin_range1)
# Plot the residuals with black markers and error bars using '.' marker format
plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.errorbar(bin_centers1, residuals, yerr=sig1, fmt='.k', label='Residuals')
plt.xlabel('Particle Energy (keV)')
plt.ylabel('Residuals')
plt.legend()
plt.savefig("new_plots/area1_residuals", dpi=200)
plt.show()