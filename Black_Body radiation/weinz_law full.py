import math
import pickle
import matplotlib.pyplot as plt
import numpy
from matplotlib import rc
import numpy as np
from numpy.ma.extras import average
from scipy.optimize import curve_fit
from scipy.stats import chi2
import csv

def calculate_temperature(voltage, current, Resist_0, Temp_0, alpha0):
    R = voltage / current
    Temp = Temp_0 + ((R / Resist_0) - 1) / alpha0
    return Temp

def weinz_law(x, A, B, C):
    return (A / (x + B)) + C

def d_lambda2_d_theta(theta_angle, A, B):
    theta_rad = np.radians(theta_angle)  # Convert degrees to radians

    # Compute k(theta)
    k_theta = (((2 / np.sqrt(3)) * np.sin(theta_rad) + (1 / 2)) ** 2) + 3 / 4

    # Compute dk/dtheta
    dk_dtheta = (4 / np.sqrt(3)) * np.cos(theta_rad) * ((2 / np.sqrt(3)) * np.sin(theta_rad) + (1 / 2))

    # Compute the derivative d(lambda^2)/d(theta)
    d_lambda2 = - (A / (np.sqrt(k_theta) - B) ** 2) * (1 / (2 * np.sqrt(k_theta))) * dk_dtheta

    return d_lambda2

# Define function to compute error propagation
def compute_lambda2_error(theta_values, theta_errors, A, B):
    d_lambda2_vals = d_lambda2_d_theta(theta_values, A, B)
    lambda2_errors = np.abs(d_lambda2_vals) * theta_errors  # Propagate error using absolute derivative
    return lambda2_errors



if __name__ == "__main__":
    voltage = []
    lambda_max = []
    temp_0 = []
    current = []
    lambda_error = []
    with open('data/weins law data.csv', mode='r') as file:
        # Create a CSV reader object
        csv_reader = csv.reader(file)
        # Skip the header row (if there is one)
        next(csv_reader, None)
        # Iterate over each row in the CSV file
        for row in csv_reader:
            voltage.append(float(row[1]))
            lambda_max.append(float(row[2]))
            temp_0.append(float(row[3]))
            current.append(float(row[5]))
            lambda_error.append(float(row[6]))

    temps = []
    for i in range(len(temp_0)):
        this_temp = calculate_temperature(voltage[i], current[i], 1.1, 293, 4.5 * (10**-3))
        temps.append(this_temp)


    avg_lambdas = []
    avg_temps = []
    avg_lambda_errors = []
    index = 0
    """Averages of lambda"""
    for i in range(7):
        this_avg_lambda = []
        this_avg_temp = []
        this_errors = []
        for j in range(3):
            this_avg_lambda.append(lambda_max[index])
            this_avg_temp.append(temps[index])
            this_errors.append(lambda_error[index])
            index += 1
        avg_lambdas.append(sum(this_avg_lambda) / len(this_avg_lambda))
        avg_temps.append(sum(this_avg_temp) / len(this_avg_temp))
        avg_lambda_errors.append(np.sqrt(sum([x**2 for x in this_errors])) / 3) #using uncertainty propagation
        # avg_lambda_errors.append(max(this_avg_lambda) - min(this_avg_lambda))
        # avg_lambda_errors.append(13.6)

    for i in range(len(avg_lambda_errors)):
        print(avg_lambdas[i], " ", avg_lambda_errors[i])

    avg_temps = np.array(temps)
    avg_lambdas = np.array(lambda_max)
    avg_lambda_errors = np.array(lambda_error)

    for i in range(len(avg_lambda_errors)):
        print(avg_lambdas[i], " ", avg_lambda_errors[i])

    plt.errorbar(avg_temps, avg_lambdas,yerr=avg_lambda_errors, fmt='.', c='k', label="Data")
    plt.ylabel("Wavelength")
    plt.xlabel("Temperature")

    popt, pcov = curve_fit(weinz_law, xdata=avg_temps, ydata=avg_lambdas, sigma=avg_lambda_errors)
    print("popt ", popt)
    y_data = weinz_law(avg_temps, *popt)
    plt.plot(avg_temps, y_data, label="Fit curve")

    """For al the temps"""
    new_temps = np.array(range(2400,3500))
    y1 = weinz_law(new_temps, *popt)
    plt.plot(new_temps, y1, label="Fit curve1")

    plt.legend()
    plt.show()

    residuals = avg_lambdas - y_data
    plt.errorbar(avg_temps, residuals, yerr=avg_lambda_errors, fmt='.', c="r")
    plt.axhline(y=0)



    # Compute chi-squared statistic
    chi_squared = np.sum(((avg_lambdas - y_data) / avg_lambda_errors) ** 2)

    # Degrees of freedom (number of data points - number of fit parameters)
    dof = len(avg_lambdas) - len(popt)

    # Compute chi-squared probability (p-value)
    chi_prob = 1 - chi2.cdf(chi_squared, dof)
    print(f"Chi-Squared Probability: {chi_prob:.5f}")
    # Compute reduced chi-squared (chi-squared per degree of freedom)
    reduced_chi_squared = chi_squared / dof
    # Print result
    print(f"Reduced Chi-Squared: {reduced_chi_squared:.2f}")


    plt.show()

