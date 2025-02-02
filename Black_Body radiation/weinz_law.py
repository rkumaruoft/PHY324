import math
import pickle
import matplotlib.pyplot as plt
import numpy
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2
import csv


def calculate_temperature(voltage, current, Resist_0, Temp_0, alpha0):
    R = voltage / current
    Temp = Temp_0 + ((R / Resist_0) - 1) / alpha0
    return Temp


if __name__ == "__main__":
    voltage = []
    lambda_max = []
    temp_0 = []
    current = []
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

    temps = []
    for i in range(len(temp_0)):
        this_temp = calculate_temperature(voltage[i], current[i], 1.1, 293, 4.5e-3)
        print(temp_0[i], " ", this_temp)
        temps.append(this_temp)

    plt.scatter(lambda_max, temps)
    plt.title("wavelength vs temps Raw")
    plt.show()
