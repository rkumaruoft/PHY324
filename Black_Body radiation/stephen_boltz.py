import math
import pickle
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2
import csv


intensity = []
temp = []
with open('data/Total Intensity vs. temperature.csv', mode='r') as file:
    # Create a CSV reader object
    csv_reader = csv.reader(file)
    # Skip the header row (if there is one)
    next(csv_reader, None)
    # Iterate over each row in the CSV file
    for row in csv_reader:
        intensity.append(float(row[5]))
        temp.append(float(row[6]) ** 4)


plt.scatter(temp, intensity)
plt.show()
