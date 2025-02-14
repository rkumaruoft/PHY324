import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import Period
from scipy.optimize import curve_fit


def convert_xy_to_theta(data):
    t = data['t'].values
    x = data['x'].values
    y = data['y'].values
    theta = np.arctan2(x, -y)
    return pd.DataFrame({'t': t, 'theta': theta})


def damped_cosine(t, theta_0, tau, T, phi_0):
    phi_0_rad = np.radians(phi_0)
    theta = theta_0 * np.exp(-t / tau) * np.cos((2 * np.pi / T) * t + phi_0_rad)
    return theta


# Load the CSV file
file_path = "t_vs_xy.csv"
data = pd.read_csv(file_path)
# Convert to theta vs. t
theta_data = convert_xy_to_theta(data)
# # Plot Theta vs Time
# plt.figure(figsize=(8, 4))
# plt.plot(theta_data['t'], theta_data['theta'], label="Theta vs Time", linestyle="none", marker=".", markersize=2
#          ,color='g')
# plt.xlabel("Time (t)")
# plt.ylabel("Theta (radians)")
# plt.title("Theta vs Time")
# plt.legend()
# plt.grid()
# plt.show()

# Convert theta from radians to degrees
theta_data['theta'] = np.degrees(theta_data['theta'])

# Plot Theta (degrees) vs Time
# plt.figure(figsize=(8, 4))
# plt.plot(theta_data['t'], theta_data['theta'], label="Theta (degrees) vs Time",
#          linestyle="none", marker=".", markersize=2,
#          color='g')
# plt.xlabel("Time (t)")
# plt.ylabel("Theta (degrees)")
# plt.title("Theta (degrees) vs Time")
# plt.grid()
# plt.show()

start_angles = []
periods = []
period_err = []
taus = []
tau_err = []
print("FUll section")
print("=============================================")
# Plot the sectioned theta data
plt.figure(figsize=(8, 4))
plt.plot(theta_data['t'], theta_data['theta'], linestyle="none", marker=".", markersize=2,
         color='g')
plt.xlabel("Time (t)")
plt.ylabel("Theta (degrees)")
plt.title("Theta vs Time Full")
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = theta_data['theta'].iloc[0]  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=theta_data["t"],ydata=theta_data['theta'],
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(theta_data['t'], *popt)
plt.plot(theta_data['t'], fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", popt[0])
print("Period = ", popt[2], "Tau = ", popt[1])
start_angles.append(popt[0])
periods.append(popt[2])
period_err.append(np.sqrt(pcov[2][2]))
taus.append(popt[1])
tau_err.append(np.sqrt(pcov[2][2]))
plt.grid()
plt.show()

print("=============================================")
theta_sec = theta_data[theta_data['t'] >= 20].reset_index(drop=True)

# Shift time so that t starts from 0 again
theta_sec['t'] -= theta_sec['t'].iloc[0]

# Plot the updated sectioned theta data
plt.figure(figsize=(8, 4))
plt.plot(theta_sec['t'], theta_sec['theta'], linestyle="none", marker=".", markersize=2,
         color='b', label="Theta Data (After 20s Cut)")
plt.xlabel("Time (s)")
plt.ylabel("Theta (°)")
plt.title("Theta vs Time (After 20s Removal)")
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = theta_sec['theta'].iloc[0]  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=theta_sec["t"],ydata=theta_sec['theta'],
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(theta_sec['t'], *popt)
plt.plot(theta_sec['t'], fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", popt[0])
print("Period = ", popt[2], "Tau = ", popt[1])
start_angles.append(popt[0])
periods.append(popt[2])
period_err.append(np.sqrt(pcov[2][2]))
taus.append(popt[1])
tau_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.grid()
plt.show()

print("=============================================")
theta_sec = theta_data[theta_data['t'] >= 30].reset_index(drop=True)

# Shift time so that t starts from 0 again
theta_sec['t'] -= theta_sec['t'].iloc[0]

# Plot the updated sectioned theta data
plt.figure(figsize=(8, 4))
plt.plot(theta_sec['t'], theta_sec['theta'], linestyle="none", marker=".", markersize=2,
         color='b', label="Theta Data (After 30s Cut)")
plt.xlabel("Time (s)")
plt.ylabel("Theta (°)")
plt.title("Theta vs Time (After 30s Removal)")
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = theta_sec['theta'].iloc[0]  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=theta_sec["t"],ydata=theta_sec['theta'],
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(theta_sec['t'], *popt)
plt.plot(theta_sec['t'], fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", popt[0])
print("Period = ", popt[2], "Tau = ", popt[1])
start_angles.append(popt[0])
periods.append(popt[2])
period_err.append(np.sqrt(pcov[2][2]))
taus.append(popt[1])
tau_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.grid()
plt.show()

print("=============================================")
this_sec = 40
theta_sec = theta_data[theta_data['t'] >= this_sec].reset_index(drop=True)

# Shift time so that t starts from 0 again
theta_sec['t'] -= theta_sec['t'].iloc[0]

# Plot the updated sectioned theta data
plt.figure(figsize=(8, 4))
plt.plot(theta_sec['t'], theta_sec['theta'], linestyle="none", marker=".", markersize=2,
         color='b', label="Theta Data (After "+ str(this_sec) +"s Cut)")
plt.xlabel("Time (s)")
plt.ylabel("Theta (°)")
plt.title("Theta vs Time (After "+ str(this_sec) +" Removal)")
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = theta_sec['theta'].iloc[0]  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=theta_sec["t"],ydata=theta_sec['theta'],
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(theta_sec['t'], *popt)
plt.plot(theta_sec['t'], fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", popt[0])
print("Period = ", popt[2], "Tau = ", popt[1])
start_angles.append(popt[0])
periods.append(popt[2])
period_err.append(np.sqrt(pcov[2][2]))
taus.append(popt[1])
tau_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.grid()
plt.show()

print("=============================================")
this_sec = 50
theta_sec = theta_data[theta_data['t'] >= this_sec].reset_index(drop=True)

# Shift time so that t starts from 0 again
theta_sec['t'] -= theta_sec['t'].iloc[0]

# Plot the updated sectioned theta data
plt.figure(figsize=(8, 4))
plt.plot(theta_sec['t'], theta_sec['theta'], linestyle="none", marker=".", markersize=2,
         color='b', label="Theta Data (After "+ str(this_sec) +"s Cut)")
plt.xlabel("Time (s)")
plt.ylabel("Theta (°)")
plt.title("Theta vs Time (After "+ str(this_sec) +" Removal)")
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = theta_sec['theta'].iloc[0]  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=theta_sec["t"],ydata=theta_sec['theta'],
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(theta_sec['t'], *popt)
plt.plot(theta_sec['t'], fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", popt[0])
print("Period = ", popt[2], "Tau = ", popt[1])
start_angles.append(popt[0])
periods.append(popt[2])
period_err.append(np.sqrt(pcov[2][2]))
taus.append(popt[1])
tau_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.grid()
plt.show()

print("=============================================")
this_sec = 10
theta_sec = theta_data[theta_data['t'] >= this_sec].reset_index(drop=True)

# Shift time so that t starts from 0 again
theta_sec['t'] -= theta_sec['t'].iloc[0]

# Plot the updated sectioned theta data
plt.figure(figsize=(8, 4))
plt.plot(theta_sec['t'], theta_sec['theta'], linestyle="none", marker=".", markersize=2,
         color='b', label="Theta Data (After "+ str(this_sec) +"s Cut)")
plt.xlabel("Time (s)")
plt.ylabel("Theta (°)")
plt.title("Theta vs Time (After "+ str(this_sec) +" Removal)")
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = theta_sec['theta'].iloc[0]  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=theta_sec["t"],ydata=theta_sec['theta'],
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(theta_sec['t'], *popt)
plt.plot(theta_sec['t'], fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", popt[0])
print("Period = ", popt[2], "Tau = ", popt[1])
start_angles.append(popt[0])
periods.append(popt[2])
period_err.append(np.sqrt(pcov[2][2]))
taus.append(popt[1])
tau_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.grid()
plt.show()

print("=============================================")
this_sec = 60
theta_sec = theta_data[theta_data['t'] >= this_sec].reset_index(drop=True)

# Shift time so that t starts from 0 again
theta_sec['t'] -= theta_sec['t'].iloc[0]

# Plot the updated sectioned theta data
plt.figure(figsize=(8, 4))
plt.plot(theta_sec['t'], theta_sec['theta'], linestyle="none", marker=".", markersize=2,
         color='b', label="Theta Data (After "+ str(this_sec) +"s Cut)")
plt.xlabel("Time (s)")
plt.ylabel("Theta (°)")
plt.title("Theta vs Time (After "+ str(this_sec) +" Removal)")
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = theta_sec['theta'].iloc[0]  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=theta_sec["t"],ydata=theta_sec['theta'],
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(theta_sec['t'], *popt)
plt.plot(theta_sec['t'], fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", popt[0])
print("Period = ", popt[2], "Tau = ", popt[1])
start_angles.append(popt[0])
periods.append(popt[2])
period_err.append(np.sqrt(pcov[2][2]))
taus.append(popt[1])
tau_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.grid()
plt.show()

print("=============================================")
this_sec = 70
theta_sec = theta_data[theta_data['t'] >= this_sec].reset_index(drop=True)

# Shift time so that t starts from 0 again
theta_sec['t'] -= theta_sec['t'].iloc[0]

# Plot the updated sectioned theta data
plt.figure(figsize=(8, 4))
plt.plot(theta_sec['t'], theta_sec['theta'], linestyle="none", marker=".", markersize=2,
         color='b', label="Theta Data (After "+ str(this_sec) +"s Cut)")
plt.xlabel("Time (s)")
plt.ylabel("Theta (°)")
plt.title("Theta vs Time (After "+ str(this_sec) +" Removal)")
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = theta_sec['theta'].iloc[0]  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=theta_sec["t"],ydata=theta_sec['theta'],
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(theta_sec['t'], *popt)
plt.plot(theta_sec['t'], fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", popt[0])
print("Period = ", popt[2], "Tau = ", popt[1])
start_angles.append(popt[0])
periods.append(popt[2])
period_err.append(np.sqrt(pcov[2][2]))
taus.append(popt[1])
tau_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.grid()
plt.show()

print("=============================================")
this_sec = 80
theta_sec = theta_data[theta_data['t'] >= this_sec].reset_index(drop=True)

# Shift time so that t starts from 0 again
theta_sec['t'] -= theta_sec['t'].iloc[0]

# Plot the updated sectioned theta data
plt.figure(figsize=(8, 4))
plt.plot(theta_sec['t'], theta_sec['theta'], linestyle="none", marker=".", markersize=2,
         color='b', label="Theta Data (After "+ str(this_sec) +"s Cut)")
plt.xlabel("Time (s)")
plt.ylabel("Theta (°)")
plt.title("Theta vs Time (After "+ str(this_sec) +" Removal)")
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = theta_sec['theta'].iloc[0]  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=theta_sec["t"],ydata=theta_sec['theta'],
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(theta_sec['t'], *popt)
plt.plot(theta_sec['t'], fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", popt[0])
print("Period = ", popt[2], "Tau = ", popt[1])
start_angles.append(popt[0])
periods.append(popt[2])
period_err.append(np.sqrt(pcov[2][2]))
taus.append(popt[1])
tau_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.grid()
plt.show()

print("=============================================")
this_sec = 90
theta_sec = theta_data[theta_data['t'] >= this_sec].reset_index(drop=True)

# Shift time so that t starts from 0 again
theta_sec['t'] -= theta_sec['t'].iloc[0]

# Plot the updated sectioned theta data
plt.figure(figsize=(8, 4))
plt.plot(theta_sec['t'], theta_sec['theta'], linestyle="none", marker=".", markersize=2,
         color='b', label="Theta Data (After "+ str(this_sec) +"s Cut)")
plt.xlabel("Time (s)")
plt.ylabel("Theta (°)")
plt.title("Theta vs Time (After "+ str(this_sec) +" Removal)")
# Set initial guesses
T_guess = 1.420695  # Replace with estimated period if known
tau_guess = T_guess * 20  # Rough estimate for damping
theta_0_guess = theta_sec['theta'].iloc[0]  # First value in the sectioned data
phi_0_guess = 0  # Assuming oscillation starts at peak
popt, pcov = curve_fit(damped_cosine, xdata=theta_sec["t"],ydata=theta_sec['theta'],
                       p0=[theta_0_guess, tau_guess, T_guess, phi_0_guess])
fit_data = damped_cosine(theta_sec['t'], *popt)
plt.plot(theta_sec['t'], fit_data, linestyle="dashed", markersize=2,
         color='g')
print("Start_theta = ", theta_sec['theta'][0])
print("Period = ", popt[2], "Tau = ", popt[1])
start_angles.append(popt[0])
periods.append(popt[2])
period_err.append(np.sqrt(pcov[2][2]))
taus.append(popt[1])
tau_err.append(np.sqrt(pcov[2][2]))
plt.legend()
plt.grid()
plt.show()


plt.errorbar(start_angles,periods, yerr=period_err, linestyle="None")
plt.show()

plt.plot(start_angles,taus, marker=".", linestyle="None")
plt.show()