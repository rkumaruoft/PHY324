import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

# Simulated experimental data (example)
t_exp = np.linspace(0, 10, 50)  # Experimental time points
theta_exp = np.radians(45) * np.cos(np.sqrt(9.81 / 1.0) * t_exp) + np.random.normal(0, 0.05, size=len(t_exp))  # Simulated noisy measurements

# Define the system of ODEs
def pendulum(t, y, g, L):
    theta, omega = y
    dtheta_dt = omega
    domega_dt = - (g / L) * np.sin(theta)
    return [dtheta_dt, domega_dt]

# Function to solve the pendulum for given g, L
def solve_pendulum(t, g, L, theta_0, omega_0):
    sol = solve_ivp(pendulum, (t[0], t[-1]), [theta_0, omega_0], t_eval=t, args=(g, L), method='RK45')
    return sol

# Function to interpolate the numerical solution for curve fitting
def interpolated_pendulum(t, g, L):
    sol = solve_pendulum(t_exp, g, L, np.radians(45), 0.0)  # Solve using guessed parameters
    interp_func = interp1d(sol.t, sol.y[0], kind='cubic', fill_value="extrapolate")  # Interpolator
    return interp_func(t)  # Return interpolated values at t

# Fit the experimental data
popt, pcov = curve_fit(interpolated_pendulum, t_exp, theta_exp, p0=[9.81, 1.0])

# Extract fitted parameters
g_fit, L_fit = popt
print(f"Fitted Parameters: g = {g_fit:.3f}, L = {L_fit:.3f}")

# Generate the fitted curve
theta_fit = interpolated_pendulum(t_exp, g_fit, L_fit)

# Plot results
plt.figure(figsize=(8, 5))
plt.scatter(t_exp, theta_exp, label="Experimental Data", color="red")
plt.plot(t_exp, theta_fit, label="Fitted Model", linestyle="--", color="blue")
plt.xlabel('Time (s)')
plt.ylabel('Angle (radians)')
plt.title('Fitting Experimental Data to Numerical Pendulum Solution')
plt.legend()
plt.grid()
plt.show()
