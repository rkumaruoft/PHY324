import numpy as np
from scipy.optimize import curve_fit


def propagate_uncertainty(x, y, sigma_x, sigma_y):
    """
    Calculate the uncertainty in the angle (theta) based on uncertainties in x and y.

    Parameters:
    x (float or np.array): The x position(s).
    y (float or np.array): The y position(s).
    sigma_x (float or np.array): The uncertainty in x.
    sigma_y (float or np.array): The uncertainty in y.

    Returns:
    sigma_theta (float or np.array): The uncertainty in the angle theta.
    """
    # Partial derivatives
    dtheta_dx = -y / (x ** 2 + y ** 2)
    dtheta_dy = x / (x ** 2 + y ** 2)

    # Propagate the uncertainty
    sigma_theta = (np.sqrt((dtheta_dx * sigma_x) ** 2 + (dtheta_dy * sigma_y) ** 2))

    return np.degrees(sigma_theta)


def calculate_theta_from_neg_y(x, y):
    # Calculate the angle in radians using np.arctan2, which is relative to the positive x-axis
    return np.degrees(np.arctan2(x, -y))


def damped_cosine(t, A, omega, phi, tau, B):
    return A * np.cos(omega * t + phi) * np.exp(-t / tau) + B

def fit_damped_cosine(section_times, section_thetas, section_uncert, p0):
    popt, pcov = curve_fit(damped_cosine, section_times, section_thetas, sigma=section_uncert,maxfev=100000000)
    return popt, pcov  # Return the optimized parameters
