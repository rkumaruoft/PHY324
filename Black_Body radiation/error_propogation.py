import numpy as np
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

def error_in_lambda(lambda_error, lambda_val):
    return lambda_error / (2 * lambda_val)
