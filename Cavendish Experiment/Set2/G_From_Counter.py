import numpy as np


def propagate_uncertainty_G(M1, M2, x1, x2, T, theta, d,
                            sigma_M1, sigma_M2, sigma_x1, sigma_x2, sigma_T, sigma_theta):
    """
    Computes the propagated uncertainty in G based on uncertainties in the measured quantities.

    Parameters:
    M1, M2: Masses of the large spheres
    x1, x2: Distances between the large and small masses
    T: Oscillation period
    theta: Angular displacement (in radians)
    d: Distance from the center of the rod to the small masses
    sigma_M1, sigma_M2: Uncertainties in M1 and M2
    sigma_x1, sigma_x2: Uncertainties in x1 and x2
    sigma_T: Uncertainty in T
    sigma_theta: Uncertainty in theta

    Returns:
    sigma_G: Propagated uncertainty in G
    """
    import numpy as np

    # Compute the denominator term (which appears in multiple places)
    denominator = T ** 3 * np.cos(theta) ** 2 * (M1 * x2 ** 2 + M2 * x1 ** 2) ** 2

    # First term: Contribution from uncertainty in theta
    term_theta = (T ** 2 * sigma_theta ** 2 * x1 ** 2 * x2 ** 2 * (M1 * x2 ** 2 + M2 * x1 ** 2) ** 2 * (
            theta * np.sin(theta) + np.cos(theta)) ** 2)

    # Second term: Contribution from uncertainties in M1, M2, x1, x2
    term_mass_distance = (T ** 2 * theta ** 2 * np.cos(theta) ** 2 *
                          (4 * M1 ** 2 * sigma_x1 ** 2 * x2 ** 6 +
                           4 * M2 ** 2 * sigma_x2 ** 2 * x1 ** 6 +
                           sigma_M1 ** 2 * x1 ** 2 * x2 ** 6 +
                           sigma_M2 ** 2 * x1 ** 6 * x2 ** 2))

    # Third term: Contribution from uncertainty in T
    term_T = (4 * sigma_T ** 2 * theta ** 2 * x1 ** 2 * x2 ** 2 * (M1 * x2 ** 2 + M2 * x1 ** 2) ** 2 * np.cos(
        theta) ** 2)

    # Compute propagated uncertainty in G
    sigma_G = (8 * np.pi ** 2 / denominator) * np.sqrt(
        d ** 2 * x1 ** 2 * x2 ** 2 * (term_theta + term_mass_distance + term_T))

    return sigma_G


if __name__ == "__main__":
    equil_empty = 0.06056076966601571
    equil_empty_err = np.sqrt((1.2037739120807117e-05 ** 2) + (0.00539 ** 2))
    T_empty = 305.4480989240732
    T_empty_err = 0.0990476308033528

    equil_counterclock = 0.16041512319248966
    equil_counterclock_err = 2.1515513960925743e-05
    T_counterclock = 304.73663015929856

    T = (T_empty + T_counterclock) / 2
    print("T = ", T)

    s = np.abs(equil_counterclock - equil_empty)
    s_err = np.sqrt((equil_empty_err ** 2) + (equil_counterclock_err ** 2))

    L = 4.47291  # meters
    L_err = 0.001

    w = 0.05

    d = 0.05

    theta = s / (2 * L)
    theta_err = theta * np.sqrt(((s_err / s) ** 2) + ((L_err / L) ** 2))

    print("Theta = ", theta, "err = ", theta_err)

    """Ball Specs"""
    m1, m1_err = 1444.5e-3, 0.1e-3
    m2, m2_err = 1467.9e-3, 0.1e-3
    r1, r1_err = 63.72e-3, 0.01e-3
    r2, r2_err = 63.89e-3, 0.01e-3

    x1 = r1 + (w / 2) - (d * np.sin(theta))
    x1_err = np.sqrt((r1_err ** 2) + ((d * np.cos(theta) * theta_err) ** 2))
    x2 = r2 + (w / 2) - (d * np.sin(theta))
    x2_err = np.sqrt((r2_err ** 2) + ((d * np.cos(theta) * theta_err) ** 2))
    print("x1 = ", x1, "x1_err =", x1_err)
    print("x2 = ", x2, "x1_err =", x2_err)

    G_1 = 2 * (np.pi ** 2) * (d ** 2) * (x1 ** 2) * theta / (m1 * np.cos(theta) * (T ** 2))
    G_2 = 2 * (np.pi ** 2) * (d ** 2) * (x2 ** 2) * theta / (m2 * np.cos(theta) * (T ** 2))
    print("G1 =", G_1)
    print("G2 =", G_2)
