import numpy as np

import numpy as np

def uncertainty_G(x, theta, T, d, m, sigma_x, sigma_theta, sigma_T):
    """
    Computes the uncertainty in G using error propagation.

    Parameters:
    x : float  -> Position value
    theta : float  -> Angle in radians
    T : float  -> Period value
    d : float  -> Constant parameter d
    m : float  -> Mass value
    sigma_x : float  -> Uncertainty in x
    sigma_theta : float  -> Uncertainty in theta
    sigma_T : float  -> Uncertainty in T

    Returns:
    sigma_G : float  -> Uncertainty in G
    """
    # Compute individual terms
    term1 = (T**2 * sigma_theta**2 * x**2 * (theta * np.sin(theta) + np.cos(theta))**2)
    term2 = (4 * T**2 * sigma_x**2 * theta**2 * np.cos(theta)**2)
    term3 = (4 * sigma_T**2 * theta**2 * x**2 * np.cos(theta)**2)

    # Compute final uncertainty
    sigma_G = (2 * np.pi**2 * d**2 * x / (T**3 * m * np.cos(theta)**2)) * np.sqrt(term1 + term2 + term3)

    return sigma_G


if __name__ == "__main__":
    equil_empty = 0.1790957545040282
    equil_empty_err= 6.707977268824151e-06
    T_empty = 304.582059688089
    T_empty_err = 1.200669150414119e-07

    # equil_clock = 0.20048868845515166
    # equil_clock_err= 8.23696432920053e-06
    # T_clock = 304.73663015929856
    # T_clock_err = 3.045829831611369e-07


    """From Averaged dataset"""
    equil_clock = 0.20072349644038168
    equil_clock_err = 2.6975467279274043e-05
    T_clock = 304.73461082572953
    T_clock_err =  2.324956493044523e-06

    T = (T_empty + T_clock) / 2
    T_err = np.sqrt(T_empty_err**2 + T_clock_err**2) / 2
    print("T = ", T, "err=", T_err)

    s = equil_clock - equil_empty
    s_err = np.sqrt((equil_empty_err ** 2) + (equil_clock_err ** 2))

    L = 4.47291 #meters
    L_err = 0.001

    w = 0.05

    d = 0.05

    theta = s / (2*L)
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

    G_1 = 2 * (np.pi**2) * (d**2) * (x1 ** 2) * theta / (m1 * np.cos(theta) * (T**2))
    G_1_err = uncertainty_G(x1, theta, T, d, m1, x1_err, theta_err,T_err)
    G_2 = 2 * (np.pi**2) * (d**2) * (x2 ** 2) * theta / (m2 * np.cos(theta) * (T**2))
    G_2_err = uncertainty_G(x2, theta, T, d, m2, x2_err, theta_err, T_err)
    print("G1 =", G_1, "err =", G_1_err)
    print("G2 =", G_2, "err =", G_2_err)


