import numpy as np



if __name__ == "__main__":
    equil_empty = 0.1790957545040282
    equil_empty_err= 6.707977268824151e-06
    T_empty = 304.582059688089
    T_empty_err = 1.200669150414119e-07

    equil_counterclock = 0.16041512319248966
    equil_counterclock_err= 2.1515513960925743e-05
    T_counterclock = 304.73663015929856

    T = (T_empty + T_counterclock) / 2
    print("T = ", T)


    s = np.abs(equil_counterclock - equil_empty)
    s_err = np.sqrt((equil_empty_err ** 2) + (equil_counterclock_err ** 2))

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
    x1_err = np.sqrt((r1_err ** 2) + ((d*np.cos(theta)*theta_err) ** 2))
    x2 = r2 + (w / 2) - (d * np.sin(theta))
    x2_err = np.sqrt((r2_err ** 2) + ((d * np.cos(theta) * theta_err) ** 2))
    print("x1 = ", x1, "x1_err =", x1_err)
    print("x2 = ", x2, "x1_err =", x2_err)

    G_1 = 2 * (np.pi**2) * (d**2) * (x1 ** 2) * theta / (m1 * np.cos(theta) * (T**2))
    G_2 = 2 * (np.pi**2) * (d**2) * (x2 ** 2) * theta / (m2 * np.cos(theta) * (T**2))
    print("G1 =", G_1)
    print("G2 =", G_2)


