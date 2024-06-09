from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize as sopt


# Define constants
PI = 2 * np.pi
C0 = 299792458 # m/s
MU0 = 4 * PI * 1e-7

def solve_symm_slab(
    d: float,
    ncl: float,
    nco: float,
    lam0: float,
    pol: Literal["TE", "TM"],
    show_plot: bool = False,
):
    if pol not in ["TE", "TM"]:
        raise ValueError("Unknown polarization. Please use 'TE' or 'TM'.")

    k0 = 2 * PI / lam0
    omega0 = 2 * PI * C0 / lam0

    # Radius of circle u**2 + v**2
    V = (PI * d / lam0) * np.sqrt(nco**2 - ncl**2)

    # Find # of (TE for now) modes based on radius
    n_modes = int(np.floor((V + PI / 2) / (PI / 2)))

    # Plotting transcendentals
    if show_plot is True:
        # Coordinates for a circle
        NUM_PTS = 201
        rad = np.linspace(0, PI/4, num=NUM_PTS)
        x_V = V * np.cos(rad)
        y_V = V * np.sin(rad)

        # v = utan(u) and v = -utan(u)
        u1 = np.linspace(0, PI/4, num=NUM_PTS)
        u2 = u1.copy() + PI / 2
        u3 = u1.copy() + PI
        u4 = u1.copy() + 3 / 2 *  PI

        v1_te = u1 * np.tan(u1)
        v2_te = -u2 / np.tan(u2)
        v3_te = u3 * np.tan(u3)
        v4_te = -u4 / np.tan(u4)

        # Plotting
        plt.figure()
        plt.plot(x_V, y_V)
        plt.plot(u1, v1_te, label='TE even')
        plt.plot(u2, v2_te, label='TE odd')
        plt.plot(u3, v3_te, label='TE even')
        plt.plot(u4, v4_te, label='TE odd')

        plt.ylim((0 , 2 * PI))
        # plt.xlim([0, np.max(u4)])
        plt.xlim(left=0)
        plt.gca().set_aspect('equal')
        plt.xlabel("u")
        plt.ylabel("u")
        plt.legend()

        plt.show()

    # pre-solve loop init
    x_lims = [0, PI / 2]
    neff = []
    beta = []
    alpha = []
    kx = []

    # Initialize field coordinates and variables
    x = np.linspace(-2 * d, 2 * d, 100)
    Ex = np.zeros((n_modes, len(x)), dtype=complex) # Field, dimensions of mode x field
    Ey = np.zeros((n_modes, len(x)), dtype=complex)
    Ez = np.zeros((n_modes, len(x)), dtype=complex)
    Hx = np.zeros((n_modes, len(x)), dtype=complex) # Field, dimensions of mode x field
    Hy = np.zeros((n_modes, len(x)), dtype=complex)
    Hz = np.zeros((n_modes, len(x)), dtype=complex)

    # # # # # # # # # # # #
    # Solve for the modes
    # # # # # # # # # # # #

    for ii in range(1, n_modes+1): # ATTENTION Python counter begin to zero
        if pol == "TE":
            # Find the intersects of the mode/wave equations
            # the intersect is "u"
            if ii % 2 == 1: # "Even" mode takes on odd mode
                def symslab_TE_even(u, V):
                    return V**2 - u**2 - (u * np.tan(u))**2
                sol = sopt.root_scalar(
                    lambda x: symslab_TE_even(x, V), bracket=x_lims
                )
                u = sol.root
            else:
                def symslab_TE_odd(u, V):
                    return V**2 - u**2 - (u / np.tan(u))**2
                sol = sopt.root_scalar(
                    lambda x: symslab_TE_odd(x, V), bracket=x_lims
                )
                u = sol.root
        elif pol == "TM":
            # Find the intersects of the mode/wave equations
            # the intersect is "u"
            if ii % 2 == 1: # "Even" mode takes on odd mode
                def symslab_TM_even(u, V, n1, n2):
                    return V**2 - u**2 - (((n1**2) / (n2**2)) * u * np.tan(u))**2
                sol = sopt.root_scalar(
                    lambda x: symslab_TM_even(x, V, ncl, nco), bracket=x_lims
                )
                u = sol.root
            else:
                def symslab_TM_odd(u, V, n1, n2):
                    return V**2 - u**2 - (((n1**2) / (n2**2)) * u / np.tan(u))**2
                sol = sopt.root_scalar(
                    lambda x: symslab_TM_odd(x, V, ncl, nco), bracket=x_lims
                )
                u = sol.root

        # Convert "u" to neff and beta
        kx.append(2 * u / d)
        beta_ii = np.sqrt(-kx[-1]**2 + (nco * k0) ** 2) # beta
        beta.append(beta_ii)
        neff_ii = beta_ii / k0 # neff
        neff.append(neff_ii)
        # if pol == "TE":
        #    # Calculate A, B, C, D
        #     if ii % 2 == 1: # "Even" mode takes on odd mode
        #         alpha.append(kx[-1] * np.tan(kx[-1] * d / 2))
        #         C = np.cos(kx[-1] * d / 2) / np.exp(-alpha[-1] * d / 2)
        #         D = C
        #     else:
        #         alpha.append(-kx[-1] / np.tan(kx[-1] * d / 2))
        #         C = np.sin(kx[-1] * d / 2) / np.exp(-alpha[-1] * d / 2)
        #         D = -C
        # elif pol == "TM":
        #    # Calculate A, B, C, D
        #     if ii % 2 == 1: # "Even" mode takes on odd mode
        #         alpha.append((ncl / nco)**2 * kx[-1] * np.tan(kx[-1] * d / 2))
        #         C = np.cos(kx[-1] * d / 2) / np.exp(-alpha[-1] * d / 2)
        #         D = C
        #     else:
        #         alpha.append(-(ncl / nco)**2 * kx[-1] / np.tan(kx[-1] * d / 2))
        #         C = np.sin(kx[-1] * d / 2) / np.exp(-alpha[-1] * d / 2)
        #         D = -C

        # Calculate A, B, C, D
        if ii % 2 == 1: # "Even" mode takes on odd mode
            alpha.append(kx[-1] * np.tan(kx[-1] * d / 2))
            C = np.cos(kx[-1] * d / 2) / np.exp(-alpha[-1] * d / 2)
            D = C
        else:
            alpha.append(-kx[-1] / np.tan(kx[-1] * d / 2))
            C = np.sin(kx[-1] * d / 2) / np.exp(-alpha[-1] * d / 2)
            D = -C
        if pol == "TM":
            # Negative sign came from the non-symmetry of the Maxwell equation
            alpha[-1] *= -(ncl / nco)**2

        # Solving for field distributions
        # Supplementary step to index on the boundaries
        idleft = np.argmin(np.abs(x + d / 2))
        idright= np.argmin(np.abs(x - d / 2))
        if pol == "TE": # Solve for the TE modes
            # Make E field
            Ey[ii-1, :idleft] = D * np.exp(alpha[-1] * x[:idleft])
            Ey[ii-1, idleft:idright+1] = np.cos(kx[-1] * x[idleft:idright+1])
            Ey[ii-1, idright+1:] = C * np.exp(-alpha[-1] * x[idright+1:])

            # Make H field
            Hx[ii-1, :] = -(beta_ii / (omega0 * MU0)) * Ey[ii-1, :]
            Hz[ii-1, :idleft] = 1j / (omega0 * MU0) * alpha[-1] * D * np.exp(alpha[-1] * x[:idleft])
            Hz[ii-1, idleft:idright+1] = -1j / (omega0 * MU0) * kx[-1] * np.sin(kx[-1] * x[idleft:idright+1])
            Hz[ii-1, idright+1:] = -1j / (omega0 * MU0) * alpha[-1] * C * np.exp(-alpha[-1] * x[idright+1:])
        elif pol == "TM":
            # Make H field
            Hy[ii-1, :idleft] = D * np.exp(alpha[-1] * x[:idleft])
            Hy[ii-1, idleft:idright+1] = np.cos(kx[-1] * x[idleft:idright+1])
            Hy[ii-1, idright+1:] = C * np.exp(-alpha[-1] * x[idright+1:])

            # Make H field
            Ex[ii-1, :] = -(beta_ii / (omega0 * MU0)) * Hy[ii-1, :]
            Ez[ii-1, :idleft] = 1j / (omega0 * MU0) * alpha[-1] * D * np.exp(alpha[-1] * x[:idleft])
            Ez[ii-1, idleft:idright+1] = -1j / (omega0 * MU0) * kx[-1] * np.sin(kx[-1] * x[idleft:idright+1])
            Ez[ii-1, idright+1:] = -1j / (omega0 * MU0) * alpha[-1] * C * np.exp(-alpha[-1] * x[idright+1:])

        x_lims = [ii * PI / 2, (ii + 1) * PI / 2] # Move to next interval
    
    Field = {
        "x": x,
        "Ex": Ex,
        "Ey": Ey,
        "Ez": Ez,
        "Hx": Hx,
        "Hy": Hy,
        "Hz": Hz,
    }

    return neff, beta, kx, alpha, Field
