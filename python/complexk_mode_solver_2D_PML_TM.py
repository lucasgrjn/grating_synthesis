from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


# To get the same flattening as Matlab use the flag (same as Fortran)
MATLAB_FLAG = 'F'


def f_bloch_complexk_mode_solver_2D_PML(
        N: np.ndarray,
        disc: float,
        k0: float,
        num_modes: int,
        guessk: float,
        BC: Literal[0, 1],
        pol: Literal["TE", "TM"],
        pml_options: list[Literal[0, 1], float, float, int]
    ):

    # Total length of unwrapped vectors
    nx, nz = N.shape
    n_elem = nx * nz

    # Relative permittivity
    er = N ** 2
    # To generate yee we use the algorithm from the following link
    # https://stackoverflow.com/questions/7656665/how-to-repeat-elements-of-an-array-along-two-axes
    er_yee_temp = np.empty((nx, 2, nz, 2))
    er_yee_temp[...] = er[:, None, :, None]
    er_yee = er_yee_temp.reshape(nx * 2, nz * 2)

    # Draw in PMLs
    if pml_options[0] == 1:
        pml_len_nm = pml_options[1] # Length of pml in nm
        pml_str = pml_options[2]    # Strength of pml in complex plane
        pml_order = pml_options[3]  # PML polynomial order

    # Setup discretizations
    nx_pml = 2 * pml_len_nm // disc # Number of discretizations that pml spans, double sampled grid
    x_indx = np.arange(1, nx_pml)

    # Using polynomial strength pml
    pml_x = (1. + 1.j * pml_str * (x_indx / nx_pml)**pml_order)

    # Draw stretched coordinate pml
    pml_x_all = np.ones((2 * nx, nz), dtype=complex)
    pml_x_all[0:nx_pml-1, :] = np.repeat(np.flip(pml_x[0:nx_pml-1])[..., np.newaxis], repeats=nz, axis=1)
    pml_x_all[-(nx_pml-1):, :] = np.repeat(pml_x[..., np.newaxis], repeats=nz, axis=1)

    # Stretched coordinate operator
    pml_x_all_vec = pml_x_all
    # TODO: Verify the indexing
    Sx_f = sp.diags(1 / pml_x_all_vec[1::2,:].flatten(MATLAB_FLAG), 0, shape=(n_elem, n_elem))
    Sx_b = sp.diags(1 / pml_x_all_vec[0:-1:2,:].flatten(MATLAB_FLAG), 0, shape=(n_elem, n_elem))


    # Generate forward Dx
    # First generate vectors corresponding to the diagonals
    # Where the suffix # of the diag = # of diagonals from the middle
    # Default BC is TE/PEC or TM/TMC
    diag0 = -np.ones((n_elem, 1)) # Diag middle
    diag1 = np.ones((n_elem, 1)) # Diag plus 1
    diag1[nx-1::nx] = 0 # Don't carry over into the next column
    diagm1 = np.zeros((n_elem, 1))
    # Boundary condition for TE/PMC or TM/PEC
    if (BC == 1 and pol == "TE") or (BC == 0 and pol == "TM"):
        diagm1[nx-2::nx] = 1
    # Shift diag1
    diag1[1:] = diag1[:-1]
    # Stitch together the diags
    # Contrarely to Matlab we need to manually flatten
    diag_all = [
        diagm1.flatten(MATLAB_FLAG),
        diag0.flatten(MATLAB_FLAG),
        diag1.flatten(MATLAB_FLAG),
    ]
    diag_indexs = [-1, 0, 1]
    # Make sparse matrix
    # Dx_f = (1 / disc) * sp.diags(diag_all, diag_indexs, shape=(n_elem, n_elem))
    Dx_f = (1 / disc) * sp.diags(diag_all, diag_indexs, shape=(n_elem, n_elem))
    # Stretched coordinate pml
    if pml_options[0] == 1:
        Dx_f = Sx_f @ Dx_f


    # Generate forward Dx on +1/2 grid
    # Default BC is TE/PEC or TM/TMC
    # The other BC is not implemented yet
    diag0 = -np.ones((n_elem, 1)) # Diag middle
    diag0[nx-1::nx] = -2
    diag1 = np.ones((n_elem, 1)) # Diag plus 1
    diag0[nx-1::nx] = 0
    diagm1 = np.zeros((n_elem, 1))
    # BC NOT IMPLEMENTED
    diag1[1:] = diag1[0:-1]
    # Stitch together the diags
    # Contrarely to Matlab we need to manually flatten
    diag_all = [
        diagm1.flatten(MATLAB_FLAG),
        diag0.flatten(MATLAB_FLAG),
        diag1.flatten(MATLAB_FLAG),
    ]
    diag_indexs = [-1, 0, 1]
    # Make sparse matrix
    Dx_f_plushalf = (1 / disc) * sp.diags(diag_all, diag_indexs, shape=(n_elem, n_elem))
    # Stretched coordinate pml
    if pml_options[0] == 1:
        Dx_f_plushalf  = Sx_b @ Dx_f_plushalf 


    # Generate backwards Dx
    diag0 = np.ones((n_elem, 1))
    diagm1 = -np.ones((n_elem, 1))
    diagm1[nx-1::nx] = 0 # No need to shift due to being in lower triangle of Dx
    diag1 = np.zeros((n_elem, 1))
    # Boundary condition for TE/PMC or TM/PEC
    if (BC == 1 and pol == "TE") or (BC == 0 and pol == "TM"):
        diag1[0::nx] = -1
    # Shift diag1
    diag1[1:] = diag1[:-1]
    # Stitch together the diags
    # Contrarely to Matlab we need to manually flatten
    diag_all = [
        diagm1.flatten(MATLAB_FLAG),
        diag0.flatten(MATLAB_FLAG),
        diag1.flatten(MATLAB_FLAG)
    ]
    diag_indexs = [-1, 0, 1]
    # Make sparse matrix
    Dx_b = (1 / disc) * sp.diags(diag_all, diag_indexs, shape=(n_elem, n_elem))
    # Stretched coordinate pml
    if pml_options[0] == 1:
        Dx_b = Sx_b @ Dx_b

    # Generate Dx squared
    Dx2 = Dx_b @ Dx_f

    # Generate Dx center
    Dx_center = (Dx_f + Dx_b) / 2


    # Generate Dz forward
    diag0 = -np.ones((n_elem, 1))
    diagP = np.ones((n_elem, 1))
    diagBC = np.ones((n_elem, 1))
    diag_all = [
        diagBC.flatten(MATLAB_FLAG),
        diag0.flatten(MATLAB_FLAG),
        diagBC.flatten(MATLAB_FLAG),
    ]
    diag_indexs = [-nx, 0, (n_elem - nx)]
    # Make sparse matrix
    Dz_f = (1 / disc) * sp.diags(diag_all, diag_indexs, shape=(n_elem, n_elem))


    # Generate Dz backward
    diag0 = np.ones((n_elem, 1))
    diagM = -np.ones((n_elem, 1))
    diagBC = -np.ones((n_elem, 1))
    diag_all = [
        diagBC.flatten(MATLAB_FLAG),
        diagM.flatten(MATLAB_FLAG),
        diagBC.flatten(MATLAB_FLAG),
    ]
    diag_indexs = [-nx, 0, (n_elem - nx)]
    # Make sparse matrix
    Dz_b = (1 / disc) * sp.diags(diag_all, diag_indexs, shape=(n_elem, n_elem))

    # Generate Dz squared
    Dz2 = Dz_b @ Dz_f

    # Generate Dz center
    Dz_center = (Dz_f + Dz_b) / 2

    # Take derivatives of permittivity
    # Wait what about PMLs
    print(er[:, -1][:, np.newaxis].shape)
    print(er.shape)
    print(er[:, 0][:, np.newaxis].shape)

    er_extended_z = np.concatenate([er[:, -1][:, np.newaxis], er, er[:, 0][:, np.newaxis]], axis=1)
    dz_er = (er_extended_z[:, 2:] - er_extended_z[:, :-2]) / (2 * disc)
    er_extended_x = np.concatenate([er[0, :][np.newaxis, :], er, er[-1, :][np.newaxis, :]], axis=0)
    dx_er = (er_extended_z[2:, :] - er_extended_z[:-2, :]) / (2 * disc)

    # epsilon = n^2 operator
    n2 = sp.diags(er.flatten(MATLAB_FLAG), 0, shape=(n_elem, n_elem))
    n2_inv = sp.diags(1 / er.flatten(MATLAB_FLAG), 0, shape=(n_elem, n_elem))

    # Identity and zeros operators
    I = sp.eye(n_elem, n_elem)
    Z = 0 * I

    if pol == "TE":
        # TE wave equation
        A = Dx2 + Dz2 + (k0**2) * n2
        B = 1.j * (Dz_b + Dz_f)
        C = -1. * I
    elif pol == "TM":
        # er_inv and n2_inv already computed above
        A = n2 @ Dz_f @ n2_inv + n2 @ Dx_f_plushalf @ n2_inv @ Dx_b + (k0**2) * n2
        B = -1.j * (Dz_b + n2 @ Dz_f @ n2_inv)
        C = -I
    else:
        raise ValueError("Polarization input must be 'TE' or 'TM'")
    
    LH = sp.vstack((
        sp.hstack((A, B)),
        sp.hstack((Z, I)),
    ))
    RH = sp.vstack((
        sp.hstack((Z, -C)),
        sp.hstack((I, Z)),
    ))
    # Solve eigs
    # Phi_out is ( Ey, Ex )(:)
    # TODO: Verify is eigs is the correct functions 
    k_all, Phi_all = spl.eigs(LH, k=num_modes, M=RH, sigma=guessk)
    k_all = np.diag(k_all)

    # Unwrap the field and stuff
    # Reshape and sort the Phis
    # When they come out raw from the modesolver, Phi_all's columns are the
    # eigenvectors
    # The eigenvectors are wrapped by column, then row
    Phi_all = Phi_all[:n_elem, :] # First remove redundant bottom half
    Phi_all = np.reshape(Phi_all, (nx, nz, Phi_all.shape[-1])) # Hopefully this is dimensions x vs. z vs. mode#

    return Phi_all, k_all, A, B
