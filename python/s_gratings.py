import matplotlib.pyplot as plt
import numpy as np

from solve_symm_slab import solve_symm_slab


PI = np.pi


# Initial settings [nm]
disc = 10
lam = 1550
k0 = 2 * PI / lam
index_clad = 1.444#1.5#1.45
index_core = 3.473#2.0#3.47
# index_clad = 1.5#1.5#1.45
# index_core = 1.6#2.0#3.47
y_size = 10000

# Grating dimensions
si_thick = 220
gap_len = 250

# First get slab mode k
neff, beta, kx, _, fields = solve_symm_slab(
    d=si_thick,
    lam0=lam,
    ncl=index_clad,
    nco=index_core,
    pol='TM',
    show_plot=False
)

# Solve for period for desired angle
theta = 15 * PI / 180
period  = (2 * PI) / (beta[0] - k0 * np.sin(theta))

# Draw index
z_coords = np.arange(0, period, step=disc)
x_coords = np.arange(0, y_size, step=disc)
N = index_clad * np.ones((len(x_coords), len(z_coords)))
x_high = x_coords > x_coords[len(x_coords)//2] - si_thick/2 - disc/2 # Waveguide
x_low  = x_coords < x_coords[len(x_coords)//2] + si_thick/2 - disc/2 # Waveguide
z_left = z_coords < period - gap_len - disc/2
mask = np.logical_and(x_high, x_low)[:, np.newaxis] & z_left[np.newaxis, :]
N[mask] = index_core

# Plot index
fig, ax = plt.subplots(1, 1)
extent = [min(x_coords), max(x_coords), min(z_coords), max(z_coords)]
ax.imshow(N.T, extent=extent, origin="lower")
ax.set_xlabel("x")
ax.set_ylabel("z")
ax.set_title("Index distribution")
fig.savefig("test.png")
plt.close(fig)

# Plot field
fig, axes = plt.subplots(2, 3, figsize=(20,13))
field_plot = [
    ["Ex", "Ey", "Ez"],
    ["Hx", "Hy", "Hz"]
]
for ix, ax in enumerate(axes.flatten()):
    ax.plot(fields["x"], np.abs(fields[field_plot[ix // 3][ix % 3]])[0]**2)
    ax.set_title(field_plot[ix // 3][ix % 3])
    ax.set_xlabel("Position [$nm$]")
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Amplitude [$a.u.$]")
fig.suptitle("Mode field")
fig.savefig("field.png")
plt.close(fig)
