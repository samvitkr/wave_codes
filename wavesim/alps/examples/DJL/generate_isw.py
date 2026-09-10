import numpy as np
import sympy as sp

# %%
from PyDJL import DJL, Diagnostic

# Specify the parameters of the problem
A = 22.5  # APE for wave (m^4/s^2)
L, Ly, H = 600.0, 60.0, 150.0  # domain width (m) and depth (m)
NX_l, NZ_l = 32, 256  # grid
NX, NY, NZ = 256, 8, 640

# The unitless density profile (normalized by a reference density rho0)
R = 0.997  # density ratio
z0_d, d_d = 7.0, 0.5
rho_2 = 2 / (1 + R)
a_d = 2 * (1 - R) / (1 + R)
z0_d = z0_d + d_d / 2
Zsym = sp.symbols("z")
# rho_expr  = (-1 - sp.tanh((Zsym + z0_d) / d_d)) * a_d / 2 + 1
rho_expr = 1 - sp.tanh((Zsym + z0_d) / d_d) * a_d / 2
rhoz_expr = sp.diff(rho_expr, Zsym)
# intrho_expr = sp.integrate(rho_expr, Zsym)
print("rho =", rho_expr)
print("rho_z =", rhoz_expr)
# print("int(rho) =", intrho_expr)
rho = sp.lambdify(Zsym, rho_expr, "numpy")
rhoz = sp.lambdify(Zsym, rhoz_expr, "numpy")
# intrho = sp.lambdify(Zsym, intrho_expr, "numpy")
# rho = lambda z: 1 - (a_d / 2) * np.tanh((z + z0_d) / d_d)
# intrho = lambda z: z - a_d * d_d * np.log(np.cosh((z + z0_d) / d_d))
# rhoz = lambda z: -(a_d / d_d) * (1.0 / np.cosh((z + z0_d) / d_d) ** 2)

print(
    "verify density ratio = rho(0) / rho(-H) = {} / {} = {}".format(
        rho(0), rho(-H), rho(0) / rho(-H)
    )
)

# %%
# Create DJL object
print(f"dx = {L / NX_l}, dz = {H / NZ_l}")
# djl = DJL(A, L, H, NX, NZ, rho, rhoz, intrho = intrho)
djl = DJL(A, L, H, NX_l, NZ_l, rho, rhoz)

# Increase the resolution, and iterate to convergence
print(f"dx = {L / NX}, dz = {H / NZ}")
# djl = DJL(A, L, H, NX, NZ, rho, rhoz, intrho = intrho, epsilon=1e-6, initial_guess=djl)
djl = DJL(A, L, H, NX, NZ, rho, rhoz, epsilon=1e-6, initial_guess=djl)

diag = Diagnostic(djl)
print("wavelength = ", diag.wavelength)

# %%
# KdV soliton solution
g = djl.g
sigma = 2 * a_d / (1 + 1 - a_d)
sigma = a_d / (1 - a_d)
b1 = z0_d
b2 = z0_d / 0.05
c0 = np.sqrt(g * sigma * (b1 * b2) / (b1 + b2))
print("c0 = ", c0)
alpha1 = 3 * c0 * (b1 - b2) / (2 * b1 * b2)
alpha2 = (
    3 * c0 / (b1 * b2) ** 2 * (7 * (b1 - b2) ** 2 / 8 - (b1**3 + b2**3) / (b1 + b2))
)
eta0 = -5.0
print("KdV c = ", c0 + eta0 / 3 * alpha1)
print("eKdV c = ", c0 + eta0 / 3 * (alpha1 + 0.5 * alpha2 * eta0))

import h5py
kx0 = 2 * np.pi / L
ky0 = 2 * np.pi / Ly

u = diag.u
ue = diag.shift_grid(u, NX, NZ, symmz = 'even')   # shift u to z endpoints
w = diag.w
we = diag.shift_grid(w, NX, NZ, symmz = 'odd')   # shift w to z endpoints
rho = diag.density
rhoe = diag.shift_grid(rho, NX, NZ, symmz = 'odd')   # shift rho to z endpoints

wout = np.zeros((NZ + 2, NY, NX))
wout[:-1, :, :] = we[:, np.newaxis, :]
uout = np.zeros((NZ + 2, NY, NX))
uout[0, :, :] = ue[0, np.newaxis, :]
uout[1:-1, :, :] = u[:, np.newaxis, :]
uout[-1, :, :] = ue[-1, np.newaxis, :]
cout = np.zeros((NZ + 2, NY, NX))
cout[0, :, :] = rhoe[0, np.newaxis, :]
cout[1:-1, :, :] = rho[:, np.newaxis, :]
cout[-1, :, :] = rhoe[-1, np.newaxis, :]

with h5py.File("grid.h5", "w") as fout:
    zw = np.linspace(0, 1, NZ + 1)
    zw = np.hstack((zw, [0.0]))
    fout["zw"] = zw
    fout["pex"] = kx0
    fout["pey"] = ky0

with h5py.File(f"restart{0:014d}.h5", "w") as fout:
    fout["time"] = 0.0
    fout["pex"] = kx0
    fout["pey"] = ky0
    fout["u"] = uout
    fout["w"] = wout
    fout["v"] = np.zeros_like(uout)
    fout["pp"] = np.zeros_like(uout)
    fout["c0"] = cout


# %%
