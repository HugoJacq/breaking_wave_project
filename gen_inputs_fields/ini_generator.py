"""
A generator of initial conditions for basilisk multilayer model

From a spectra, we compute eta (surface elevation) and the ocean currents (3D).

Hypothesis :
    - linear theory
    - modes have random phases


Author : Hugo Jacquet, october 2025
    from a code by Jiarong Wu (@jiarong-wu)
"""

# import pandas as pd
import numpy as np
import toml
from matplotlib import pyplot as plt
from spec_gen import spectrum_PM, spectrum_gen_linear
from progvar_gen import eta_random, gen_eta_velocities_at_z

file_config = "../namlist.toml"

"""
    Reading config file
"""
with open(file_config, "r") as stream:
    config = toml.load(stream)

P = config["P"]
L0 = config["L0"]
kp = config["coeff_kpL0"] * np.pi / L0
N_mode = config["N_mode"]
N_power = config["N_power"]
F_shape = config["F_shape"]  # not used for now
N_grid = config["N_grid"]


def shape(kmod):
    """Choose values here"""
    global P, kp
    F_kmod = spectrum_PM(P=P, kp=kp, kmod=kmod)
    return F_kmod


kmod, F_kmod, kx, ky, F_kxky_tile = spectrum_gen_linear(
    shape, N_mode=N_mode, L=L0, N_power=N_power
)

""" Generate a grid in x-y to visualize random eta """
x = np.linspace(-L0 / 2, L0 / 2, N_grid)
y = np.linspace(-L0 / 2, L0 / 2, N_grid)
x_tile, y_tile = np.meshgrid(x, y)
kx_tile, ky_tile = np.meshgrid(kx, ky)
t = 0
eta_tile, phase_tile = eta_random(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile)
print("kpHs = %g" % (kp * np.std(eta_tile) * 4))
print(eta_tile.shape)

# my test
z = np.zeros(x_tile.shape)
(eta_tile, phase_tile), u_tile, v_tile, w_tile = gen_eta_velocities_at_z(
    t, 0, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile
)
print("kpHs = %g" % (kp * np.std(eta_tile) * 4))
print(eta_tile.shape)


""" Visualization """
plt.imshow(eta_tile, cmap="RdBu_r", vmax=eta_tile.max(), vmin=-eta_tile.max())
plt.colorbar()
plt.show()
