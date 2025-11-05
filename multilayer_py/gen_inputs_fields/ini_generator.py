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
import xarray as xr

# from matplotlib import pyplot as plt
from spec_gen import spectrum_PM, spectrum_gen_linear
from progvar_gen import eta_random, gen_velocities

file_config = "../namlist.toml"

"""
    Reading config file
"""
with open(file_config, "r") as stream:
    config = toml.load(stream)

P = config["P"]
L0 = config["L"]
kp = config["coeff_kpL0"] * np.pi / L0
N_mode = config["N_mode"]
N_power = config["N_power"]
F_shape = config["F_shape"]  # not used for now
N_grid = config["N_grid"]
N_layer = config["N_layer"]
H0 = 2 * np.pi / kp  # depth = wavelength, should be enough for deep water assumption

beta=1.29 # for geometric progression of layer height

def shape(kmod):
    """Choose values here"""
    global P, kp
    F_kmod = spectrum_PM(P=P, kp=kp, kmod=kmod)
    return F_kmod


kmod, F_kmod, kx, ky, F_kxky_tile = spectrum_gen_linear(
    shape, N_mode=N_mode, L=L0, N_power=N_power
)


""" Visualization of surface"""
# TO DO


""" Produce eta and 3D velocity fields"""
x = np.linspace(-L0 / 2, L0 / 2, N_grid)
y = np.linspace(-L0 / 2, L0 / 2, N_grid)
z = np.linspace(-H0, 0, N_layer)
# raise Exception("I need to know how it is done in basilisk : before or after remap ?")
x_tile, y_tile = np.meshgrid(x, y)
kx_tile, ky_tile = np.meshgrid(kx, ky)
t = 0
z = np.zeros(x_tile.shape)
eta_tile, phase_tile = eta_random(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile)

h = np.zeros((N_layer, N_grid, N_grid))
H = eta_tile - H0
for i3 in range(N_layer):
    h[i3,:,:] = H*BETA ?
    raise Exception('Im here')



z_array = np.zeros((N_layer, N_grid, N_grid)) + (-H)
for i3 in range(N_layer):
    z_array[i3:,:,:] += 
u, v, w = gen_velocities(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile, z)
""" Show surface fields """

""" Save the file as netcdf"""
# save also input namlist values
attrs = {
    "P": P,
    "L0": L0,
    "kp": kp,
    "N_mode": N_mode,
    "N_power": N_power,
    "F_shape": F_shape,
    "N_grid": N_grid,
    "H": H,
}
coords = {"x": x, "y": y, "z": z, "kx": kx, "ky": ky}
data_vars = {
    "eta": (["y", "x"], eta_tile, {"long_name": "sea surface elevation", "units": "m"}),
    "u": (["z", "y", "x"], u, {"long_name": "Zonal current", "units": "m.s-2"}),
    "v": (["z", "y", "x"], v, {"long_name": "Meridional current", "units": "m.s-2"}),
    "w": (["z", "y", "x"], w, {"long_name": "Vertical current", "units": "m.s-2"}),
    "F_kxky": (
        ["kx", "ky"],
        F_kxky_tile,
        {"long_name": "Initial wave spectrum", "units": "m4"},
    ),
}
ds = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)
ds.to_netcdf(path="./initial.nc", mode="w")
ds.close()
