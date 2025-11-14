import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from multilayer_py.tools.io import read_data
# import from io, diags, plot


filename="/home/jacqhugo/Debut_these/breaking_wave_project/ml_breaking/out.nc"

ds, grid = read_data(filename)





