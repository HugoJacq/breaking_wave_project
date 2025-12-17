import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from multilayer_py.tools.io import read_data
from multilayer_py.tools.spec import Spectra2D
from multilayer_py.tools.diags import interpz, vorticity, dissipation
# import from io, diags, plot




### --------------------
#
# General parameters
#
### --------------------

#filename="/home/jacqhugo/Debut_these/breaking_wave_project/ml_breaking/out.nc"
filename="/home/jacqhugo/breaking_wave_project/ml_breaking_save/out.nc"
# getting back Jiarong's data
Jpath = "/home/jacqhugo/breaking_wave_project/multilayer_py/data_Jiarong/"

dpi=200

# TODO: use the namlist to get the parameters

g = 9.81
L0 = 200.0 # m
kp = 10  * np.pi / L0
wp = np.sqrt(g*kp)









# opening file
ds, grid = read_data(filename)



### --------------
### surface elevation
### --------------
fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=dpi)
s=ax.pcolormesh(ds.x, ds.y, ds.eta.isel(time=-1), cmap='Greys', vmin=-2, vmax=2)
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
plt.colorbar(s,ax=ax,orientation='vertical',label=r'$\eta$ (m)')
fig.savefig('eta.png')

print('Depth always under the surface: %f (m)' %ds.z.isel(zl=-1).min().values)

### --------------
### eta spectra
### --------------
# Jiarong's data
Js_0 = np.loadtxt(Jpath+'2023_fig3c_0.txt',skiprows=1,delimiter=",")
Js_20 = np.loadtxt(Jpath+'2023_fig3c_25.txt',skiprows=1,delimiter=",")
Js_100 = np.loadtxt(Jpath+'2023_fig3c_124.txt',skiprows=1,delimiter=",")
Js_120 = np.loadtxt(Jpath+'2023_fig3c_149.txt',skiprows=1,delimiter=",")

Js = [Js_0, Js_20, Js_100, Js_120]

at_t = [0,2, 7] # [0, 20, 100, 120]
colors = plt.get_cmap('plasma')(np.linspace(0, 1, len(at_t)))
# Computing the spectra at some timesteps
s_eta = Spectra2D(ds.eta, compute=False)
# Plotting
fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=dpi)
for k in range(len(at_t)):
    ax.loglog(s_eta.freq_r*2*np.pi/kp, s_eta.sel(time=at_t[k])*kp**3, 
              c=colors[k],
              label=r'$\omega_p t=$'+str(int(wp*at_t[k])))
    ax.set_xlabel(r'$k/k_p$')
    ax.set_ylabel(r'$F_{\eta}(k).k_p^3$')
    ax.set_ylim([1e-7,1e-2])

ax.vlines(1,1e-8,1,ls='--', colors='gray') # TODO: modify this into 1/dx
ax.set_ylim([1e-8, 1e-1])
ax.set_xlim([4e-1,2e1])
ax.legend()
fig.savefig('eta_spectra_evolution.svg')


# vs Jiarong
for k,ttime in enumerate(at_t):
    fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=dpi)
    ax.loglog(s_eta.freq_r*2*np.pi/kp, s_eta.sel(time=ttime)*kp**3, 
              c=colors[k],
              label=r'$\omega_p t=$'+str(int(wp*at_t[k])))
    ax.loglog(Js[k][:,0]/L0/kp,Js[k][:,1], c=colors[k], ls='--')
    ax.set_ylim([1e-7,1e-2])
    fig.savefig(r"eta_spectra_vs_J_%d.svg" % int(wp*at_t[k]))


### --------------
### profiles
### --------------
it = -1
znew = np.arange(-8,0.1,0.1)
z_lagr = ds.z.isel(time=-1).mean(['x','y'])

ux_interp1 = interpz(znew, ds.z.isel(time=-1), ds['u.x'].isel(time=-1), fill_value=np.nan).mean(dim=['x','y'])
ux_interp2 = interpz(znew, ds.z.isel(time=-1), ds['u.x'].isel(time=-1), fill_value=0.).mean(dim=['x','y'])
ux_lagr = ds['u.x'].isel(time=it).mean(['x','y'])

ds = vorticity(ds, grid)
ens_interp1 = interpz(znew, ds.z.isel(time=-1), ds['enstrophy'].isel(time=-1), fill_value=np.nan).mean(dim=['x','y'])
ens_interp2 = interpz(znew, ds.z.isel(time=-1), ds['enstrophy'].isel(time=-1), fill_value=0.).mean(dim=['x','y'])
ens_lagr = ds['enstrophy'].isel(time=it).mean(['x','y'])

ds = dissipation(ds, grid)
diss_interp1 = interpz(znew, ds.z.isel(time=-1), ds['epsilon'].isel(time=-1), fill_value=np.nan).mean(dim=['x','y'])
diss_interp2 = interpz(znew, ds.z.isel(time=-1), ds['epsilon'].isel(time=-1), fill_value=0.).mean(dim=['x','y'])
diss_lagr = ds['epsilon'].isel(time=it).mean(['x','y'])



Jux_interp1 = np.loadtxt(Jpath+'2025_fig7a_abs1.txt',skiprows=1,delimiter=",")
Jux_interp2 = np.loadtxt(Jpath+'2025_fig7a_abs2.txt',skiprows=1,delimiter=",")
Jux_lagr = np.loadtxt(Jpath+'2025_fig7a_layer.txt',skiprows=1,delimiter=",")

Jens_interp1 = np.loadtxt(Jpath+'2025_fig7b_abs1.txt',skiprows=1,delimiter=",")
Jens_interp2 = np.loadtxt(Jpath+'2025_fig7b_abs2.txt',skiprows=1,delimiter=",")
Jens_lagr = np.loadtxt(Jpath+'2025_fig7b_layer.txt',skiprows=1,delimiter=",")

Jdiss_interp1 = np.loadtxt(Jpath+'2025_fig7c_abs1.txt',skiprows=1,delimiter=",")
Jdiss_interp2 = np.loadtxt(Jpath+'2025_fig7c_abs2.txt',skiprows=1,delimiter=",")
Jdiss_lagr = np.loadtxt(Jpath+'2025_fig7c_layer.txt',skiprows=1,delimiter=",")


fig, ax = plt.subplots(1,3,figsize = (9,3),constrained_layout=True,dpi=dpi)
# U
ax[0].set_xlabel('<u> (m/s)')
ax[0].plot(ux_lagr, z_lagr, c='k', ls='-', label='layer')
ax[0].plot(ux_interp1, znew, c='b', ls='-', label='abs 1')
ax[0].plot(ux_interp2, znew, c='b', ls='--', label='abs 2')

# Vorticity
ax[1].set_xlabel('<enstrophy>')
ax[1].semilogx(ens_lagr, z_lagr, c='k', ls='-', label='layer')
ax[1].semilogx(ens_interp1, znew, c='b', ls='-', label='abs 1')
ax[1].semilogx(ens_interp2, znew, c='b', ls='--', label='abs 2')

# Dissipation
ax[2].set_xlabel('<diss>')
ax[2].semilogx(diss_lagr, z_lagr, c='k', ls='-', label='layer')
ax[2].semilogx(diss_interp1, znew, c='b', ls='-', label='abs 1')
ax[2].semilogx(diss_interp2, znew, c='b', ls='--', label='abs 2')

for axe in ax:
    axe.legend(frameon=False)
    axe.set_ylim([-8,0])
    axe.set_ylabel('z (m)')
fig.savefig('avg_profiles.svg')



# Profiles comparison vs Jiarong's papers
fig, ax = plt.subplots(1,3,figsize = (9,3),constrained_layout=True,dpi=dpi)
# U
ax[0].set_xlabel('<u> (m/s)')
ax[0].plot(ux_lagr, z_lagr, c='k', ls='-', label='layer')
ax[0].plot(ux_interp1, znew, c='b', ls='-', label='abs 1')
ax[0].plot(ux_interp2, znew, c='b', ls='--', label='abs 2')
# Vorticity
ax[1].set_xlabel('<enstrophy>')
ax[1].semilogx(ens_lagr, z_lagr, c='k', ls='-', label='layer')
ax[1].semilogx(ens_interp1, znew, c='b', ls='-', label='abs 1')
ax[1].semilogx(ens_interp2, znew, c='b', ls='--', label='abs 2')
# Dissipation
ax[2].set_xlabel('<diss>')
ax[2].semilogx(diss_lagr, z_lagr, c='k', ls='-', label='layer')
ax[2].semilogx(diss_interp1, znew, c='b', ls='-', label='abs 1')
ax[2].semilogx(diss_interp2, znew, c='b', ls='--', label='abs 2')



#plt.show()

