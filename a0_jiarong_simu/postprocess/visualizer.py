"""
# Description

    This script compares a Basilisk multi-layer simulation to [1,2].

[1] Wu, J., Popinet, S., & Deike, L. (2023). Breaking wave field statistics 
with a multi-layer model. Journal of Fluid Mechanics, 968, A12. 
https://doi.org/10.1017/jfm.2023.522

[2] Wu, J., Popinet, S., Chapron, B., Farrar, J. T., & Deike, L. (2025). 
Turbulence and Energy Dissipation from Wave Breaking. Journal of Physical 
Oceanography, 55(9), 1521–1534. https://doi.org/10.1175/JPO-D-25-0052.1

# Conda Environment

    conda activate ml_dec2025

# How to run

    python3 visualizer.py
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os.path
import time

from data_reader import read_data, build_grid
from tools import Spectra2D
from diags import interpz, grad_velocities, vorticity, dissipation
# import from io, diags, plot




### --------------------
#
# General parameters
#
### --------------------

# My data
filename="/home/jacqhugo/basilisk/wiki/sandbox/hugoj/reproducing_jiarongs_plots/N10_P0.02/out.nc"
# getting back Jiarong's data
Jpath = "data_Jiarong/"
save_nc = './data.nc'

dpi=200

g = 9.81
L0 = 200.0 # m
kp = 10  * np.pi / L0
wp = np.sqrt(g*kp)
fp = wp/(2*np.pi)
Tp = 1/fp

print("--------------------")
print("Peak values:")
print("kp",kp)
print("wp",wp)
print("fp",fp)
print("TP",Tp)
print("--------------------\n")


### --------------
### surface elevation
### --------------
if True:
    print("* Surface elevation, spectra")
    # opening file
    ds, grid = read_data(filename)
    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=dpi)
    s=ax.pcolormesh(ds.x, ds.y, ds.eta.isel(time=-1), cmap='Greys_r', vmin=-1, vmax=1)
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

    at_t = [0, 20, 100, 120]
    colors = plt.get_cmap('plasma')(np.linspace(0, 1, len(at_t)))
    # Computing the spectra at some timesteps
    s_eta = Spectra2D(ds.eta, compute=False)
    # Plotting
    fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=dpi)
    for k in range(len(at_t)):
        ax.loglog(s_eta.freq_r*2*np.pi*L0, s_eta.sel(time=at_t[k])*kp**3,
        #ax.loglog(s_eta.freq_r, s_eta.sel(time=at_t[k])*kp**3,
                    c=colors[k],
                #label=r'$t/T_p=$'+str(int(at_t[k]/Tp)))
                label=r'$w_p t=$'+str(int(wp*at_t[k])))
        
    ax.set_xlabel(r'$kL$')
    #ax.set_xlabel(r'$k$')
    ax.set_ylabel(r'$F_{\eta}(k).k_p^3$')
    ax.vlines(1,1e-8,1,ls='--', colors='gray') # TODO: modify this into 1/dx
    ax.vlines(kp*L0,1e-8,1,ls='--', colors='gray')
    ax.set_ylim([1e-7, 1e-2])
    # ax.set_xlim([4e-1,2e1])
    ax.set_xlim([1e1,1e3])
    ax.legend()
    fig.savefig('eta_spectra_evolution.svg')


    # vs Jiarong
    for k,ttime in enumerate(at_t):
        fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=dpi)
        ax.loglog(s_eta.freq_r*2*np.pi*L0, s_eta.sel(time=ttime)*kp**3, 
                c=colors[k],
                #label=r'$t/Tp=$'+str(int(at_t[k]/Tp)))
                label=r'$w_p t=$'+str(int(wp*at_t[k])))
        ax.loglog(Js[k][:,0],Js[k][:,1], c=colors[k], ls='--')
        ax.set_ylim([1e-7,1e-2])
        ax.set_xlim([1e1,1e3])
        ax.set_xlabel(r'$kL$')
        ax.set_ylabel(r'$F_{\eta}(k).k_p^3$')
        ax.legend()
        fig.savefig(r"eta_spectra_vs_J_%d.svg" % int(wp*at_t[k]))

    # if not os.path.isfile(save_nc):
    #     ds.close()
    #     ds.to_netcdf(save_nc)
    

### --------------
### profiles
### --------------

if False:
    print("* Profiles")
        
    print('Computing profiles ...')



    it = -1
    znew = np.arange(-8,0.1,0.1)
    
    
    # gradients
    if not os.path.isfile(save_nc):
        print('I need to compute diags ...')
        ds, grid = read_data(filename)
        
        ds, update = grad_velocities(ds, grid)
        ds, update = vorticity(ds, grid)
        ds, update = dissipation(ds, grid)
        if update:
            print(f'saving diags to {save_nc}')
            ds.to_netcdf(save_nc)
    else:
        ds = xr.open_dataset(save_nc)
        grid = build_grid(ds)


    
    
    ds = ds.chunk({'time':5, 'x':128, 'y':128, 'zl':-1})

    # print('starting interp1')
    # t1 = time.time()
    # print(ds.z.isel(time=-1))
    # print(ds['u.x'].isel(time=-1))
    # u_interp = interpz(ds.z.isel(time=-1), ds['u.x'].isel(time=-1), znew).compute()
    # u_interp = u_interp.assign_coords(znew=znew) 
    # t2 = time.time()
    # print(f'interp1 done, {int(t2-t1)} s')

    z_lagr = ds.z.isel(time=-1).mean(['x','y'])
    ux_interp1 = interpz(ds.z.isel(time=-1), ds['u.x'].isel(time=-1), znew, fill_value=np.nan).mean(dim=['x','y']).compute()
    ux_interp2 = interpz(ds.z.isel(time=-1), ds['u.x'].isel(time=-1), znew, fill_value=0.).mean(dim=['x','y']).compute()
    ux_lagr = ds['u.x'].isel(time=it).mean(['x','y'])
    
    ens_interp1 = interpz(ds.z.isel(time=-1), ds['enstrophy'].isel(time=-1), znew, fill_value=np.nan).mean(dim=['x','y']).compute()
    ens_interp2 = interpz(ds.z.isel(time=-1), ds['enstrophy'].isel(time=-1), znew, fill_value=0.).mean(dim=['x','y']).compute()
    ens_lagr = ds['enstrophy'].isel(time=it).mean(['x','y'])

    diss_interp1 = interpz(ds.z.isel(time=-1), ds['epsilon'].isel(time=-1), znew, fill_value=np.nan).mean(dim=['x','y']).compute()
    diss_interp2 = interpz(ds.z.isel(time=-1), ds['epsilon'].isel(time=-1), znew, fill_value=0.).mean(dim=['x','y']).compute()
    diss_lagr = ds['epsilon'].isel(time=it).mean(['x','y'])

    # z_lagr = ds.z.isel(time=-1).mean(['x','y'])
    # ux_interp1 = interpz(znew, ds.z.isel(time=-1), ds['u.x'].isel(time=-1), fill_value=np.nan).mean(dim=['x','y']).compute()
    # ux_interp2 = interpz(znew, ds.z.isel(time=-1), ds['u.x'].isel(time=-1), fill_value=0.).mean(dim=['x','y'])
    # ux_lagr = ds['u.x'].isel(time=it).mean(['x','y'])
    #
    # ens_interp1 = interpz(znew, ds.z.isel(time=-1), ds['enstrophy'].isel(time=-1), fill_value=np.nan).mean(dim=['x','y'])
    # ens_interp2 = interpz(znew, ds.z.isel(time=-1), ds['enstrophy'].isel(time=-1), fill_value=0.).mean(dim=['x','y'])
    # ens_lagr = ds['enstrophy'].isel(time=it).mean(['x','y'])
    #
    # diss_interp1 = interpz(znew, ds.z.isel(time=-1), ds['epsilon'].isel(time=-1), fill_value=np.nan).mean(dim=['x','y'])
    # diss_interp2 = interpz(znew, ds.z.isel(time=-1), ds['epsilon'].isel(time=-1), fill_value=0.).mean(dim=['x','y'])
    # diss_lagr = ds['epsilon'].isel(time=it).mean(['x','y'])

    print(' Done !')


    Jux_interp1 = np.loadtxt(Jpath+'2025_fig7a_abs1.txt',skiprows=1,delimiter=",")
    Jux_interp2 = np.loadtxt(Jpath+'2025_fig7a_abs2.txt',skiprows=1,delimiter=",")
    Jux_lagr = np.loadtxt(Jpath+'2025_fig7a_layer.txt',skiprows=1,delimiter=",")
    
    Jens_interp1 = np.loadtxt(Jpath+'2025_fig7b_abs1.txt',skiprows=1,delimiter=",")
    Jens_interp2 = np.loadtxt(Jpath+'2025_fig7b_abs2.txt',skiprows=1,delimiter=",")
    Jens_lagr = np.loadtxt(Jpath+'2025_fig7b_layer.txt',skiprows=1,delimiter=",")

    Jdiss_interp1 = np.loadtxt(Jpath+'2025_fig7c_abs1.txt',skiprows=1,delimiter=",")
    Jdiss_interp2 = np.loadtxt(Jpath+'2025_fig7c_abs2.txt',skiprows=1,delimiter=",")
    Jdiss_lagr = np.loadtxt(Jpath+'2025_fig7c_layer.txt',skiprows=1,delimiter=",")
    
    if True:
        # Profiles comparison vs Jiarong's papers
        fig, ax = plt.subplots(1,3,figsize = (9,3),constrained_layout=True,dpi=dpi)
        # U
        ax[0].set_xlabel('<u> (m/s)')
        ax[0].plot(ux_lagr, z_lagr, c='k', ls='-', label='layer',  marker='s', markerfacecolor='None')
        ax[0].plot(ux_interp1, znew, c='b', ls='-', label='abs 1')
        ax[0].plot(ux_interp2, znew, c='b', ls='--', label='abs 2')
        ax[0].set_xlim([0,0.2])

        # Vorticity
        ax[1].set_xlabel('<enstrophy>')
        ax[1].semilogx(ens_lagr, z_lagr, c='k', ls='-', label='layer',  marker='s', markerfacecolor='None')
        ax[1].semilogx(ens_interp1, znew, c='b', ls='-', label='abs 1')
        ax[1].semilogx(ens_interp2, znew, c='b', ls='--', label='abs 2')
        ax[1].set_xlim([1e-4,1])

        # Dissipation
        ax[2].set_xlabel('<diss>')
        ax[2].semilogx(diss_lagr, z_lagr, c='k', ls='-', label='layer', marker='s', markerfacecolor='None')
        ax[2].semilogx(diss_interp1, znew, c='b', ls='-', label='abs 1')
        ax[2].semilogx(diss_interp2, znew, c='b', ls='--', label='abs 2')
        ax[2].set_xlim([1e-3,10])

        for axe in ax:
            axe.legend(frameon=False)
            axe.set_ylim([-10,0])
            axe.set_ylabel('z (m)')
        fig.savefig('avg_profiles.svg')


    if True:
        # Profiles comparison vs Jiarong's papers
        fig, ax = plt.subplots(1,3,figsize = (9,3),constrained_layout=True,dpi=dpi)
        # U
        ax[0].set_xlabel('<u> (m/s)')
        ax[0].plot(Jux_lagr[:,0], Jux_lagr[:,1], c='gray', label='J layer', marker='x', markerfacecolor='None')
        ax[0].plot(ux_lagr, z_lagr, c='k', ls='-', label='layer', marker='s', markerfacecolor='None')
        ax[0].set_xlim([0,0.2])
        # ax[0].plot(ux_interp1, znew, c='b', ls='-', label='abs 1')
        # ax[0].plot(ux_interp2, znew, c='b', ls='--', label='abs 2')
        # Vorticity
        ax[1].set_xlabel('<enstrophy>')
        ax[1].semilogx(Jens_lagr[:,0], Jens_lagr[:,1], c='gray', label='J layer', marker='x', markerfacecolor='None')
        ax[1].semilogx(ens_lagr, z_lagr, c='k', ls='-', label='layer', marker='s', markerfacecolor='None')
        ax[1].set_xlim([1e-4,1])

        # ax[1].semilogx(ens_interp1, znew, c='b', ls='-', label='abs 1')
        # ax[1].semilogx(ens_interp2, znew, c='b', ls='--', label='abs 2')
        # Dissipation
        ax[2].set_xlabel('<diss>')
        ax[2].semilogx(Jdiss_lagr[:,0], Jdiss_lagr[:,1], c='gray', label='J layer', marker='x', markerfacecolor='None' )
        ax[2].semilogx(diss_lagr, z_lagr, c='k', ls='-', label='layer', marker='s', markerfacecolor='None')
        # ax[2].semilogx(diss_interp1, znew, c='b', ls='-', label='abs 1')
        # ax[2].semilogx(diss_interp2, znew, c='b', ls='--', label='abs 2')
        ax[2].set_xlim([1e-3,10])

        for axe in ax:
            axe.legend(frameon=False)
            axe.set_ylim([-10,0])
            axe.set_ylabel('z (m)')
        fig.savefig('avg_profiles_vsJiarong.png')


plt.show()

