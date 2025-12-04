'''
This file gather diagnostics for a Basilisk simulation.

References (and associated codes): 

Popinet, S. (2020). A vertically-Lagrangian, non-hydrostatic, multilayer model 
for multiscale free-surface flows. Journal of Computational Physics, 418, 
109609. https://doi.org/10.1016/j.jcp.2020.109609

Wu, J., Popinet, S., & Deike, L. (2023). Breaking wave field statistics with a 
multi-layer model. Journal of Fluid Mechanics, 968, A12. 
https://doi.org/10.1017/jfm.2023.522

Wu, J., Popinet, S., Chapron, B., Farrar, J. T., & Deike, L. (2025). Turbulence
and Energy Dissipation from Wave Breaking. Journal of Physical Oceanography, 
55(9), 1521‑1534. https://doi.org/10.1175/JPO-D-25-0052.1
'''

import numpy as np
import gc


def grad_velocities(ds, grid):
    '''
    Computes the gradient of velocities from a velocity field named u.x u.y u.z
    '''
    delta = ds.x[1]-ds.x[0]
    dudx = grid.interp(grid.diff(ds['u.x'], 'X'), 'X')/delta
    dudy = grid.interp(grid.diff(ds['u.x'], 'Y'), 'Y')/delta
    dudzl = grid.interp(grid.diff(ds['u.x'], 'Z'), 'Z')
    dvdx = grid.interp(grid.diff(ds['u.y'], 'X'), 'X')/delta
    dvdy = grid.interp(grid.diff(ds['u.y'], 'Y'), 'Y')/delta
    dvdzl = grid.interp(grid.diff(ds['u.y'], 'Z'), 'Z')
    dwdx = grid.interp(grid.diff(ds['u.z'], 'X'), 'X')/delta
    dwdy = grid.interp(grid.diff(ds['u.z'], 'Y'), 'Y')/delta
    dwdzl = grid.interp(grid.diff(ds['u.z'], 'Z'), 'Z')
    
    dzdx = grid.interp(grid.diff(ds.z, 'X'), 'X')/delta
    dzdy = grid.interp(grid.diff(ds.z, 'Y'), 'Y')/delta
    dzdzl = grid.interp(grid.diff(ds.z, 'Z'), 'Z')

    ds['dudz'] = (dudzl/dzdzl).compute()
    ds['dudy'] = (dudy - ds['dudz']*dzdy).compute()
    ds['dudx'] = (dudx - ds['dudz']*dzdx).compute()
    ds['dvdz'] = (dvdzl/dzdzl).compute()
    ds['dvdy'] = (dvdy - ds['dvdz']*dzdy).compute()
    ds['dvdx'] = (dvdx - ds['dvdz']*dzdx).compute()
    ds['dwdz'] = (dwdzl/dzdzl).compute()
    ds['dwdy'] = (dwdy - ds['dwdz']*dzdy).compute()
    ds['dwdx'] = (dwdx - ds['dwdz']*dzdx).compute()
    
    print('Field gradient computed!')
    
    del(dudx, dudy, dudzl, dvdx, dvdy, dvdzl, dwdx, dwdy, dwdzl, dzdx, dzdy, dzdzl)
    gc.collect()

    return ds

def vorticity(ds, grid):   
    '''
    Computes vorticity from a velocity field named u.x u.y u.z and a xgcm grid
    '''

    ds = grad_velocities(ds, grid)

    ds['omegaxp'] = ds.dwdy - ds.dvdz
    ds['omegayp'] = ds.dudz - ds.dwdx
    ds['omegazp'] = ds.dvdx - ds.dudy

def enstrophy(u, v, w, dx, dy, dz):
    a=1
