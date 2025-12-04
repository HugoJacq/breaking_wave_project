''' This file gather diagnostics for a Basilisk simulation, With a focus on
spectral analysis.

References:

Lecture F. Ardhuin 2025 "Waves in geoscience", 
https://github.com/ardhuin/waves in geosciences

'''
import xarray as xr


def Spectra2D(eta: xr.DataArray, compute=False):

    window = 'hann'
    nfactor = 4
    truncate = True
    detrend = None
    window_correction = True

    s = xrft.isotropic_power_spectrum(ds.z.isel(zl=-1).drop(['z']), 
                                dim=('x','y'), 
                                window=window,
                                nfactor=nfactor, 
                                truncate=truncate, 
                                detrend=detrend,
                                window_correction=window_correction)


    return s




def multi_Spectra2D(A: np.ndarray, dx: float, dy: float, N: int): 
    ''' 
    This function computes a 2D spectra on N x N tiles of the field A, then 
    average the result to give a less spiky spetra.
    Also the 95% interval confidence is given as E1 and E2.
    '''

    a=1
