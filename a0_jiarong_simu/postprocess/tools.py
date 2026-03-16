import xarray as xr
import xrft

def Spectra2D(eta: xr.DataArray, compute=False):

    window = 'hann'
    nfactor = 4
    truncate = True
    detrend = None
    window_correction = True

    s = xrft.isotropic_power_spectrum(eta, 
                                dim=('x','y'), 
                                window=window,
                                nfactor=nfactor, 
                                truncate=truncate, 
                                detrend=detrend,
                                window_correction=window_correction)
    return s
