import numpy as np
from scipy.interpolate import griddata
from tools import pol2cart

""" 
Some shape of omni-directional spectra
"""


def spectrum_PM(P, kp, kmod):
    """Pierson-Moskowitz spectrum"""
    F_kmod = P * kmod ** (-2.5) * np.exp(-1.25 * (kp / kmod) ** 2)
    return F_kmod


def spectrum_JONSWAP(alpha, kp, kmod):
    """JONSWAP spectrum"""
    F_kmod = alpha * kmod ** (-3) * np.exp(-1.25 * (kp / kmod) ** 2)
    return F_kmod


def spectrum_Gaussian(G, span, kp, kmod):
    """Gaussian spectrum"""
    F_kmod = (G / span) * np.exp(-0.5 * ((kmod - kp) ** 2 / span**2))
    return F_kmod


def spectrum_gen_linear(shape, N_mode=32, N_power=5, L=200):
    """
    Function to generate a kx-ky spectrum based on uni-directional spectrum and a
    cos^N (theta) directional spreading.
    Arguments:
        shape: the spectrum shape (a function)
        N_mode: # of modes (Right now it's N_mode for kx and N_mode+1 for ky;
                has to match what's hard-coded in the spectrum.h headerfile.
        L: physical domain size
    """

    # of modes for the uni-directional spectrum
    #    (doesn't matter as much because of interpolation anyway)
    N_kmod = 64
    N_theta = 64  # Uniform grid in kmod and ktheta, can be finer than N_mode
    thetam = 0  # midline direction
    kmod = np.linspace(2 * np.pi / L, 1.41 * 100 * 2 * np.pi / L, N_kmod)  # Change
    theta = (
        np.linspace(-0.5 * np.pi, 0.5 * np.pi, N_theta) + thetam
    )  # Centered around thetam
    kmod_tile, theta_tile = np.meshgrid(kmod, theta)

    """ Pick the spectrum shape """
    F_kmod = shape(
        kmod
    )  # includes the spectral shape, peak, and energy level of choice
    D_theta = np.abs(np.cos(theta - thetam) ** N_power)
    D_theta = D_theta / np.trapz(D_theta, theta)  # Normalize so the sum equals one
    F_kmod_tile, D_theta_tile = np.meshgrid(F_kmod, D_theta)
    F_kmodtheta_tile = F_kmod_tile * D_theta_tile / kmod_tile  # Notice!! Normalize by k

    """ Uniform grid in kx,ky """
    kx = (
        np.arange(1, N_mode + 1) * 2 * np.pi / L
    )  # based on the grid, interval can't go smaller then pi/L
    ky = np.arange(-N_mode / 2, N_mode / 2 + 1) * 2 * np.pi / L
    kx_tile, ky_tile = np.meshgrid(kx, ky)
    kxp_tile, kyp_tile = pol2cart(kmod_tile, theta_tile)

    """ Project from uniform k to uniform kx,ky """
    F_kxky_tile = griddata(
        (kxp_tile.ravel(), kyp_tile.ravel()),
        F_kmodtheta_tile.ravel(),
        (kx_tile, ky_tile),
        method="linear",
        fill_value=0,
    )

    return kmod, F_kmod, kx, ky, F_kxky_tile
