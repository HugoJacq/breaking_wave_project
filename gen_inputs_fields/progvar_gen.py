import numpy as np
from joblib import Parallel, delayed
from tools import grav, cart2pol


def eta_random(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile):
    """Function to generate a random field given kx-ky spectrum and the x-y-array. Only
    works with uniformly spaced kx-ky so far.
    """

    np.random.seed(0)
    # It is only for focusing case but we don't set tb, xb (time and location of breaking)
    #     tb = 40; xb = 0; yb = 0
    #        phase_tile = -kx_tile*xb-ky_tile*yb+np.random.random_sample(kx_tile.shape)*2*np.pi
    phase_tile = np.random.random_sample(kx_tile.shape) * 2 * np.pi
    eta_tile = np.zeros(x_tile.shape)

    kmod_cart_tile, theta_cart_tile = cart2pol(kx_tile, ky_tile)
    # frequency based on kx or kmod (might need to change for capillary waves)
    omega_tile = (grav * kmod_cart_tile) ** 0.5
    dkx = kx_tile[0, 1] - kx_tile[0, 0]
    dky = ky_tile[1, 0] - ky_tile[0, 0]
    N_grid = x_tile.shape[0]

    if True:  # parallelized version with joblib
        ampl = (2 * F_kxky_tile * dkx * dky) ** 0.5

        def eta_mode(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile):
            a = (
                kx_tile * x_tile[i1, i2]
                + ky_tile * y_tile[i1, i2]
                - omega_tile * t
                + phase_tile
            )
            mode = ampl * np.cos(a)
            return i1, i2, np.sum(mode)

        eta_tile = np.zeros((N_grid, N_grid))
        results = Parallel(n_jobs=-1)(  # -1 means use all available CPU cores
            delayed(eta_mode)(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile)
            for i1 in range(N_grid)
            for i2 in range(N_grid)
        )
        # Assign results back into eta_tile
        for i1, i2, eta_val in results:
            eta_tile[i1, i2] = eta_val

    else:  # Simple summation. To-do: parallelize this
        for i1 in range(0, N_grid):
            for i2 in range(0, N_grid):
                ampl = (2 * F_kxky_tile * dkx * dky) ** 0.5
                a = (
                    (kx_tile * x_tile[i1, i2] + ky_tile * y_tile[i1, i2])
                    - omega_tile * t
                    + phase_tile
                )
                mode = ampl * (np.cos(a))  # uniform spacing in kx and ky
                eta_tile[i1, i2] = np.sum(mode)

    return eta_tile, phase_tile


def gen_velocities_at_z(t, z, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile):
    """
    Function generating the surface elevation and the velocities from
        the spectra for deep water linear waves

        z = 0 at surface, positive upward
    """
    # Common
    np.random.seed(0)
    phase_tile = np.random.random_sample(kx_tile.shape) * 2 * np.pi  # random phase

    kmod_cart_tile, theta_cart_tile = cart2pol(kx_tile, ky_tile)
    # frequency based on kx or kmod (might need to change for capillary waves)
    omega_tile = np.sqrt(grav * kmod_cart_tile)
    dkx = kx_tile[0, 1] - kx_tile[0, 0]
    dky = ky_tile[1, 0] - ky_tile[0, 0]

    N_grid = x_tile.shape[0]  # number of physical grid points

    # initialisation
    u_tile = np.zeros(x_tile.shape)
    v_tile = np.zeros(x_tile.shape)
    w_tile = np.zeros(x_tile.shape)

    # how to compute at one poisition i1, i2
    A = (2 * F_kxky_tile * dkx * dky) ** 0.5  # amplitude of eta
    kmod = np.sqrt(kx_tile**2 + ky_tile**2)  # module of vector k
    B = np.sqrt(grav * kmod)  # coeff for velocities
    z_actual = np.where(z < A, z, A)
    C = np.exp(kmod * z_actual)

    def u_xy(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile):
        a = (
            kx_tile * x_tile[i1, i2]
            + ky_tile * y_tile[i1, i2]
            - omega_tile * t
            + phase_tile
        )
        mode = A * B * C * kx_tile * np.cos(a)
        return i1, i2, np.sum(mode)

    def v_xy(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile):
        a = (
            kx_tile * x_tile[i1, i2]
            + ky_tile * y_tile[i1, i2]
            - omega_tile * t
            + phase_tile
        )
        mode = A * B * C * ky_tile * np.cos(a)
        return i1, i2, np.sum(mode)

    def w_xy(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile):
        a = (
            kx_tile * x_tile[i1, i2]
            + ky_tile * y_tile[i1, i2]
            - omega_tile * t
            + phase_tile
        )
        mode = A * B * C * np.sin(a)
        return i1, i2, np.sum(mode)

    """//compute """

    results = Parallel(n_jobs=-1)(  # -1 means use all available CPU cores
        delayed(u_xy)(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile)
        for i1 in range(N_grid)
        for i2 in range(N_grid)
    )
    for i1, i2, u_val in results:
        u_tile[i1, i2] = u_val

    results = Parallel(n_jobs=-1)(  # -1 means use all available CPU cores
        delayed(v_xy)(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile)
        for i1 in range(N_grid)
        for i2 in range(N_grid)
    )
    for i1, i2, v_val in results:
        v_tile[i1, i2] = v_val

    results = Parallel(n_jobs=-1)(  # -1 means use all available CPU cores
        delayed(w_xy)(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile)
        for i1 in range(N_grid)
        for i2 in range(N_grid)
    )
    for i1, i2, w_val in results:
        w_tile[i1, i2] = w_val

    # print(eta_tile.shape, eta_tile[1, 1], eta_tile.nbytes / 1024**2, "MB")
    # print(u_tile.shape, u_tile[1, 1], u_tile.nbytes / 1024**2, "MB")
    # print(v_tile.shape, v_tile[1, 1])
    # print(w_tile.shape, w_tile[1, 1])

    return u_tile, v_tile, w_tile


def gen_velocities(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile, z_array):
    """Generate 3D velocities from a spectrum, using linear wave theory"""

    Nx = x_tile.shape[1]
    Ny = x_tile.shape[0]
    Nz = z_array.shape[0]
    u = np.zeros((Nz, Ny, Nx))
    v = np.zeros(u.shape)
    w = np.zeros(u.shape)

    for i3, depth in enumerate(z_array):
        u[i3, :, :], v[i3, :, :], w[i3, :, :] = gen_velocities_at_z(
            t, depth, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile
        )

    return u, v, w
