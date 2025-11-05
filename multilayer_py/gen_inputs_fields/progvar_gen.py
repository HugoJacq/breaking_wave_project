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

    else:  
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

def u_x(t, x, y, z, F_kxky_tile, kx_tile, ky_tile, omega_tile, phase_tile):
    u_x = 0.
    dkx = kx_tile[0, 1] - kx_tile[0, 0]
    dky = ky_tile[1, 0] - ky_tile[0, 0]
    kmod = np.sqrt(kx_tile**2 + ky_tile**2)     # module of vector k
    A = (2 * F_kxky_tile * dkx * dky) ** 0.5    # amplitude of eta
    B = np.sqrt(grav * kmod)                    # coeff for velocities
    theta = np.arctan(ky_tile/kx_tile)
    # loop on modes
    for j in range(kx_tile.shape[0]):
        for i in range(kx_tile.shape[1]):
            z_actual = np.where(z < A[j,i], z, A[j,i])
            C = np.exp(kmod[j,i] * z_actual)
            a = (kx_tile * x
            + ky_tile * y
            - omega_tile * t
            + phase_tile)
            u_x += A[j,i] * B[j,i] * C * np.cos(theta) * np.cos(a)
    return u_x

def u_y(t, x, y, z, F_kxky_tile, kx_tile, ky_tile, omega_tile, phase_tile):
    u_y = 0.
    dkx = kx_tile[0, 1] - kx_tile[0, 0]
    dky = ky_tile[1, 0] - ky_tile[0, 0]
    kmod = np.sqrt(kx_tile**2 + ky_tile**2)     # module of vector k
    A = (2 * F_kxky_tile * dkx * dky) ** 0.5    # amplitude of eta
    B = np.sqrt(grav * kmod)                    # coeff for velocities
    theta = np.arctan(ky_tile/kx_tile)
    # loop on modes
    for j in range(kx_tile.shape[0]):
        for i in range(kx_tile.shape[1]):
            z_actual = np.where(z < A[j,i], z, A[j,i])
            C = np.exp(kmod[j,i] * z_actual)
            a = (kx_tile * x
            + ky_tile * y
            - omega_tile * t
            + phase_tile)
            u_y += A[j,i] * B[j,i] * C * np.sin(theta) * np.cos(a)
    return u_y

def u_z(t, x, y, z, F_kxky_tile, kx_tile, ky_tile, omega_tile, phase_tile):
    u_z = 0.
    dkx = kx_tile[0, 1] - kx_tile[0, 0]
    dky = ky_tile[1, 0] - ky_tile[0, 0]
    kmod = np.sqrt(kx_tile**2 + ky_tile**2)     # module of vector k
    A = (2 * F_kxky_tile * dkx * dky) ** 0.5    # amplitude of eta
    B = np.sqrt(grav * kmod)                    # coeff for velocities
    # loop on modes
    for j in range(kx_tile.shape[0]):
        for i in range(kx_tile.shape[1]):
            z_actual = np.where(z < A[j,i], z, A[j,i])
            C = np.exp(kmod[j,i] * z_actual)
            a = (kx_tile * x
            + ky_tile * y
            - omega_tile * t
            + phase_tile)
            u_y += A[j,i] * B[j,i] * C  * np.sin(a)
    return u_y




def gen_velocities_at_z(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile, z_tile):
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
    #
    # def u_xy_old(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile):
    #     for j in range(kx_tile.shape[0]):
    #         for i in range(kx_tile.shape[1]):
    #             z_actual = np.where(z_tile[i1,i2] < A[j,i], z_tile[i1,i2], A[j,i])
    #             C = np.exp(kmod[j,i] * z_actual)
    #
    #             a = (
    #                 kx_tile * x_tile[i1, i2]
    #                 + ky_tile * y_tile[i1, i2]
    #                 - omega_tile * t
    #                 + phase_tile
    #             )
    #             mode = A[j,i] * B[j,i] * C * kx_tile * np.cos(a)
    #     return i1, i2, np.sum(mode)

    def u_xy(i1, i2):
        # wrapper for // compute on i1, i2
        ux = u_x(t, x_tile[i1,i2], y_tile[i1,i2], z_tile[i1,i2], F_kxky_tile, kx_tile, ky_tile, omega_tile, phase_tile)
        return i1, i2, ux
    
    def v_xy(i1, i2):
        # wrapper for // compute on i1, i2
        v = u_y(t, x_tile[i1,i2], y_tile[i1,i2], z_tile[i1,i2], F_kxky_tile, kx_tile, ky_tile, omega_tile, phase_tile)
        return i1, i2, v

    def w_xy(i1, i2):
        # wrapper for // compute on i1, i2
        w = u_z(t, x_tile[i1,i2], y_tile[i1,i2], z_tile[i1,i2], F_kxky_tile, kx_tile, ky_tile, omega_tile, phase_tile)
        return i1, i2, w


    # def v_xy(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile):
    #     a = (
    #         kx_tile * x_tile[i1, i2]
    #         + ky_tile * y_tile[i1, i2]
    #         - omega_tile * t
    #         + phase_tile
    #     )
    #     mode = A * B * C * ky_tile * np.cos(a)
    #     return i1, i2, np.sum(mode)

    # def w_xy(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile):
    #     a = (
    #         kx_tile * x_tile[i1, i2]
    #         + ky_tile * y_tile[i1, i2]
    #         - omega_tile * t
    #         + phase_tile
    #     )
    #     mode = A * B * C * np.sin(a)
    #     return i1, i2, np.sum(mode)

    """//compute """
    print('     u')
    results = Parallel(n_jobs=-1)(  # -1 means use all available CPU cores
        #delayed(u_xy)(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile)
        delayed(u_xy)(i1, i2)
        for i1 in range(N_grid)
        for i2 in range(N_grid))
    for i1, i2, u_val in results:
        u_tile[i1, i2] = u_val
    raise Exception

    print('     v')
    results = Parallel(n_jobs=-1)(  # -1 means use all available CPU cores
        #delayed(v_xy)(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile)
        delayed(v_xy)(i1, i2)
        for i1 in range(N_grid)
        for i2 in range(N_grid))
    for i1, i2, v_val in results:
        v_tile[i1, i2] = v_val
    
    print('     w')
    results = Parallel(n_jobs=-1)(  # -1 means use all available CPU cores
        #delayed(w_xy)(t, i1, i2, x_tile, y_tile, kx_tile, ky_tile)
        delayed(w_xy)(i1, i2)
        for i1 in range(N_grid)
        for i2 in range(N_grid))
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

    for i3 in range(1): #Nx):
        print("level={%d}",i3)
        u[i3, :, :], v[i3, :, :], w[i3, :, :] = gen_velocities_at_z(
            t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile, z_array[i3]
        )

    return u, v, w
