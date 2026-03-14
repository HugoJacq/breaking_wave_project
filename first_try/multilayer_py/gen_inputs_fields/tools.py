"""
Tool box

Author : Hugo Jacquet, october 2025
    based on tools by Jiarong Wu (@jiarong-wu)

"""

import numpy as np

grav = 9.81  # m/s2


# Function to convert polar to cartesian and interpolate
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)
