# -*- coding: utf-8 -*-
"""
This package contains routines for calculating atmospheric parameters related to rayleigh scattering.

Most functions are taken from the following reference:

Bucholtz, A. (1995). Rayleigh-scattering calculations for the terrestrial atmosphere. Applied Optics, 34(15), 2765–2773. doi:10.1364/AO.34.002765

@author: hagne, mtat76
"""

import numpy as _np
import pylab as _plt
from scipy import integrate as _integrate

from atmPy.atmosphere import standards as _ats
from atmPy.tools import array_tools as _array_tools


def rayleigh_phase_function(theta, wl):
    """Bucholtz 95 eq(12&13)

    Parameters
    ----------
    theta: float
        Scattering angle in radiant
    wl: float
        Wavelength in nm
    """

    rho = depolarization_factor(wl)  # Depolarization factor at 550 nm
    gamma = rho / (2 - rho)
    P = 3 / (4 * (1 + 2 * gamma)) * ((1 + 3 * gamma) + (1 - gamma) * _np.cos(theta) ** 2)
    return P


def rayleigh_volume_scattering_coeff(P, T, wl):
    """Bucholtz 95 eq(9 + 10)"""
    N_s = 2.54743e19  # cm^{-3}
    beta_s = N_s * scatt_cross(wl) * 1e2  # eqn. 9 unit: 1/m

    P_s = 1013.25  # mbars
    T_s = 288.15  # K
    beta = beta_s * P * T_s / (P_s * T)
    return beta


def scatt_cross(wl):
    """Bucholtz 95 eq(1)

    Arguments
    ---------
    wl in nm

    Returns
    -------
    cross section in cm^2 (watch out mie_scattering calculations are in um^2)"""
    n_s = index_of_refraction(wl)
    N_s = 2.54743e19  # cm^{-3}
    rho = depolarization_factor(wl)
    wl *= 1e-7
    sigma = ((24 * _np.pi ** 3 * (n_s ** 2 - 1) ** 2) / (wl ** 4 * N_s ** 2 * (n_s + 2) ** 2)) * (
    (6 + (3 * rho)) / (6 - (7 * rho)))
    return sigma


def index_of_refraction(wl):
    """Bucholtz 95 eq(4 + 5)
    wl in nm"""

    wl *= 1e-3
    if wl < 0.23:
        n = ((8060.51 + (2480990 / (132.274 - (1 / wl) ** 2)) + (17455.7 / (39.32957 - (1 / wl) ** 2))) * 1e-8) + 1
    else:
        n = (((5791817 / (238.0185 - (1 / wl) ** 2)) + (167909 / (57.362 - (1 / wl) ** 2))) * 1e-8) + 1
    return n


def depolarization_factor(wl):
    """Bucholtz 95 Table 1
    Depolarization factor as fct of wavelength in nm"""
    rho = _np.array([[0.200, 4.545],
                     [0.205, 4.384],
                     [0.210, 4.221],
                     [0.215, 4.113],
                     [0.220, 4.004],
                     [0.225, 3.895],
                     [0.230, 3.785],
                     [0.240, 3.675],
                     [0.250, 3.565],
                     [0.260, 3.455],
                     [0.270, 3.400],
                     [0.280, 3.289],
                     [0.290, 3.233],
                     [0.300, 3.178],
                     [0.310, 3.178],
                     [.32, 3.122],
                     [.33, 3.066],
                     [.34, 3.066],
                     [.35, 3.01],
                     [.37, 3.01],
                     [.38, 2.955],
                     [0.4, 2.955],
                     [0.45, 2.899],
                     [0.5, 2.842],
                     [0.55, 2.842],
                     [0.6, 2.786],
                     [0.75, 2.786],
                     [0.8, 2.73],
                     [1., 2.73]
                     ]).transpose()
    rho[0] *= 1000
    rho[1] *= 0.01
    idx = _array_tools.find_closest(rho[0], wl)
    return rho[1][idx]


def rayleigh_optical_depth(alt, P, T, wl):
    """Bucholtz 95 eq(15). This function describes the optical depth at a particular
    wavelength (wl) due to molecular scattering.

    Parameters
    ----------
    alt : np.array
        Atltitude array in meter.
    P : np.array
        pressure in hPa.
    T : np.array
        Temperature in K.
    wl : float
        wavelength in nm.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    """
    alt in meter
    P in mbar
    T in K
    wl in nm"""

    # when alt has two identical values after another sims will return nan... this checks for it
    dist = alt[1:] - alt[:-1]
    if (dist == 0).sum():
        txt = 'The altitude value has two neiboring elements with the same value, this will case the integration to be nan ... '
        raise ValueError(txt)
    rvsc = rayleigh_volume_scattering_coeff(P, T, wl)
    return _integrate.simps(rvsc, alt)


def rayleigh_angular_scattering_intensity(alt, P, T, wl, theta):
    """This is the integrated intensity of light scattered into the given angle. Note, this is assuming that the
    orientation of the scattered light is vertical. In reality this is only the case if one looks straight up. In most
    cases a slant angle needs to be considered

    Arguments
    ---------
    theta: float or array.
        Scattering angle in radian.

    alt: float or ndarray.
        Altitude in meter.

    P: float or ndarray (same type and shape as alt).
        Pressure in mbar.

    T: float or ndarray (sama type and shape as alt).
        Temperatur in K

    wl: float.
        Wavelength of the considered light in nm


    Retruns
    -------
    ndarray
        this has the unit cm^2 (what out, mie_scattering calculations have um^2
        )
    """

    #     if type(alt).__name__ in ['int','float']:
    #         alt = np.linspace(71000,alt,200)
    if type(alt).__name__ == 'ndarray':
        pass
    else:
        raise TypeError('Sorry only array is allowed currently. "alt" has to be an ndarray')

    rvsc = angular_volume_scattering_coeff(P, T, wl, theta)
    integ = lambda rvsc1D: _integrate.simps(rvsc1D, alt)
    out = _np.apply_along_axis(integ, 0, rvsc)
    #     return integrate.simps(rvsc,alt)
    return out


def angular_volume_scattering_coeff(P, T, wl, theta):
    """Bucholtz 95 eq(14)
    calcultes the volume scattering coeff for all P and Ts. Than it repeats this 1d array by the number of angles theta.
    Result is a 2d array.

    Parameters
    ----------
    P: float or array
        Pressure in hPa
    T: float or array
        Temperatur in K
    wl: float
        Wavelength in nanometer
    theta: float or array
        Angle in radiant

    Returns
    -------
    angular scattering coefficient in 1/m
    """

    if type(theta) != _np.ndarray:
        theta = _np.array([theta], dtype=float)

    if type(P) != _np.ndarray:
        P = _np.array([P], dtype=float)

    if type(T) != _np.ndarray:
        T = _np.array([T], dtype=float)

    b = rayleigh_phase_function(theta, wl)
    a = _np.matrix(_np.repeat(_np.array([rayleigh_volume_scattering_coeff(P, T, wl) / (4 * _np.pi)]), b.shape[0], axis=0)).T
    beta = _np.multiply(a, b)
    return beta


def example():
    """
    wl = 500
    alt = np.linspace(0,71000,20)
    P,T = ats.standard_atmosphere(alt)
    theta = np.linspace(0,np.pi, 50)
    out = rayleigh_angular_scattering_intensity(alt,P,T,wl,theta)
    plt.plot(theta,out)

    rod = rayleigh_optical_depth(alt,P,T,wl)
    print('rod: ', rod)
    rod_alt = integrate.simps(out * np.sin(theta) ,theta) * 2 * np.pi #this is exactly the AOD. the 2pi nesessary because we need to aditionally integrate over phi (from 0 to 2pi), since the function is independend from phi we simply need to multiply by 2pi
    print('rod, when integrating over all scattered light:', rod_alt)"""

    wl = 500
    alt = _np.linspace(0, 71000, 20)
    P, T = _ats.standard_atmosphere(alt)
    theta = _np.linspace(0, _np.pi, 50)
    out = rayleigh_angular_scattering_intensity(alt, P, T, wl, theta)
    _plt.plot(theta, out[0])

    rod = rayleigh_optical_depth(alt, P, T, wl)
    print('rod: ', rod)
    rod_alt = _integrate.simps(out * _np.sin(theta),
                               theta) * 2 * _np.pi  # this is exactly the AOD. the 2pi nesessary because we need to aditionally integrate over phi (from 0 to 2pi), since the function is independend from phi we simply need to multiply by 2pi
    print('rod, when integrating over all scattered light:', rod_alt)
