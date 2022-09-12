#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 11:27:53 2022

@author: hagen
"""
# import numpy as np

def rayleigh_od_johnsmethod(pressure_surface, wavelength, pressure_standard = 1013.25):
    """
    This is the method to estimate total collum OD from rayleigh scattering 
    that John is using in his AOD analysis (aod_analysis.f)

    Parameters
    ----------
    pressure_surface : TYPE
        DESCRIPTION.
    wavelength : TYPE
        DESCRIPTION.
    pressure_standard : TYPE, optional
        DESCRIPTION. The default is 1013.25.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    p_s = pressure_surface
    p_o = pressure_standard
    # wave_l = 500
    wave_l = wavelength
    exponent = -4.15 + (0.2*(wave_l/1000.))
    exp_term = (wave_l/1000.)**exponent
    rayleigh = (0.0088*exp_term*p_s)/p_o
    # trans = np.e**(-rayleigh)
    return rayleigh