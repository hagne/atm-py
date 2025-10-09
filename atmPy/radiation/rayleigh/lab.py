#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 11:27:53 2022

@author: hagen
"""
import numpy as np

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


def rayleigh_od_first_principles(wl, pressure=1013.25, g=9.80665,
                                      depol_ratio=0.0279, xco2_ppm=420.0):
    """
    Rayleigh optical depth from first principles with adjustable CO2.
    - Adjusts mean molecular mass for a given CO2 mole fraction (ppm).
    - Leaves refractivity n(λ) as in Peck & Reeder (1972); the CO2-induced
      change in n is ~0.01–0.02% and usually negligible for τ.

    Physics
    -------
    σ_R(λ) = [24 π^3 / (λ^4 N_s^2)] * [ (n^2 - 1)^2 / (n^2 + 2)^2 ] * F_K
    τ_R(λ, P) = σ_R(λ) * N_col,
    with N_col ≈ P / (m_air * g) for dry air (hydrostatic + ideal gas).

    Where:
      - λ is vacuum wavelength,
      - n(λ) is refractive index of standard dry air,
      - N_s is Loschmidt’s number at STP,
      - F_K is the King correction factor = (6 + 3δ) / (6 - 7δ),
      - m_air is mean molecular mass per molecule of dry air,
      - g is gravitational acceleration.

    Parameters
    ----------
    wl : float
        Vacuum wavelength in nanonometers (nm). Valid ~200 nm – 2 µm for the
        Peck & Reeder refractivity used here.
    pressure : float
        Surface pressure in hPa. Defaults to 1013.25 hPa.
    g : float
        Gravitational acceleration [m s^-2]. Default 9.80665.
    depol_ratio : float
        Molecular depolarization ratio δ for dry air (weakly λ-dependent).
        Using a common value δ = 0.0279 yields F_K ≈ 1.061.

    Returns
    -------
    numpy.ndarray or float
        Rayleigh optical depth τ_R at each wavelength.

    Notes
    -----
    - Refractivity uses Peck & Reeder (1972) standard dry air (300 ppm CO2).
    - Loschmidt’s number N_s = 2.54743e25 m^-3 (1013.25 hPa, 288.15 K).
    - Mean molecular mass of dry air m_air = 28.9655e-3 kg/mol / N_A.
    - For sub-percent work, consider wavelength-dependent δ(λ) and a modern
      refractivity (e.g., Ciddor 1996) and adjust for CO2/H2O.

    References
    ----------
    - Bucholtz, A. (1995), J. Atmos. Oceanic Technol., 12, 1043–1052.
    - Bodhaine et al. (1999), J. Atmos. Oceanic Technol., 16, 1854–1861.
    - Peck & Reeder (1972), JOSA, 62, 958–962.
    """

    lam_um = np.asarray(wl * 1e-3, dtype=float)
    lam_m = lam_um * 1e-6  # meters

    # --- Refractive index n(λ): Peck & Reeder (1972), λ in µm (unchanged) ---
    sigma_um_inv = 1.0 / lam_um
    sigma2 = sigma_um_inv**2
    n_minus_1_1e8 = (
        8060.51
        + 2480990.0 / (132.274 - sigma2)
        + 17455.7   / (39.32957 - sigma2)
    )
    n = 1.0 + n_minus_1_1e8 * 1e-8

    # --- King factor from depolarization ratio δ ---
    delta = float(depol_ratio)
    F_K = (6.0 + 3.0 * delta) / (6.0 - 7.0 * delta)

    # --- Constants ---
    N_s = 2.54743e25  # Loschmidt number at STP [m^-3]
    N_A = 6.02214076e23  # Avogadro [mol^-1]

    # --- Mean molecular mass with adjustable CO2 (dry air only) ---
    xco2 = float(xco2_ppm) * 1e-6
    # Base dry-air molar fractions (close to US Std Atmosphere), then replace CO2 and
    # take the difference out of N2 (the largest component) to keep sum ≈ 1.
    x_O2 = 0.20946
    x_Ar = 0.00934
    x_N2 = 1.0 - x_O2 - x_Ar - xco2  # ignore other traces at this precision

    M_N2  = 28.0134e-3
    M_O2  = 31.9988e-3
    M_Ar  = 39.9480e-3
    M_CO2 = 44.0095e-3

    M_dry = x_N2*M_N2 + x_O2*M_O2 + x_Ar*M_Ar + xco2*M_CO2  # kg/mol
    m_air = M_dry / N_A                                     # kg per molecule

    # --- Rayleigh cross-section σ_R(λ) [m^2] ---
    n2 = n * n
    term = ((n2 - 1.0) ** 2) / ((n2 + 2.0) ** 2)
    sigma_R = (24.0 * np.pi**3 / (lam_m**4 * N_s**2)) * term * F_K

    # --- Molecular column N_col [m^-2] from hydrostatic balance ---
    P_pa = pressure * 100.0
    N_col = P_pa / (m_air * g)

    tau = sigma_R * N_col
    return tau
