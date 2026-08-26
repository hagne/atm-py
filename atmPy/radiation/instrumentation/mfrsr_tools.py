"""
Tools for MFRSR instrumentation data processing and analysis.
"""

import scipy as sp
import pandas as pd
import numpy as np
import xarray as xr
from .spectral import Mfrsr



def check_shadowband_misalignment(data,wl = 870):
    """
    Check for shadowband misalignment in the spectral data.
    Parameters
    ----------
    data : xarray.Dataset
        Spectral data from the MFRSR instrument.
    wl : int, optional
        Wavelength channel to analyze. The default is 870 nm.
    Returns
    -------
    dict
        dict['detected'] : bool
            True if misalignment is detected, False otherwise.
    """
    ds = data
    out = {}
    # do it on single channel
    athis = ds.direct_normal.sel(channel= wl).dropna('datetime')
    
    #remove background to highlight pattern
    athis_rm = athis.rolling(datetime =  11, center = True, ).mean()
    athis = athis - athis_rm
    
    # cut of noise, even if sawtooth peaks are cut of this should still be bettern than having all noise
    lim = 40
    athis = athis.dropna('datetime')
    athis = athis.where(athis < lim).interpolate_na('datetime')
    athis = athis.where(athis > -lim).interpolate_na('datetime')
    
    # get a rolling power spectrum; mor relyable than doing it for the entire day.
    periods = 300
    roll = athis.rolling(datetime = periods, center = True, min_periods=None)
    
    ii = []
    wdata = []
    index = []
    for idx, dat in roll:
        # break
        # pass
        if dat.dropna('datetime').shape[0] != periods:
            continue
        x,y = sp.signal.welch(dat, fs = 1/20, nperseg=periods,)
        wdata.append(y)
        index.append(idx.values)
        ii.append(dat)
    
    wdata = np.array(wdata)
    wdata = pd.DataFrame(wdata, columns=x, index= index)
    wdata.index.name = 'datetime'
    wdata.columns.name = 'freq'
    wdata = xr.DataArray(wdata)
    
    # fit and remove background from power spectrum
    S = wdata.where(wdata.freq > 0, drop=True)
    
    log_f = np.log10(S.freq)
    log_p = np.log10(S.clip(min=np.finfo(float).tiny))
    
    def fit_powerlaw(y, x):
        slope, intercept, *_ = sp.stats.theilslopes(y, x)
        return np.array([slope, intercept])
    
    fit = xr.apply_ufunc(
        fit_powerlaw,
        log_p,
        log_f,
        input_core_dims=[["freq"], ["freq"]],
        output_core_dims=[["parameter"]],
        vectorize=True,
        output_dtypes=[float],
    )
    
    fit = fit.assign_coords(parameter=["slope", "intercept"])
    
    slope = fit.sel(parameter="slope")
    intercept = fit.sel(parameter="intercept")
    
    
    background = 10 ** (
        intercept + slope * np.log10(S.freq)
    )
    
    excess_db = 10 * np.log10(S / background)
    excess = S / background
    
    # detect if powerspectrum shows peaks at expected frequencies and overtone
    # note, this is currently not actually looking for a peak ... potential for improvement if needed
    f0 = 0.009
    tol = 0.001
    
    fundamental = excess.where(
        abs(excess.freq - f0) < tol
    ).max("freq")
    
    harmonic = excess.where(
        abs(excess.freq - 2 * f0) < tol
    ).max("freq")
    
    # determine when peak is present
    detected = (fundamental > 100) & (harmonic > 20)
    out['cleaned_powerspectrum_ts'] = excess
    out['detected_ts'] = detected
    out['detected'] = bool(detected.sum() > 10)
    return out
