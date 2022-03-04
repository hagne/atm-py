import xarray as xr
import atmPy.precipitation.radar as _radar
import pandas as _pd
import numpy as np
import atmPy.general.timeseries as _timeseries
# def read_ceilometer_nc(fname):
#     ds = xr.open_dataset(fname)
#     ceil = _ceilometry.Ceilometer()
#     ceil.backscatter = _ceilometry.Backscatter(ds.backscatter.to_pandas())
#     return ceil

def read_kazr_nc(fname, timezone = None, keep_xr_dataset = False):
    kazr = _radar.Kazr()
    if type(fname) == str:
        fname = [fname]

    reflects = []
    stds = []
    for fn in fname:
        ds = xr.open_dataset(fn)
        stds.append(ds.snr_copol.to_pandas())
        reflects.append(ds.reflectivity.to_pandas())

    stds = _pd.concat(stds).sort_index()
    reflects =_pd.concat(reflects).sort_index()
    if not isinstance(timezone, type(None)):
        stds.index += np.timedelta64(timezone, 'h')
        reflects.index += np.timedelta64(timezone, 'h')
    kazr.signal2noise_ratio = stds
    kazr.reflectivity = reflects

    if keep_xr_dataset:
        kazr.xr_dataset = ds
    return kazr