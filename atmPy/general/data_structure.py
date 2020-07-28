import pandas as _pd
import numpy as _np
import warnings as _warnings
try:
    from statsmodels import robust as _robust
except ModuleNotFoundError:
    _warnings.warn('statsmodels not installed. You might encounter some functionality limitations.')
from functools import wraps as _wraps

def close_gaps(ts, verbose = False):
    """This is an older version to deal with gaps ... rather consider using the ones below"""
    ts = ts.copy()
    ts.data = ts.data.sort_index()
    if type(ts.data).__name__ == 'Panel':
        data = ts.data.items.values
        index = ts.data.items
    else:
        data = ts.data.index.values
        index = ts.data.index
    index_df = _pd.DataFrame(index = index)

    dt = data[1:] - data[:-1]
    dt = dt / _np.timedelta64(1,'s')

    median = _np.median(dt)
    if median > (1.1 * ts._data_period) or median < (0.9 * ts._data_period):
        _warnings.warn('There is a periode and median missmatch (%0.1f,%0.1f), this is either due to an error in the assumed period or becuase there are too many gaps in the _timeseries.'%(median,ts._data_period))

    point_dist = (index.values[1:] - index.values[:-1]) / _np.timedelta64(1, 's')
    where = point_dist > 2 * ts._data_period
    off_periods = _np.array([index[:-1][where], index[1:][where]]).transpose()
    if verbose:
        print('found %i gaps'%(off_periods.shape[0]))
    for i, op in enumerate(off_periods):
        no_periods = round((op[1] - op[0])/ _np.timedelta64(1,'s')) / ts._data_period
        out = _pd.date_range(start = op[0], periods= no_periods, freq= '%i s'%ts._data_period)
        out = out[1:]
        out = _pd.DataFrame(index = out)
        index_df = _pd.concat([index_df, out])
    index_df.sort_index(inplace=True)
    ts.data = ts.data.reindex(index_df.index)
    return ts

def detect_gaps(ts, toleranz=1.95, return_all=False):
    if type(ts.data).__name__ == "Panel":
        idx = ts.data.items
    else:
        idx = ts.data.index
    dt = (idx[1:] - idx[:-1]) / _np.timedelta64(1, 's')

    med = _np.median(dt)
    mad = _robust.mad(dt)

    if mad == 0:
        noofgaps = 0
        dt = 0
        period = int(med)
    else:
        hist, edges = _np.histogram(dt[_np.logical_and((med - mad) < dt, dt < (med + mad))], bins=100)
        period = int(round((edges[hist.argmax()] + edges[hist.argmax() + 1]) / 2))
        noofgaps = dt[dt > toleranz * period].shape[0]

    if return_all:
            return {'index':idx,
                    'number of gaps':noofgaps,
                    'dt': dt,
                    'period (s)': period}
    else:
        return noofgaps

def fill_gaps_with(ts, what=0, toleranz=1.95, inplace = True):
    # if type(ts).__name__ == 'DataStructure':
    #     ts = ts.parent_ts
    gaps = ts.data_structure.detect_gaps(toleranz=toleranz, return_all=True)
    idx = gaps['index']
    # noofgaps = gaps['number of gaps']
    dt = gaps['dt']
    # period = gaps['period (s)']
    period = ts._data_period
    # print(type(toleranz),toleranz, type(period), period)
    gap_data_list = [ts.data]
    for idxf, idxn, dtt in zip(idx[:-1][dt > toleranz * period], idx[1:][dt > toleranz * period],
                               dt[dt > toleranz * period]):
        no2add = int(round(((idxn - idxf) / _np.timedelta64(1, 's')) / period)) - 1


        # print(no2add)
        newidx = _pd.to_datetime((_np.arange(no2add) + 1) * _np.timedelta64(period, 's') + idxf.to_datetime64())
        dft = _pd.DataFrame(index=newidx, columns=ts.data.columns)
        if what != _np.nan:
            dft[:] = what
        gap_data_list.append(dft)


        # return ts, no2add, idxf, period, what
        # for i in range(no2add):
        #     newidx = idxf + _np.timedelta64((i + 1) * period, 's')
        #     ts.data.loc[newidx, :] = what
    new_data = _pd.concat(gap_data_list)
    new_data.sort_index(inplace=True)
    if inplace:
        ts.data = new_data
        # ts.data.sort_index(inplace=True)
        return
    else:
        tst = ts.copy()
        tst.data = new_data
        return tst

def estimate_sampling_period(ts, toleranz=1.95):
    gaps = ts.data_structure.detect_gaps(toleranz=toleranz, return_all=True)
    return gaps['period (s)']

def _adjusttype(fkt):
    @_wraps(fkt)
    def wrapper(ts, *args, **kwargs):
        if type(ts).__name__ == 'DataStructure':
            ts = ts.parent_ts
        return fkt(ts, *args, **kwargs)
    return wrapper


class DataStructure(object):
    def __init__(self, ts):
        self.parent_ts = ts

    detect_gaps = _adjusttype(detect_gaps)
    fill_gaps_with = _adjusttype(fill_gaps_with)
    estimate_sampling_period = _adjusttype(estimate_sampling_period)
