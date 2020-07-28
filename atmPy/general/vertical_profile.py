import warnings as _warnings
import pandas as _pd
import numpy as _np
import matplotlib.pylab as _plt
from copy import deepcopy as _deepcopy
try:
    from netCDF4 import Dataset as _Dataset
except ModuleNotFoundError:
    _warnings.warn('netCDF4 not installed. You might encounter some functionality limitations.')
# from netCDF4 import date2num as _date2num
import atmPy.general.timeseries
from atmPy.tools import plt_tools as _plt_tools
from atmPy.tools import pandas_tools as _pandas_tools
import os as _os
from atmPy.tools import git as _git_tools
from atmPy.tools import array_tools as _array_tools

# _unit_time = 'days since 1900-01-01'

def none2nan(var):
    if type(var).__name__ == 'NoneType':
        var = _np.nan
    return var

def save_netCDF(vp, fname, leave_open = False):

    # if ts._time_format == 'timedelta':
    #     ts.timed

    file_mode = 'w'
    try:
        ni = _Dataset(fname, file_mode)
    except RuntimeError:
        if _os.path.isfile(fname):
            _os.remove(fname)
            ni = _Dataset(fname, file_mode)

    time_dim = ni.createDimension('altitude', vp.data.shape[0])
    dim_data_col = ni.createDimension('data_columns', vp.data.shape[1])

    # ts_time_num = _date2num(ts.data.index.to_pydatetime(), _unit_time)#.astype(float)
    altitude = vp.data.index
    altitude_var = ni.createVariable('altitude', altitude.dtype, 'altitude')
    altitude_var[:] = altitude.values
    altitude_var.units = 'meters'

    var_data = ni.createVariable('data', vp.data.values.dtype, ('altitude', 'data_columns'))
    var_data[:] = vp.data.values

    vp_columns = vp.data.columns.values.astype(str)
    var_data_collumns = ni.createVariable('data_columns', vp_columns.dtype, 'data_columns')
    var_data_collumns[:] = vp_columns

    ni._type = type(vp).__name__
    # ni._data_period = none2nan(vp._data_period)
    ni._x_label = none2nan(vp._x_label)
    ni._y_label =  none2nan(vp._y_label)
    # ni.info = none2nan(vp.info)
    ni._atm_py_commit = _git_tools.current_commit()

    if leave_open:
        return ni
    else:
        ni.close()

def correlate(data,correlant, data_column = False, correlant_column = False, remove_zeros=True, data_lim = None, correlant_lim = None):
    data = data.copy()
    correlant = correlant.copy()

    if not _np.array_equal(data.data.index, correlant.data.index):
        raise ValueError('The indexes of the two columns are not identical, there is no align for vertical profiles yet ... programm an align function')

    if data_column:
        data_values = data.data[data_column].values
    elif data.data.shape[1] > 1:
        raise ValueError('Data contains more than 1 column. Specify which to correlate. Options: %s'%(list(data.data.keys())))
    else:
        data_values = data.data.iloc[:,0].values

    if correlant_column:
        correlant_values = correlant.data[correlant_column].values
    elif correlant.data.shape[1] > 1:
        raise ValueError('''Correlant contains more than 1 column. Specify which to correlate. Options:
%s'''%(list(correlant.data.keys())))
    else:
        correlant_values = correlant.data.iloc[:,0].values


    if data_lim:
        if data_lim[0]:
            data_values[data_values < data_lim[0]] = _np.nan
        if data_lim[1]:
            data_values[data_values > data_lim[1]] = _np.nan

    if correlant_lim:
        if correlant_lim[0]:
            correlant_values[correlant_values < correlant_lim[0]] = _np.nan
        if correlant_lim[1]:
            correlant_values[correlant_values > correlant_lim[1]] = _np.nan

    # import pdb
    # pdb.set_trace()
    out = _array_tools.Correlation(data_values, correlant_values, remove_zeros=remove_zeros, index = data.data.index)
#     out._x_label_orig = 'DataTime'
    return out


class VerticalProfile(object):
    def __init__(self, data):
        data.sort_index(inplace=True)
        data = data[~data.index.duplicated()]
        self.data = data
        self._x_label = None
        self._y_label = 'Altitude'

    ###########################################
    def __sub__(self, other):
        vp = self.copy()
        vp.data = _pd.DataFrame(vp.data.iloc[:, 0] - other.data.iloc[:, 0])
        return vp
    ###########################################

    def align_to(self, ts_other):
        return align_to(self, ts_other)

    def merge(self, ts):
        return merge(self,ts)

    def plot(self, ax=False, **kwargs):
        if not ax:
            f, a = _plt.subplots()
        else:
            a = ax

        for e,k in enumerate(self.data.keys()):
            kwt = kwargs.copy()
            if 'label' not in kwt.keys():
                kwt['label'] = k
            a.plot(self.data[k].values, self.data.index, **kwt)

        if len(self.data.keys()) > 1:
            a.legend(loc = 'best')
        a.set_ylabel(self._y_label)
        a.set_xlabel(self._x_label)
        a.set_ylim((self.data.index.min(), self.data.index.max()))
        return a

    save_netCDF = save_netCDF

    def save(self, fname):
        self.data.to_csv(fname)

    def copy(self):
        return _deepcopy(self)

    def convert2timeseries(self, ts):
        """merges a vertical profile with a timeseries that contains height data
        and returns the a time series where the data of the vertical profile is interpolated
        along the time of the timeseries. ...

        Arguments
        ---------
        ts: timeseries"""
        hk_tmp = ts.convert2verticalprofile()
        data = hk_tmp.data[['TimeUTC']]
        cat_sort_int = _pd.concat([data, self.data]).sort_index().interpolate()
        cat_sort_int = cat_sort_int.dropna()
        cat_sort_int.index = cat_sort_int.TimeUTC
        cat_sort_int = cat_sort_int.drop('TimeUTC', axis=1)
        return atmPy.general.timeseries.TimeSeries(cat_sort_int)

    correlate_to = correlate

    def drop_all_columns_but(self, keep, inplace = False):
        if inplace:
            ts = self
        else:
            ts = self.copy()
        all_keys = ts.data.keys()
        del_keys = all_keys.drop(keep)
        ts.data = ts.data.drop(labels=del_keys, axis=1)
        if inplace:
            return
        else:
            return ts

class VerticalProfile_2D(VerticalProfile):
    def plot(self, xaxis = 0, ax = None, autofmt_xdate = True, cb_kwargs = {}, pc_kwargs = {},  **kwargs):
        if 'cb_kwargs' in kwargs.keys():
            cb_kwargs = kwargs['cb_kwargs']
        if 'pc_kwargs' in kwargs.keys():
            pc_kwargs = pc_kwargs
        f, a, pc, cb = _pandas_tools.plot_dataframe_meshgrid(self.data, xaxis=xaxis, ax=ax, pc_kwargs=pc_kwargs, cb_kwargs=cb_kwargs)
        if autofmt_xdate:
            f.autofmt_xdate()
        return f, a, pc, cb

#### Tools
def merge(ts, ts_other):
    """ Merges current with other timeseries. The returned timeseries has the same time-axes as the current
    one (as opposed to the one merged into it). Missing or offset data points are linearly interpolated.

    Argument
    --------
    ts_orig: the other time series will be merged to this, therefore this timeseries
    will define the time stamps.
    ts: timeseries or one of its subclasses.
        List of TimeSeries objects.

    Returns
    -------
    TimeSeries object or one of its subclasses

    """
    ts_this = ts.copy()
    # ts_data_list = [ts_this.data, ts_other.data]
    # catsortinterp = _pd.concat(ts_data_list).sort_index().interpolate()
    # merged = catsortinterp.groupby(catsortinterp.index).mean().reindex(ts_data_list[0].index)
    # ts_this.data = merged

    ts_data_list = [ts_this.data, ts_other.data]
    bla = _pd.concat(ts_data_list).sort_index()
    catsortinterp = bla.interpolate().where(bla.bfill().notnull())
    merged = catsortinterp.groupby(catsortinterp.index).mean().reindex(ts_data_list[0].index)
    ts_this.data = merged
    return ts_this

def align_to(ts, ts_other):
    """
    Align the TimeSeries ts to another time_series by interpolating (linearly).

    Parameters
    ----------
    ts: original time series
    ts_other: timeseries to align to

    Returns
    -------
    timeseries eqivalent to the original but with an index aligned to the other
    """
    ts = ts.copy()
    ts_other = ts_other.copy()
    ts_other.data = ts_other.data.loc[:,[]]
    ts_t =  merge(ts_other, ts)
    ts.data = ts_t.data
    return ts
