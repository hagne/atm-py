__author__ = 'htelg'

from copy import deepcopy as _deepcopy
import atmPy.general.vertical_profile 

import pandas as _pd
import numpy as _np
import matplotlib.pylab as _plt

from atmPy.tools import pandas_tools as _pandas_tools
from atmPy.tools import time_tools as _time_tools
from atmPy.tools import array_tools as _array_tools
from atmPy.tools import plt_tools as _plt_tools
from netCDF4 import Dataset as _Dataset
from netCDF4 import num2date as _num2date
from netCDF4 import date2num as _date2num

from atmPy.tools import git as _git_tools

import warnings as _warnings
import datetime
from matplotlib.ticker import FuncFormatter as _FuncFormatter
from matplotlib.dates import DayLocator as _DayLocator
import os as _os

unit_time = 'days since 1900-01-01'



def load_csv(fname):
    """Loads the dat of a saved timesereis instance and creates a new TimeSeries instance

    Arguments
    ---------
    fname: str.
        Path to the file to load"""
    data = _pd.read_csv(fname, index_col=0)
    data.index = _pd.to_datetime(data.index)
    return TimeSeries(data)

def none2nan(var):
    if type(var).__name__ == 'NoneType':
        var = _np.nan
    return var

# def nan2none(var):
#     if type(var).__name__ == 'str':
#         pass
#     elif _np.isnan(var):
#         var = None
#     return var

def save_netCDF(ts,fname, leave_open = False):

    # if ts._time_format == 'timedelta':
    #     ts.timed

    file_mode = 'w'
    try:
        ni = _Dataset(fname, file_mode)
    except RuntimeError:
        if _os.path.isfile(fname):
            _os.remove(fname)
            ni = _Dataset(fname, file_mode)

    time_dim = ni.createDimension('time', ts.data.shape[0])
    dim_data_col = ni.createDimension('data_columns', ts.data.shape[1])

    ts_time_num = _date2num(ts.data.index.to_pydatetime(), unit_time)#.astype(float)
    time_var = ni.createVariable('time', ts_time_num.dtype, 'time')
    time_var[:] = ts_time_num
    time_var.units = 'days since 1900-01-01'

    var_data = ni.createVariable('data', ts.data.values.dtype, ('time', 'data_columns'))
    var_data[:] = ts.data.values

    ts_columns = ts.data.columns.values.astype(str)
    var_data_collumns = ni.createVariable('data_columns', ts_columns.dtype, 'data_columns')
    var_data_collumns[:] = ts_columns

    ni._ts_type = type(ts).__name__
    ni._data_period = none2nan(ts._data_period)
    ni._x_label = none2nan(ts._x_label)
    ni._y_label =  none2nan(ts._y_label)
    ni.info = none2nan(ts.info)
    ni._atm_py_commit = _git_tools.current_commit()

    if leave_open:
        return ni
    else:
        ni.close()

def load_netCDF(fname):

    ni = _Dataset(fname, 'r')

    # load time
    time_var = ni.variables['time']
    time_var.units
    ts_time = _num2date(time_var[:], time_var.units)
    timestamp = _pd.DatetimeIndex(ts_time)

    # load  data
    var_data = ni.variables['data']
    ts_data = _pd.DataFrame(var_data[:], index=timestamp)

    # load column names
    var_data_col = ni.variables['data_columns']
    ts_data = _pd.DataFrame(var_data[:], index=timestamp,
                           columns=var_data_col[:])

    # test which type of timeseries (1D, 2D, 3D)
    try:
        ts_type = ni.getncattr('_ts_type')
    except AttributeError:
        ts_type = 'TimeSeries'
    # create time series
    if ts_type == 'TimeSeries':
        ts_out = TimeSeries(ts_data)
    elif ts_type == 'TimeSeries_2D':
        ts_out = TimeSeries_2D(ts_data)
    elif ts_type == 'TimeSeries_3D':
        ts_out = TimeSeries_3D(ts_data)

    # load attributes and attach to time series
    for atr in ni.ncattrs():
        value = ni.getncattr(atr)
        # there is a bug in pandas where it does not like numpy types ->
        if type(value).__name__ == 'str':
            pass
        elif 'float' in value.dtype.name:
            value = float(value)
        elif 'int' in value.dtype.name:
            value = int(value)
        # netcdf did not like NoneType so i converted it to np.nan. Here i am converting back.
        elif _np.isnan(value):
            value = None

        setattr(ts_out, atr, value)

    ni.close()
    return ts_out


#### Tools
def close_gaps(ts, verbose = False):
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


def align_to(ts, ts_other, verbose= False):
    """
    Main change, timestamp at beginning!
    Align the TimeSeries ts to another time_series by interpolating (linearly). If
    data periods differe by at least a factor of 2 a rolling mean is calculated
    with a window size equal to the ratio (if ratio is positive !!!).

    Notes
    -----
    For the obvious reason it is recommended to align the time series with the smaller data period to that
    with the larger period.

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
    if verbose:
        print('=================================')
        print('=====  perform alignment ========')

    if _np.array_equal(ts.data.index, ts_other.data.index):
        if verbose:
            print('indeces are identical, returning original time series.')
        return ts

    window = ts_other._data_period / ts._data_period
    if window < 0.5:
        _warnings.warn('Time period of other time series is smaller (ratio: %s). You might want '
                      'align the other time series with this one instead?'%window)
    window = int(round(window))

    if window > 2:
        if verbose:
            print('Data period difference larger than a factor of 2 -> performing rolling mean')

        roll = ts.data.rolling(window,
                               min_periods=1,
                               # center=True
                               )
        dfrm = roll.mean()
        dfrm = dfrm.shift( - window) # rolling puts the timestamp to the end of the window, this shifts it to the beginning of the window
        tsrm = TimeSeries(_pd.DataFrame(dfrm))
        # tsrm._data_period = ts._data_period
    else:
        if verbose:
            print('Data period difference smaller than a factor of 2 -> do nothing')
        tsrm = ts

    ts_other.data = ts_other.data.loc[:,[]]
    ts_other.data.columns.name = None # if this is not empty it will give an error
    if verbose:
        print('performing merge with empty index of other time series')
    ts_t =  merge(ts_other, tsrm, verbose = verbose)
    tsrm.data = ts_t.data

    tsrm._data_period = ts_other._data_period
    if verbose:
        print('=====  alignment done ========')
        print('=================================')

    return tsrm

def align_to_old(ts, ts_other, verbose= False):
    """
    Align the TimeSeries ts to another time_series by interpolating (linearly). If
    data periods differe by at least a factor of 2 a rolling mean is calculated
    with a window size equal to the ratio (if ratio is positive !!!).

    Notes
    -----
    For the obvious reason it is recommended to align the time series with the smaller data period to that
    with the larger period.

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
    if verbose:
        print('=================================')
        print('=====  perform alignment ========')
    # if _np.all(ts.data.index == ts_other.data.index):
    if _np.array_equal(ts.data.index, ts_other.data.index):
        if verbose:
            print('indeces are identical, returning original time series.')
        return ts

    window = ts_other._data_period / ts._data_period
    if window < 0.5:
        _warnings.warn('Time period of other time series is smaller (ratio: %s). You might want '
                      'align the other time series with this one instead?'%window)
    window = int(round(window))

    if window > 2:
        if verbose:
            print('Data period difference larger than a factor of 2 -> performing rolling mean')

        roll = ts.data.rolling(window,
                               min_periods=1,
                               center=True)
        dfrm = roll.mean()

        tsrm = TimeSeries(_pd.DataFrame(dfrm))
        # tsrm._data_period = ts._data_period
    else:
        if verbose:
            print('Data period difference smaller than a factor of 2 -> do nothing')
        tsrm = ts

    ts_other.data = ts_other.data.loc[:,[]]
    ts_other.data.columns.name = None # if this is not empty it will give an error
    if verbose:
        print('performing merge with empty index of other time series')
    ts_t =  merge(ts_other, tsrm, verbose = verbose)
    tsrm.data = ts_t.data
    tsrm._data_period = ts_other._data_period
    if verbose:
        print('=====  alignment done ========')
        print('=================================')

    return tsrm


def merge(ts, ts_other, verbose = False):
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
    if verbose:
        print('=================================')
        print('=====  perform merge ========')

    ts_this = ts.copy()

    # if _np.all(ts_this.data.index == ts_other.data.index):
    if _np.array_equal(ts_this.data.index, ts_other.data.index):
        ts_this.data = _pd.concat([ts_this.data, ts_other.data], axis=1)

    else:
        ts_data_list = [ts_this.data, ts_other.data]
        # catsortinterp = _pd.concat(ts_data_list).sort_index().interpolate(method='index')
        # merged = catsortinterp.groupby(catsortinterp.index).mean().reindex(ts_data_list[0].index)
        # ts_this.data = merged

        # ts_data_list = [data.data, correlant_nocol]
        catsort = _pd.concat(ts_data_list).sort_index()  # .interpolate(method='index')

        mask = catsort.copy()
        grp = ((mask.notnull() != mask.shift().notnull()).cumsum())
        grp['ones'] = 1
        for i in catsort.columns:
            mask[i] = (grp.groupby(i)['ones'].transform('count') < 3) | catsort[i].notnull()

        catsortinterp = catsort.interpolate(method='index')
        catsortinterpmasked = catsortinterp[mask]

        merged = catsortinterpmasked.groupby(catsortinterpmasked.index).mean().reindex(ts_data_list[0].index)
        ts_this.data = merged

    if verbose:
        print('=====  merge done ========')
        print('==========================')
    return ts_this

def concat(ts_list):
    for ts in ts_list:
        if type(ts).__name__ != 'TimeSeries':
            raise TypeError('Currently works only with TimeSeries not with %s'%(type(ts).__name__))
    ts = ts_list[0].copy()
    ts.data = _pd.concat([i.data for i in ts_list])

    return ts


def correlate(data,correlant, data_column = False, correlant_column = False, remove_zeros=True):
    """Correlates data in correlant to that in data. In the process the data in correlant
    will be aligned to that in data. Make sure that data has the lower period (less data per period of time)."""
    if data_column:
        data_values = data.data[data_column].values
    elif data.data.shape[1] > 1:
        raise ValueError('Data contains more than 1 column. Specify which to correlate. Options: %s'%(list(data.data.keys())))
    else:
        data_values = data.data.iloc[:,0].values
    correlant_aligned = correlant.align_to(data)
    if correlant_column:
        correlant_values = correlant_aligned.data[correlant_column].values
    elif correlant.data.shape[1] > 1:
        raise ValueError('''Correlant contains more than 1 column. Specify which to correlate. Options:
%s'''%(list(correlant_aligned.data.keys())))
    else:
        correlant_values = correlant_aligned.data.iloc[:,0].values

    out = _array_tools.Correlation(data_values, correlant_values, remove_zeros=remove_zeros, index = data.data.index)
    out._x_label_orig = 'DataTime'
    return out

def rolling_correlation(data, correlant, window, data_column = False, correlant_column = False,  min_good_ratio = 0.67, verbose = True):
    "time as here: http://docs.scipy.org/doc/numpy/reference/arrays.datetime.html#datetime-units"

    if correlant_column:
        correlant = correlant._del_all_columns_but(correlant_column)

    if data_column:
        data = data._del_all_columns_but(data_column)

    correlant = correlant.align_to(data) # I do align before merge, because it is more suffisticated!
    merged = data.copy()
    merged.data['correlant'] = correlant.data

    data_period = _np.timedelta64(int(merged._data_period), 's')
    window = _np.timedelta64(window[0], window[1])
    window = int(window/data_period)
    if verbose:
        print('Each window contains %s data points of which at least %s are not nan.'%(window, int(window * min_good_ratio)))

    min_good = window * min_good_ratio
    size = merged.data.shape[0]-window + 1
    timestamps = _pd.TimeSeries(_pd.to_datetime(_pd.Series(_np.zeros(size))))
    pear_r = _np.zeros(size)
    for i in range(size):
        secment = TimeSeries(merged.data.iloc[i:i+window,:])
        secment._data_period = merged._data_period
    #     print(secment.data.dropna().shape[0] < min_good)
        if secment.data.dropna().shape[0] < min_good:
            pear_r[i]= _np.nan
        else:
            corr = secment.correlate_to(secment, data_column=merged.data.columns[0], correlant_column=merged.data.columns[1])
            pear_r[i] = corr.pearson_r[0]
        timestamps.iloc[i] = secment.data.index[0] + ((secment.data.index[-1] - secment.data.index[0])/2.)
    #     break

    pear_r_ts = TimeSeries(_pd.DataFrame(pear_r, index = timestamps, columns = ['pearson_r']))
    pear_r_ts._data_period = merged._data_period
    pear_r_ts._y_label = 'r'
    return pear_r_ts


class Rolling(_pd.core.window.Rolling):
    def __init__(self, obj, window, min_good_ratio=0.67,
                 verbose=True,center = True,
                 **kwargs):
        self.data = obj
        window = _np.timedelta64(window[0], window[1])
        window = int(window / _np.timedelta64(int(obj._data_period), 's'))
        min_periods = int(window * min_good_ratio)
        if obj.data.columns.shape[0] == 1:
            self._data_column = obj.data.columns[0]
        else:
            txt = 'please make sure the timeseries has only one collumn'
            raise ValueError(txt)

        if verbose:
            print('Each window contains %s data points of which at least %s are not nan.' % (window,
                                                                                             min_periods))
        super().__init__(obj.data[self._data_column],
                         window,
                         min_periods=min_periods,
                         center = center,
                         **kwargs)

    #         print(self.window)

    def corr(self, other, *args, **kwargs):
        if other.data.columns.shape[0] == 1:
            other_column = other.data.columns[0]
        else:
            txt = 'please make sure the timeseries has only one collumn'
            raise ValueError(txt)
        # other_column = 'Bs_G_Dry_1um_Neph3W_1'
        other = other.align_to(self.data)
        other = other.data[other_column]
        corr_res = super().corr(other, *args, **kwargs)
        corr_res_ts = TimeSeries(_pd.DataFrame(corr_res))
        corr_res_ts._data_period = self.data._data_period
        return corr_res_ts

    def corr_timelag(self, other, dt=(5, 'm'), no_of_steps=10, center=0, normalize=True, **kwargs):
        """dt: tuple
                first arg of tuple can be int or array-like of dtype int. Second arg is unit. if array-like no_of... is ignored"""
        if other.data.columns.shape[0] == 1:
            other_column = other.data.columns[0]
        else:
            txt = 'please make sure the timeseries has only one collumn'
            raise ValueError(txt)
        # other =  acsm.copy()
        #         dt = (5, 'm')
        #         no_of_seps = 10
        #         center = 0


        if hasattr(dt[0], '__len__'):
            if type(dt[0]).__name__ == 'list':
                dt_array = _np.array(dt[0])
            else:
                dt_array = dt[0]
            no_of_steps = len(dt_array)
            center = 0
        else:
            dt_array = _np.arange(0, dt[0] * no_of_steps, dt[0]) - int(no_of_steps * dt[0] / 2)

        if center:
            dt_array += int(center)
        out = False
        for dtt in dt_array:
            tst = other.copy()
            tst.data.index += _np.timedelta64(int(dtt), dt[1])
            corr = self.corr(tst)
            if not out:
                out = corr
                out.data.columns = [dtt]
            else:
                out.data[dtt] = corr.data  # [self._data_column]

        if normalize:
            out = TimeSeries_2D(out.data.apply(lambda line: line / line.max(), axis=1))
        else:
            out = TimeSeries_2D(out.data)

        def aaa(line):
            if (_np.isnan(line)).sum() > ((~_np.isnan(line)).sum() * 0.3):
                return _np.nan
            col_t = cols.copy()
            lt = line[~ _np.isnan(line)]

            col_t = col_t[~ _np.isnan(line)]
            argmax = lt.argmax()
            realMax = col_t[argmax]
            return realMax

        cols = out.data.columns
        dt_max = _np.apply_along_axis(aaa, 1, out.data.values)



        # cols = out.data.columns
        # dt_max = _np.apply_along_axis(lambda arr: cols[arr.argmax()], 1, out.data.values)
        dt_max = TimeSeries(_pd.DataFrame(dt_max, index=out.data.index))
        if dt[1] == 'm':
            ylt = 'min.'
        else:
            ylt = dt[1]
        dt_max._y_label = 'Time lag (%s)' % ylt
        return out, dt_max


class Rolling_old(object):
    # def __init__(self, data):
    #     self.data = data

    def __init__(self, data, correlant, window, data_column=False,
                 correlant_column=False, min_good_ratio=0.67, verbose=True):
        "time as here: http://docs.scipy.org/doc/numpy/reference/arrays.datetime.html#datetime-units"

        if correlant_column:
            correlant = correlant._del_all_columns_but(correlant_column)

        if data_column:
            data = data._del_all_columns_but(data_column)


        correlant = correlant.align_to(data)  # I do align before merge, because it is more suffisticated!
        merged = data.copy()
        merged.data.columns = ['data']
        merged.data['correlant'] = correlant.data

        #         self.data = data
        #         self.correlant = correlant
        self._merged = merged
        self._data_column = data_column
        self._correlant_column = correlant_column
        self._min_good_ratio = min_good_ratio
        self._verbose = verbose
        self._data_period = _np.timedelta64(int(merged._data_period), 's')
        window = _np.timedelta64(window[0], window[1])
        self._window = int(window / self._data_period)
        if verbose:
            print('Each window contains %s data points of which at least %s are not nan.' % (self._window,
                                                                                             int(self._window * min_good_ratio)))
        self._min_good = self._window * min_good_ratio
        self._size = merged.data.shape[0] - self._window + 1
        self._timestamps = _pd.TimeSeries(_pd.to_datetime(_pd.Series(_np.zeros(self._size))))
        # self.dt_min = dt_min
        # self.center = center
        # self.no_of_steps = no_of_steps

    def correlation(self):
        # self._prepare(correlant, window, data_column=data_column, correlant_column=correlant_column, min_good_ratio=min_good_ratio, verbose=verbose)
        out = self._loop_over_segments_and_apply_funk('correlation')
        out._y_label = 'r'
        return out

    def timelaps_correlation(self,
                             # correlant, window,
                             dt_min=10, center=0,
                             no_of_steps=10,
                             # data_column=False, correlant_column=False,
                             # min_good_ratio=0.67, verbose=True
                             # output_dimentions = 1
                             ):
        # """"""
        # self._prepare(correlant, window, data_column=data_column, dt_min=dt_min, center=center,
        #               no_of_steps=no_of_steps,
        #               correlant_column=correlant_column, min_good_ratio=min_good_ratio, verbose=verbose)
        self.dt_min = dt_min
        self.center = center
        self.no_of_steps = no_of_steps
        out = self._loop_over_segments_and_apply_funk('timelaps_correlation')

        return out

    def _timelap_correlation(self, secment):
        corcoeffs = _np.zeros((self.no_of_steps))
        dataVcorr_mean = _np.zeros((self.no_of_steps))
        dataVcorr_std = _np.zeros((self.no_of_steps))
        for e, i in enumerate(self.timelaps_shifts):
            datatmp = secment._del_all_columns_but('data')
            datatmp.data.index += _pd.Timedelta('%s minutes' % i)
            out = datatmp.correlate_to(secment._del_all_columns_but('correlant'))
            corcoeffs[e] = out.pearson_r[0]
            dataVcorr_mean[e] = (out._data / out._correlant).mean()
            dataVcorr_std[e] = (out._data / out._correlant).std()
        time_at_max = self.timelaps_shifts[corcoeffs.argmax()]
        # corframe = _pd.DataFrame(corcoeffs, index = dt_array)
        # corframe = corcoeffs
        out = {'time_at_max': time_at_max,
               'time_lapse_corr': corcoeffs / corcoeffs.max()}
        return out

    def _loop_over_segments_and_apply_funk(self, what):
        out = _np.zeros(self._size)
        out[:] = _np.nan
        if what == 'correlation':
            colname = 'pearson_r'
        elif what == 'timelaps_correlation':
            colname = 'dt at max correlations'
            time_lapse_corr = _np.zeros((self._size, self.no_of_steps))
            time_lapse_corr[:] = _np.nan
            start_dt = (- int(self.no_of_steps / 2) + self.center) * self.dt_min
            end_dt = (int(self.no_of_steps / 2) + self.center) * self.dt_min
            dt_array = _np.linspace(start_dt, end_dt, self.no_of_steps)
            self.timelaps_shifts = dt_array
        else:
            corname = 'blablabla'
        for i in range(self._size):
            secment = TimeSeries(self._merged.data.iloc[i:i + self._window, :])
            secment._data_period = self._merged._data_period
            #     print(secment.data.dropna().shape[0] < min_good)
            if secment.data.dropna().shape[0] < self._min_good:
                pass
            else:
                if what == 'correlation':
                    corr = secment.correlate_to(secment, data_column=self._merged.data.columns[0],
                                                correlant_column=self._merged.data.columns[1])
                    out[i] = corr.pearson_r[0]
                if what == 'timelaps_correlation':
                    tlc = self._timelap_correlation(secment)
                    out[i] = tlc['time_at_max']
                    time_lapse_corr[i] = tlc['time_lapse_corr']

            self._timestamps.iloc[i] = secment.data.index[0] + ((secment.data.index[-1] - secment.data.index[0]) / 2.)
        # break

        out_ts = TimeSeries(_pd.DataFrame(out, index=self._timestamps, columns=[colname]))
        out_ts._y_label = colname
        out_ts._data_period = self._merged._data_period
        if what == 'timelaps_correlation':
            time_lapse_corr_ts = TimeSeries_2D(_pd.DataFrame(time_lapse_corr, index = self._timestamps, columns=self.timelaps_shifts))
            time_lapse_corr_ts._data_period = self._merged._data_period
            time_lapse_corr_ts._y_label = 'Time shift (min.)'

            # out = {'time_at_max': time_at_max_ts,
            #        'time_lapse_corr': time_lapse_corr_ts}
            return out_ts, time_lapse_corr_ts
        else:
            return out_ts


class WrappedPlot(list):
    def __init__(self, axes_list):
        for a in axes_list:
            self.append(a)


def plot_wrapped(ts,periods = 1, frequency = 'h', ylabel = 'auto', max_wraps = 10, ylim = None, ax = None, twin_x = None, **plot_kwargs):
    """frequency: http://docs.scipy.org/doc/numpy/reference/arrays.datetime.html#datetime-units

    if ax is set, all other parameters will be ignored
    ylim: set to False if you don't want """
    if 'cb_kwargs' in plot_kwargs.keys():
        cb_kwargs = plot_kwargs.pop('cb_kwargs')
    else:
        cb_kwargs = False

    if ax:
        periods = ax._periods
        frequency =  ax._frequency
        ylabel = ax._ylabel
        max_wraps = ax._max_wraps
        ylim_old = ax._ylim
        col_no = _np.array([len(a.get_lines()) for a in ax]).max()

    if periods >1:
        raise ValueError('Sorry periods larger one is not working ... consider fixing it?!?')
    start, end = ts.get_timespan()
    length = end - start

    if frequency == 'Y':
        start = _np.datetime64(str(start.year))
        end = _np.datetime64(str(end.year))
        periods_no = int((end - start + 1) / _np.timedelta64(1, 'Y'))
        autofmt_xdate = False
        xlabel = 'Month of year'

    elif frequency == 'M':
        start = _np.datetime64(str(start.year)) + _np.timedelta64(start.month - 1, 'M')
        end = _np.datetime64(str(end.year)) + _np.timedelta64(end.month - 1, 'M')
        periods_no = int((end - start + 1) / _np.timedelta64(1, 'M'))
        autofmt_xdate = False
        xlabel = 'Day of month'


    elif frequency == 'D':
        start = _np.datetime64(str(start.year)) + _np.timedelta64(start.month - 1, 'M') + _np.timedelta64(start.day - 1, 'D')
        end = _np.datetime64(str(end.year)) + _np.timedelta64(end.month - 1, 'M') + _np.timedelta64(end.day - 1, 'D')
        periods_no = int((end - start + 1) / _np.timedelta64(1, 'D'))
        autofmt_xdate = True
        xlabel = 'Hour of day'

    else:
        periods_no = int(_np.ceil(length/_np.timedelta64(periods, frequency))) - 1
        autofmt_xdate = False
        xlabel = 'Time'

    start_t = start
    if periods_no > max_wraps:
        raise ValueError("To many wraps (%i). Change frequency or max_wraps."%periods_no)
    if ax:
        a = ax
        f = a[0].get_figure()
    else:
        height_ratios = [1] * periods_no
        # if cb_kwargs:
        #     no_of_ax = periods_no + 1
        #     height_ratios = [0.08] + height_ratios
        # else:
        no_of_ax = periods_no
        f,a = _plt.subplots(no_of_ax, sharex=True, gridspec_kw={'hspace': 0,
                                                                'height_ratios': height_ratios})
        # twins_x = []
        # for at in a:
        #     twins_x.append(at.twinx())
        # if cb_kwargs:
        #     a_cb = a[0]
        #     a = a[1:]

        f.set_figheight(3*periods_no)
        col_no = 0
    bbox_props = dict(boxstyle="round,pad=0.3", fc=[1,1,1,0.8], ec="black", lw=1)

    if twin_x:
        if not hasattr(a,'twins_x'):
            twins_x = []
            for at in a:
                twins_x.append(at.twinx())
        else:
            twins_x = a.twins_x

    if not ylim:
        ylim = [ts.data.min().min(), ts.data.max().max()]
        if ax:
            if ylim_old[0] < ylim[0]:
                ylim[0] = ylim_old[0]
            if ylim_old[1] > ylim[1]:
                ylim[1] = ylim_old[1]

    for i in range(int(periods_no)):
        end_t = start_t + _np.timedelta64(periods, frequency)
        try:
            tst = ts.zoom_time(start_t, end_t)
        except IndexError:
            # print('voll der index error')
            tst = False

        if twin_x:
            at = twins_x[i]
        else:
            at = a[i]
        if 1:
            txtpos = (0.05,0.8)

            text = str(start_t).split(' ')
            if frequency == 'h':
                text = text[1]
            elif frequency == 'M':
                def day_formatter(x, pos):
                    dt = _np.datetime64('0') + _np.timedelta64(int(x), 'D')
                    dt = _pd.Timestamp(dt)
                    out = dt.day
                    return out
                text = start_t

                df = _FuncFormatter(day_formatter)
                at.xaxis.set_major_formatter(df)
                at.xaxis.set_major_locator(_DayLocator(interval=4))

            elif frequency == 'Y':
                text = text[0].split('-')[0]
            elif frequency == 'D':
                text = start_t
            else:
                text = 'not set'
            at.text(txtpos[0], txtpos[1], text, transform=at.transAxes, bbox=bbox_props)


        if tst:
            tst.data.index  = tst.data.index - start_t
            tst.data.index += _np.datetime64('1900')
            if type(tst).__name__ == 'TimeSeries':
                if twin_x:
                    # plt_out = tst.plot(ax=at, color=_plt_tools.color_cycle[col_no])
                    plt_out = tst.plot(ax=at, autofmt_xdate=autofmt_xdate, color=_plt_tools.color_cycle[col_no], **plot_kwargs)
                else:
                    plt_out = tst.plot(ax=at, autofmt_xdate = autofmt_xdate, color = _plt_tools.color_cycle[col_no], **plot_kwargs)
            if type(tst).__name__ == 'TimeSeries_2D':
                # if 'cb_kwargs' in plot_kwargs.keys():
                #     cb_kwargs = plot_kwargs['cb_kwargs']
                # else:
                #     cb_kwargs = False
                # if i == len(a) - 1:
                #     cb_kwargs_t = cb_kwargs
                # else:
                #     cb_kwargs_t = False
                plt_out = tst.plot(ax=at, autofmt_xdate=autofmt_xdate, color=_plt_tools.color_cycle[col_no], cb_kwargs = False, **plot_kwargs)
                plt_out[2].set_clim(ylim)

        if type(tst).__name__ == 'TimeSeries':
            if not twin_x:
                at.set_ylim(ylim)
        # formatter = FuncFormatter(timeTicks)
        # at.xaxis.set_major_formatter(formatter)



        # from matplotlib.dates import DayLocator
        # if frequency == 'M':
        #     at.xaxis.set_major_locator(DayLocator())


        start_t = end_t
        at.set_ylabel('')

        # if an existing set of axes was provided the following will prevent an IndexError
        if i == len(a) - 1:
            # print('did the break')
            break

    if cb_kwargs:
        if 'shrink' not in cb_kwargs.keys():
            cb_kwargs['shrink'] = 0.5
        if 'anchor' not in cb_kwargs.keys():
            cb_kwargs['anchor'] = (0, 1)
        if 'pad' not in cb_kwargs.keys():
            cb_kwargs['pad'] = 0.01
            
        cb = f.colorbar(plt_out[2], ax=a.ravel().tolist(),
                        # shrink=cb_kwargs['shrink'], anchor=cb_kwargs['anchor'],
                        **cb_kwargs)

    # at.set_xlim(left=0, right = _np.timedelta64(periods,frequency)/_np.timedelta64(1,'ns'))
    # at.set_xlim(left=0, right=_np.timedelta64(periods, frequency) / _np.timedelta64(1, 'ns'))
    if ylabel == 'auto':
        ylabel = ts._y_label
    _plt_tools.set_shared_ylabel(a,ylabel)
    at.set_xlabel(xlabel)
    out = WrappedPlot(a)
    out._periods = periods
    out._frequency = frequency
    out._ylabel = ylabel
    out._max_wraps = max_wraps
    out._ylim = ylim
    if cb_kwargs:
        out.cb = cb
    if twin_x:
        if not hasattr(a,'twins_x'):
            out.twins_x = twins_x

    return out


class TimeSeries(object):
    """
    TODO: depersonalize!!! Try composition approach.

    This class simplifies the handling of housekeeping information from measurements.
    Typically this class is created by a housekeeping function of the particular instruments.

    Notes
    -----
    Make sure the passed pandas DataFrame includes the following column names:
     - barometric_pressure
     - altitude

    Attributes
    ----------
    data:  pandas DataFrame with index=DateTime and columns = housekeeping parameters
    """

    def __init__(self, data, info=None):
        # if not type(data).__name__ == 'DataFrame':
        #     raise TypeError('Data has to be of type DataFrame. It currently is of type: %s'%(type(data).__name__))
        self._data_period = None
        self.data = data
        self.info = info
        self._y_label = ''
        self._x_label = 'Time'
        self._time_format = 'datetime' #'timedelta'
        if hasattr(self.data, 'index'):
            self._start_time = self.data.index[0]
        elif hasattr(self.data, 'major_axis'):
            self._start_time = self.data.major_axis[0]
        else:
            txt = 'Looks like data is neither series nor panel ... what is it!?!'
            raise AttributeError(txt)

    def __str__(self):
        return self.data.__str__()

    def __repr__(self):
        return self.data.__repr__()

    def __truediv__(self,other):
        self_t = self.copy()

        if type(other).__name__ in ['int', 'float']:
            self_t.data /= other
            return self_t

        else:
            other = other.copy()
            if self_t._data_period > other._data_period:
                other = other.align_to(self_t)
                # other._data_period = self_t._data_period
            else:
                self_t = self_t.align_to(other)
                # self_t._data_period = other._data_period
            # return self_t,other
            if other.data.shape[1] == 1:
                out = self_t.data.divide(other.data.iloc[:,0], axis = 0)
            elif self_t.data.shape[1] == 1:
                out = other.data.divide(self_t.data.iloc[:,0], axis = 0)
                out = 1/out
            else:
                txt = 'at least one of the dataframes have to have one column only'
                raise ValueError(txt)

            ts = TimeSeries(out)
            ts._data_period = self_t._data_period
            return ts

    def __add__(self,other):
        self = self.copy()
        other = other.copy()
        if self._data_period > other._data_period:
            other = other.align_to(self)
            other._data_period = self._data_period
        else:
            self = self.align_to(other)
            self._data_period = other._data_period

        if other.data.shape[1] == 1:
            out = self.data.add(other.data.iloc[:,0], axis = 0)
        elif self.data.shape[1] == 1:
            out = other.data.add(self.data.iloc[:,0], axis = 0)
        else:
            txt = 'at least one of the dataframes have to have one column only'
            raise ValueError(txt)

        ts = TimeSeries(out)
        ts._data_period = self._data_period
        return ts

    def __sub__(self,other):
        self = self.copy()
        other = other.copy()
        if self._data_period > other._data_period:
            other = other.align_to(self)
            other._data_period = self._data_period
        else:
            self = self.align_to(other)
            self._data_period = other._data_period

        if other.data.shape[1] == 1:
            out = self.data.sub(other.data.iloc[:,0], axis = 0)
        elif self.data.shape[1] == 1:
            out = other.data.sub(self.data.iloc[:,0], axis = 0)
            out = - out
        else:
            txt = 'at least one of the dataframes have to have one column only'
            raise ValueError(txt)

        ts = TimeSeries(out)
        ts._data_period = self._data_period
        return ts

    def __mul__(self,other):
        self = self.copy()
        other = other.copy()
        if self._data_period > other._data_period:
            other = other.align_to(self)
            other._data_period = self._data_period
        else:
            self = self.align_to(other)
            self._data_period = other._data_period

        if other.data.shape[1] == 1:
            out = self.data.multiply(other.data.iloc[:,0], axis = 0)
        elif self.data.shape[1] == 1:
            out = other.data.multiply(self.data.iloc[:,0], axis = 0)
        else:
            txt = 'at least one of the dataframes have to have one column only'
            raise ValueError(txt)

        ts = TimeSeries(out)
        ts._data_period = self._data_period
        return ts

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, data):
        if not type(data).__name__ == 'DataFrame':
            raise TypeError('Data has to be of type DataFrame. It currently is of type: %s'%(type(data).__name__))
        self.__data = data



    def convert2verticalprofile(self, alt_label = None, alt_timeseries = None):
        ts_tmp = self.copy()
    #     hk_tmp.data['Time'] = hk_tmp.data.index
    #     if alt_label:
    #         label = alt_label
    #     else:
    #         label = 'Altitude'
        if alt_timeseries:
            alt_timeseries = alt_timeseries.align_to(ts_tmp)
            _pandas_tools.ensure_column_exists(alt_timeseries.data, 'Altitude', col_alt=alt_label)
            ts_tmp.data.index = alt_timeseries.data['Altitude']
        else:
            _pandas_tools.ensure_column_exists(ts_tmp.data, 'Altitude', col_alt=alt_label)
            ts_tmp.data.index = ts_tmp.data['Altitude']
        out = atmPy.general.vertical_profile.VerticalProfile(ts_tmp.data)
        out._x_label = self._y_label
        return out

    def _del_all_columns_but(self, keep, inplace = False):
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

    def datetime2timedelta(self):
        """Sets the time index so that it starts at zero"""
        if self._time_format == 'timedelta':
            return self

        ts = self.copy()
        data = ts.data
        time_from_start = data.index - ts._start_time
        ts.data.index = time_from_start
        ts._time_format = 'timedelta'
        return ts

    def timedelta2datetime(self):
        """Sets the time index so that it starts at zero"""
        if self._time_format == 'datetime':
            return self

        ts = self.copy()
        data = ts.data
        time = data.index + ts._start_time
        ts.data.index = time
        ts._time_format = 'datetime'
        return ts

    def average_time(self, window, std = False, envelope = False):
        """Massive change: time stamp at beginning! returns a copy of the sizedistribution_TS with reduced size by averaging over a given window.
        The difference to panda's resample is that it takes a time window instead of a point window.

        Arguments
        ---------
        window: tuple
            tuple[0]: periods
            tuple[1]: frequency (Y,M,W,D,h,m,s...) according to:
                http://docs.scipy.org/doc/numpy/reference/arrays.datetime.html#datetime-units
                if error also check out:
                http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases

        Returns
        -------
        SizeDistribution_TS instance
            copy of current instance with resampled data frame
        """

        ts = self.copy()
        # ts.data = ts.data.resample(window, closed='right', label='right').mean() #old

        # determine offset, so the label is in the center
        # toff = _np.timedelta64(window[0], window[1]) / _np.timedelta64(2, 's')
        # toff = '%iS' % toff

        ts._data_period = _np.timedelta64(window[0], window[1]) / _np.timedelta64(1, 's')

        resample = ts.data.resample(window,
                            label = 'left',
                            # loffset = toff,
                           )

        ts.data = resample.mean()
        if std or envelope:
            std_tmp = resample.std()
            if std:
                ts.data['std'] = std_tmp
            if envelope:
                ts.data['envelope_low'] = ts.data.iloc[:,0] - std_tmp.iloc[:,0]
                ts.data['envelope_high'] = ts.data.iloc[:,0] + std_tmp.iloc[:,0]

        ts._start_time = ts.data.index[0]
        return ts


    def average_time_old(self, window, std = False, envelope = False):
        """returns a copy of the sizedistribution_TS with reduced size by averaging over a given window.
        The difference to panda's resample is that it takes a time window instead of a point window.

        Arguments
        ---------
        window: tuple
            tuple[0]: periods
            tuple[1]: frequency (Y,M,W,D,h,m,s...) according to:
                http://docs.scipy.org/doc/numpy/reference/arrays.datetime.html#datetime-units
                if error also check out:
                http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases

        Returns
        -------
        SizeDistribution_TS instance
            copy of current instance with resampled data frame
        """

        ts = self.copy()
        # ts.data = ts.data.resample(window, closed='right', label='right').mean() #old

        # determine offset, so the label is in the center
        toff = _np.timedelta64(window[0], window[1]) / _np.timedelta64(2, 's')
        ts._data_period = _np.timedelta64(window[0], window[1]) / _np.timedelta64(1, 's')
        toff = '%iS' % toff

        resample = ts.data.resample(window,
                            label = 'left',
                            loffset = toff,
                           )

        ts.data = resample.mean()
        if std or envelope:
            std_tmp = resample.std()
            if std:
                ts.data['std'] = std_tmp
            if envelope:
                ts.data['envelope_low'] = ts.data.iloc[:,0] - std_tmp.iloc[:,0]
                ts.data['envelope_high'] = ts.data.iloc[:,0] + std_tmp.iloc[:,0]

        ts._start_time = ts.data.index[0]
        return ts

    align_to = align_to
    align_to_new = align_to_old
    # def align_to(self, ts_other):
    #     return align_to(self, ts_other)

    close_gaps = close_gaps

    correlate_to = correlate

    rolling_correlation = rolling_correlation

    merge = merge

    copy = _deepcopy

    # rollingR = Rolling

    def rolling(self, window, data_column=False,
                 correlant_column=False, min_good_ratio=0.67, verbose=True):
        return Rolling(self, window, data_column=data_column,
                 correlant_column=correlant_column, min_good_ratio=min_good_ratio, verbose=verbose)

    def rolling_old(self, correlant, window, data_column=False,
                correlant_column=False, min_good_ratio=0.67, verbose=True):
        return Rolling_old(self, correlant, window, data_column=data_column,
               correlant_column=correlant_column, min_good_ratio=min_good_ratio, verbose=verbose)

    def plot(self, ax = None, legend = True, label = None, autofmt_xdate = True, **kwargs):
        """Plot each parameter separately versus time
        Arguments
        ---------
        same as pandas.plot

        Returns
        -------
        list of matplotlib axes object """

        def timeTicks(x, pos):
            d = datetime.timedelta(seconds=x * 1e-9)
            return str(d)

        # a = self.data.plot(**kwargs)
        if not ax:
            f,ax = _plt.subplots()
        else:
            f = ax.get_figure()

        did_plot = False  # had to implement that since matploglib cept crashing when all where nan
        for k in self.data.keys():
            if not label:
                label_t = k
            else:
                label_t = label

            if _np.all(_np.isnan(self.data[k].values)):
                continue

            ax.plot(self.data.index, self.data[k].values, label = label_t, **kwargs)

            if self._time_format == 'timedelta':
                formatter = _FuncFormatter(timeTicks)
                ax.xaxis.set_major_formatter(formatter)

            did_plot = True

        ax.set_xlabel(self._x_label)
        ax.set_ylabel(self._y_label)
        if len(self.data.keys()) > 1:
            ax.legend()

        if did_plot:
            if autofmt_xdate:
                f.autofmt_xdate()

        return ax

    plot_wrapped = plot_wrapped

    def zoom_time(self, start=None, end=None, copy=True):
        """ Selects a strech of time from a housekeeping instance.

        Arguments
        ---------
        start (optional):   string - Timestamp of format '%Y-%m-%d %H:%M:%S.%f' or '%Y-%m-%d %H:%M:%S'
        end (optional):     string ... as start
        copy (optional):    bool - if False the instance will be changed. Else, a copy is returned

        Returns
        -------
        If copy is True:  housekeeping instance
        else:             nothing (instance is changed in place)


        Example
        -------
        >>> from atmPy.aerosols.instruments.piccolo import piccolo
        >>> launch = '2015-04-19 08:20:22'
        >>> landing = '2015-04-19 10:29:22'
        >>> hk = piccolo.read_file(filename) # create housekeeping instance
        >>> hk_zoom = zoom_time(hk, start = launch, end= landing)
        """

        if copy:
            housek = self.copy()
        else:
            housek = self

        if start:
            start = _time_tools.string2timestamp(start)
        if end:
            end = _time_tools.string2timestamp(end)

        try:
            housek.data = housek.data.truncate(before=start, after=end)
        except KeyError:
            txt = '''This error is most likely related to the fact that the index of the timeseries is not in order.
                  Run the sort_index() attribute of the DataFrame'''
            raise KeyError(txt)

        housek._start_time = housek.data.index[0]

        if copy:
            return housek
        else:
            return






    # def plot_versus_altitude_all(self):
    # axes = []
    #     for key in self.data.keys():
    #         f, a = plt.subplots()
    #         a.plot(self.data.altitude, self.data[key], label=key)
    #         a.legend()
    #         a.grid(True)
    #         axes.append(a)
    #     return axes

    def get_timespan(self, verbose = False):
        """
        Returns the first and last value of the index, which should be the first and last timestamp

        Returns
        -------
        tuple of timestamps
        """
        start = self.data.index[0]
        end = self.data.index[-1]
        if verbose:
            print('start: %s' % start.strftime('%Y-%m-%d %H:%M:%S.%f'))
            print('end:   %s' % end.strftime('%Y-%m-%d %H:%M:%S.%f'))
        return start, end

    def save(self, fname):
        """currently this simply saves the data of the timeseries using pandas
        to_csv

        Arguments
        ---------
        fname: str.
            Path to the file."""

        self.data.to_csv(fname)

    save_netCDF = save_netCDF


class TimeSeries_2D(TimeSeries):
    """
    experimental!!
    inherits TimeSeries

    differences:
        plotting
    """
    def __init__(self, *args):
        super(TimeSeries_2D,self).__init__(*args)

    def plot(self, xaxis = 0, ax = None, autofmt_xdate = True, cb_kwargs = {}, pc_kwargs = {},  **kwargs):
        if 'cb_kwargs' in kwargs.keys():
            cb_kwargs = kwargs['cb_kwargs']
        if 'pc_kwargs' in kwargs.keys():
            pc_kwargs = pc_kwargs
        f, a, pc, cb = _pandas_tools.plot_dataframe_meshgrid(self.data, xaxis=xaxis, ax=ax, pc_kwargs=pc_kwargs, cb_kwargs=cb_kwargs)
        if autofmt_xdate:
            f.autofmt_xdate()
        return f, a, pc, cb


class TimeSeries_3D(TimeSeries):
    """
    experimental!!
    inherits TimeSeries

    differences:
        plotting
    """
    def __init__(self, *args):
        super(TimeSeries_3D,self).__init__(*args)


    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, data):
        if not type(data).__name__ == 'Panel':
            raise TypeError('Data has to be of type DataFrame. It currently is of type: %s'%type(data).__name__)
        self.__data = data

    def plot(self, xaxis = 0, yaxis = 1, sub_set = 0, ax = None, kwargs = {}):

        f,a,pc,cb =  _pandas_tools.plot_panel_meshgrid(self.data, xaxis = xaxis,
                                                       yaxis = yaxis,
                                                       sub_set = sub_set,
                                                       ax = ax,
                                                       kwargs = kwargs)
        return f,a,pc,cb



# Todo: revive following as needed
# def get_sun_position(self):
#     """read docstring of solar.get_sun_position_TS"""
#     out = solar.get_sun_position_TS(self)
#     return out
#
# def convert2verticalprofile(self):
#     hk_tmp = self.copy()
#     hk_tmp.data['TimeUTC'] = hk_tmp.data.index
#     hk_tmp.data.index = hk_tmp.data.Altitude
#     return hk_tmp
#
# def plot_map(self, resolution = 'c', three_d=False):
#     """Plots a map of the flight path
#
#     Note
#     ----
#     packages: matplotlib-basemap,
#
#     Arguments
#     ---------
#     three_d: bool.
#         If flight path is plotted in 3D. unfortunately this does not work very well (only costlines)
#     """
#
#     data = self.data.copy()
#     data = data.loc[:,['Lon','Lat']]
#     data = data.dropna()
#
#     lon_center = (data.Lon.values.max() + data.Lon.values.min()) / 2.
#     lat_center = (data.Lat.values.max() + data.Lat.values.min()) / 2.
#
#     points = np.array([data.Lat.values, data.Lon.values]).transpose()
#     distances_from_center_lat = np.zeros(points.shape[0])
#     distances_from_center_lon = np.zeros(points.shape[0])
#     for e, p in enumerate(points):
#         distances_from_center_lat[e] = vincenty(p, (lat_center, p[1])).m
#         distances_from_center_lon[e] = vincenty(p, (p[0], lon_center)).m
#
#     lat_radius = distances_from_center_lat.max()
#     lon_radius = distances_from_center_lon.max()
#     scale = 1
#     border = scale * 2 * np.array([lat_radius, lon_radius]).max()
#
#     height = border + lat_radius
#     width = border + lon_radius
#     if not three_d:
#         bmap = Basemap(projection='aeqd',
#                        lat_0=lat_center,
#                        lon_0=lon_center,
#                        width=width,
#                        height=height,
#                        resolution=resolution)
#
#         # Fill the globe with a blue color
#         wcal = np.array([161., 190., 255.]) / 255.
#         boundary = bmap.drawmapboundary(fill_color=wcal)
#
#         grau = 0.9
#         continents = bmap.fillcontinents(color=[grau, grau, grau], lake_color=wcal)
#         costlines = bmap.drawcoastlines()
#         x, y = bmap(data.Lon.values, data.Lat.values)
#         path = bmap.plot(x, y,
#                          color='m')
#         return bmap
#
#     else:
#         bmap = Basemap(projection='aeqd',
#                    lat_0=lat_center,
#                    lon_0=lon_center,
#                    width=width,
#                    height=height,
#                    resolution=resolution)
#
#         fig = plt.figure()
#         ax = Axes3D(fig)
#         ax.add_collection3d(bmap.drawcoastlines())
#         x, y = bmap(self.data.Lon.values, self.data.Lat.values)
#         # ax.plot(x, y,self.data.Altitude.values,
#         #           color='m')
#         N = len(x)
#         for i in range(N - 1):
#             color = plt.cm.jet(i / N)
#             ax.plot(x[i:i + 2], y[i:i + 2], self.data.Altitude.values[i:i + 2],
#                     color=color)
#         return bmap, ax
#
#
# def plot_versus_pressure_sep_axes(self, what):
#     what = self.data[what]
#
#     ax = self.data.barometric_pressure.plot()
#     ax.legend()
#     ax.set_ylabel('Pressure (hPa)')
#
#     ax2 = ax.twinx()
#     what.plot(ax=ax2)
#     g = ax2.get_lines()[0]
#     g.set_color('red')
#     ax2.legend(loc=4)
#     return ax, ax2
#
# def plot_versus_pressure(self, what, ax=False):
#     what = self.data[what]
#
#     if ax:
#         a = ax
#     else:
#         f, a = plt.subplots()
#     a.plot(self.data.barometric_pressure.values, what)
#     a.set_xlabel('Barometric pressure (mbar)')
#
#     return a
#
#
#
# def plot_versus_altitude_sep_axes(self, what):
#     what = self.data[what]
#
#     ax = self.data.altitude.plot()
#     ax.legend()
#     ax.set_ylabel('Altitude (m)')
#
#     ax2 = ax.twinx()
#     what.plot(ax=ax2)
#     g = ax2.get_lines()[0]
#     g.set_color('red')
#     ax2.legend(loc=4)
#     return ax, ax2
#
# def plot_versus_altitude(self, what, ax=False, figsize=None):
#     """ Plots selected columns versus altitude
#
#     Arguments
#     ---------
#     what: {'all', key, list of keys}
#
#     Returns
#     -------
#     matplotlib.axes instance
#     """
#
#     allowed = ['altitude', 'Altitude', 'Height']
#
#     if what == 'all':
#         what = self.data.keys()
#
#     found = False
#     for i in allowed:
#         try:
#             x = self.data[i]
#             found = True
#             # print('found %s'%i)
#             break
#         except KeyError:
#             continue
#
#     if not found:
#         txt = 'TimeSeries instance has no attribute associated with altitude (%s)' % allowed
#         raise AttributeError(txt)
#     if ax:
#         f = ax[0].get_figure()
#     else:
#         f, ax = plt.subplots(len(what), sharex=True, gridspec_kw={'hspace': 0.1})
#         if len(what) == 1:
#             ax = [ax]
#
#     if not figsize:
#         f.set_figheight(4 * len(what))
#     else:
#         f.set_size_inches(figsize)
#
#     for e, a in enumerate(ax):
#         a.plot(x, self.data[what[e]], label=what[e])
#         a.legend()
#         a.set_xlabel('Altitude')
#
#     return ax
