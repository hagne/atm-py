__author__ = 'htelg'
import warnings as _warnings
from copy import deepcopy as _deepcopy
import atmPy.general.vertical_profile 

import pandas as _pd
import numpy as _np
import matplotlib.pylab as _plt

from atmPy.tools import pandas_tools as _pandas_tools
from atmPy.tools import time_tools as _time_tools
from atmPy.tools import array_tools as _array_tools
from atmPy.tools import plt_tools as _plt_tools
from atmPy.tools import git as _git_tools
from atmPy.general import data_structure

try:
    from netCDF4 import Dataset as _Dataset
    # from netCDF4 import num2date as _num2date
    from netCDF4 import date2num as _date2num
except ModuleNotFoundError:
    _warnings.warn('netCDF4 not installed. You might encounter some functionality limitations.')
import warnings as _warnings
import datetime
from matplotlib.ticker import FuncFormatter as _FuncFormatter
from matplotlib.dates import DayLocator as _DayLocator
from matplotlib.dates import MonthLocator as _MonthLocator
import os as _os
from matplotlib import dates as _dates
import atmPy.general.statistics as _statistics

from atmPy.atmosphere import standards as _atm_std
from atmPy.atmosphere import atmosphere as _atmosphere


_unit_time = 'days since 1900-01-01'

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

    ts_time_num = _date2num(ts.data.index.to_pydatetime(), _unit_time)#.astype(float)
    time_var = ni.createVariable('time', ts_time_num.dtype, 'time')
    time_var[:] = ts_time_num
    time_var.units = 'days since 1900-01-01'

    var_data = ni.createVariable('data', ts.data.values.dtype, ('time', 'data_columns'))
    var_data[:] = ts.data.values

    ts_columns = ts.data.columns.values.astype(str)
    var_data_collumns = ni.createVariable('data_columns', ts_columns.dtype, 'data_columns')
    var_data_collumns[:] = ts_columns

    ni._type = type(ts).__name__
    ni._data_period = none2nan(ts._data_period)
    ni._x_label = none2nan(ts._x_label)
    ni._y_label =  none2nan(ts._y_label)
    ni.info = none2nan(ts.info)
    ni._atm_py_commit = _git_tools.current_commit()

    if leave_open:
        return ni
    else:
        ni.close()

#### Tools

def close_gaps(ts, verbose = False):
    """This is an older version to deal with gaps ... rather consider using the ones in the data_structure module"""
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


def align_to(ts, ts_other, how = 'linear_avg', tolerance = (5,'m'), verbose= False):
    """Similar to merge but without having all columns together. Also the function is more sufisticated with respect
    to up an downsampling (a rolling mean is used ... see how argument). Two different alignment methods are
    distinguished (see how argument). One that interpolates and averages, which assumes the data to be interpolated
    over the time periods. And that finds the closest value, which should be used if the the data is measured exactly
    at the given Timestamp.


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
    how: string ['linear_avg', 'find_closest']
        linear_avg: Align the TimeSeries ts to another by interpolating (linearly). If data periods differe by at least
            a factor of 2 a rolling mean is calculated with a window size equal to the ratio (if ratio is positive !!!).
        find_closest: the closest index in the timestamp is used as long as it is within the tolerance (see tolerance
            kwarg). This can be see as down sampling of a higher frequency data to a lower frequency data, by picking
            the closest value rather then interpolating as done by the align function.
    tolerance: tuple of value and unit
        Ignored if how is not find_closest. This is the maximum time difference that is allowed when finding the closest
        value. For details see pandas.Timedelta.

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

    if how == 'linear_avg':
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
    elif how == 'find_closest':
        ts_other.data = ts_other.data.iloc[:,[]]
        ts_other.data.columns.name = None
        # ts.data.columns = ['dummy']
        tsrm = merge(ts_other, ts, how = how, tolerance = tolerance)

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


def merge(ts, ts_other, how = 'interpolate_linear', tolerance = (5,'m'), recognize_gaps = True, verbose = False):
    """ Merges current with other timeseries. The returned timeseries has the same time-axes as the current
    one (as opposed to the one merged into it).
    There are two options on how the timeseries are merged. Either missing or offset data points are linearly
    interpolated. Or the closest value is found in a certein tolerance (see how kewarg)

    Argument
    --------
    ts: the other time series will be merged to this, therefore this timeseries
        will define the time stamps.
    ts_other: timeseries or one of its subclasses.
        List of TimeSeries objects.
    how: str
        interpolate_linear: missing or offset data is linearly interpolated. This can be seen as an upsampling of a
            lower frequency data to a higher frequency data.
        find_closest: the closest index in the timestamp is used as long as it is within the tolerance (see tolerance
            kwarg). This can be see as down sampling of a higher frequency data to a lower frequency data, by picking
            the closest value rather then interpolating as done by the align function.
    tolerance: tuple of value and unit
        Ignored if how is not find_closest. This is the maximum time difference that is allowed when finding the closest
        value. For details see pandas.Timedelta.
    recognize_gaps: bool
        ignored if how is not interpolate_linear. If ts_other has data gaps (missing) setting recognize_gaps to True
        will result in an attempt to detect the gaps and fill them with nans. False will lead to a linear interpolation
        of any type of gap (this includs streches of nan)

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
        if how == 'interpolate_linear':
            ts_data_list = [ts_this.data, ts_other.data]
            if not recognize_gaps:
                catsortinterp = _pd.concat(ts_data_list).sort_index().interpolate(method='index')
                merged = catsortinterp.groupby(catsortinterp.index).mean().reindex(ts_data_list[0].index)
                ts_this.data = merged

                # ts_data_list = [data.data, correlant_nocol]

            else:
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
        elif how == 'find_closest':
            merged = _pd.merge_asof(ts_this.data, ts_other.data, left_index=True, right_index=True,
                                   tolerance=_pd.Timedelta(*tolerance),
                                   #               on = 500, right_on=498.6,
                                   direction='nearest'
                                   )
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


def correlate(data, correlant, data_column = False, correlant_column = False, differenciate = None, remove_zeros=True, data_lim = None,
              correlant_lim = None,
              align_timeseries = True,
              weights = 'relative'):
    """Correlates data in correlant to that in data. In the process the data in correlant
    will be aligned to that in data. Make sure that data has the lower period (less data per period of time).

    Args:
        data:
        correlant:
        data_column:
        correlant_column:
        remove_zeros:
        data_lim: tuple
            lower and upper limit of data values
        correlant_lim:
            lower and upper limit of correlant values
        weights: 'str' ([relative], absolute)
            odr only!
            relative: the weights will scale with the x-data. Since we usually assume uncertainties to be relative to
                    the data rather then absolute this is the standard setting
            absolute:   bsically this means not weighted ... in principle you could apply different for x and y ... not
                        implemented yet

    Returns:

    """
    data = data.copy()
    correlant = correlant.copy()
    assert(type(differenciate).__name__ in ['Series', 'NoneType'])
    if type(differenciate).__name__ == 'Series':
        differenciate = differenciate.values
        if any([align_timeseries, data_lim, correlant_lim]):
            raise ValueError('Sorry differencieate and any of align_timeseries, data_lim, or correlant_lim do not work at the same time .... programming required')


    if data_column:
        data_values = data.data[data_column].values
    elif data.data.shape[1] > 1:
        raise ValueError('Data contains more than 1 column. Specify which to correlate. Options: %s'%(list(data.data.keys())))
    else:
        data_values = data.data.iloc[:,0].values

    if align_timeseries:
        correlant_aligned = correlant.align_to(data)
    else:
        correlant_aligned = correlant

    if correlant_column:
        correlant_values = correlant_aligned.data[correlant_column].values
    elif correlant.data.shape[1] > 1:
        raise ValueError('''Correlant contains more than 1 column. Specify which to correlate. Options:
%s'''%(list(correlant_aligned.data.keys())))
    else:
        correlant_values = correlant_aligned.data.iloc[:,0].values


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
    out = _array_tools.Correlation(data_values, correlant_values, differenciate, remove_zeros=remove_zeros, index = data.data.index, weights = weights)
    out._x_label_orig = 'DataTime'
    return out

def rolling_correlation_old(data, correlant, window, data_column = False, correlant_column = False,  min_good_ratio = 0.67, verbose = True):
    """This is an old version ... better use rolling(...).corr(...)
    time as here: http://docs.scipy.org/doc/numpy/reference/arrays.datetime.html#datetime-units"""

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
    # timestamps = _pd.DataFrame(_pd.to_datetime(_pd.Series(_np.zeros(size))))
    timestamps = _pd.to_datetime(_pd.Series(_np.zeros(size)))
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


def correlate_rolling(ts_data, ts_correlant,
                         data_column=None,
                         correlant_column=None,
                         window=(30, 'D'),
                         steps=30,
                         minvalidvalues=400,
                         what_to_calc=['r', 'odr'],
                         ):
    """As long as pandas still has problems to make a rolling correlation with a sliding time-window I came up with
    this. Rolling correlation of two timeseries.

    Parameters
    ----------
    ts_data: TimeSeries
    ts_correlant: TimeSeries
    data_column: str
        name of column in the first Timeseries to use in correlation.
    correlant_column: str
        name of column in the second Timeseries to use in correlation.
    window: tuple
        Value, unit
    steps: int
        To increase resolution: how many Timestamps within one window
    minvalidvalues: int
        When number of not-nan datapoints in window is below this value nan is returned
    what_to_calc: list of str
        r: Pearson r
        odr: results from a orthogonal distance regression

    Returns
    -------
    TimeSeries
    """


    # translate dictionary for frequency and date units
    freq2time = {'H': 'h',
                 'D': 'D'}

    # make one dataframe from the two inputs
    df = _pd.DataFrame()
    if data_column:
        df['data'] = ts_data.data.loc[:, data_column]
    else:
        df['data'] = ts_data.data.iloc[:,0]
    if correlant_column:
        df['correlant'] = ts_correlant.data.loc[:, correlant_column]
    else:
        df['correlant'] = ts_correlant.data.iloc[:, 0]

    # dropna because we want only to count valid datapoints
    df.dropna(inplace=True)

    # the function that is applied

    def fct(df):
        if df.size < minvalidvalues:
            return _np.nan
        else:
            return df.corr().values[0, 1]

    def get_odr(df):
        """function to retrieve orthogonal distance regression results"""
        if df.size < minvalidvalues:
            return _pd.Series((_np.nan, _np.nan), index=['c', 'm'])
        else:
            ts = TimeSeries(df)
            corr = ts.correlate_to(ts, data_column='data', correlant_column='correlant')
            odrout = corr.orthogonla_distance_regression['output']
            return _pd.Series(odrout.beta, index=['c', 'm'])

    df_list = []
    for e in range(steps):
        freq = '{}{}'.format(-1 * e * window[0] / (steps), window[1])
        dft = df.shift(freq=freq)

        # add a fixed value at window width below start of df to make the grouping is always at the same start
        dft.loc[df.index[0] - _pd.Timedelta(window[0], freq2time[window[1]])] = _np.nan
        dft.sort_index(inplace=True)

        # group and apply
        group = dft.groupby(_pd.Grouper(freq=window,
                                       label='left',
                                       ))
        out = group.apply(fct)

        dfout = _pd.DataFrame(out, columns=['r'])
        dfout_list = [dfout]
        if 'odr' in what_to_calc:
            dfout = group.apply(get_odr)
            dfout_list.append(dfout)

        dfout = _pd.concat(dfout_list, axis=1)
        # move label to center
        dfout.index += _pd.Timedelta(window[0] / 2, freq2time[window[1]])

        # adjust the time stamp according to the shift we applied
        freq = '{}{}'.format(e * window[0] / (steps), window[1])
        dfout = dfout.shift(freq=freq)

        df_list.append(dfout)

    # concat and sort the list of results to single dataframe
    out = _pd.concat(df_list)
    out.sort_index(inplace=True)
    out = TimeSeries(out)

    return out


def add_time_of_interest2axes(ax, times_of_interest):
    if type(times_of_interest) == dict:
        times_of_interest = [times_of_interest]
    for toi in times_of_interest:
        ts = toi.pop('datetime')
        if 'color' not in toi.keys():
            toi['color'] = 'black'
        try:
            annotate = toi.pop('annotate')
            annotate_kwargs = toi.pop('annotate_kwargs')
        except:
            annotate = None

        if 'vline_kwargs' not in toi.keys():
            toi['vline_kwargs'] = {}

        if 'color' not in toi['vline_kwargs'].keys():
            toi['vline_kwargs']['color'] = toi['color']

        ax.vlines(ts, -300, 40000, **toi['vline_kwargs'])

        if annotate:

            if 'bbox' not in annotate_kwargs:
                annotate_kwargs['bbox'] = dict(boxstyle="round,pad=0.3", fc=[1, 1, 1, 0.8], ec=toi['color'], lw=1)
            if 'ha' not in annotate_kwargs:
                annotate_kwargs['ha'] = 'center'
            # pos_y = annotate_kwargs.pop('pos_y')
            ax.annotate(annotate[0], (ts, annotate[1]), **annotate_kwargs)

def corr_timelag(ts, other, dt=(5, 'm'), no_of_steps=10, center=0, direction=None, min_good_ratio = 0, normalize=True, **kwargs):
    """
    Parameters
    ----------
    dt: tuple
        first arg of tuple can be int or array-like of dtype int. Second arg is unit. if array-like no_of... is ignored
    direction: bool or string
        if direction is set the center parameter will be ignored
        p for positive, n for negative
    min_good_ratio: float
        the minimum ratio of points that are not nan and the total number of data points. E.g. if there is a section of
        spotty data (lots of nan). The few resulting valid points would give a poor correlation which would be excluded"""


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
        dt_array = _np.arange(0, dt[0] * no_of_steps, dt[0])

        if direction:
            if direction == 'p':
                pass
            elif direction == 'n':
                dt_array *= -1
            else:
                txt = 'Direction has to be None, "p", or "n" ... it is %s' % direction
                raise ValueError(txt)

        else:
            dt_array += int(center) - int(no_of_steps * dt[0] / 2)

    # if center:
    #     dt_array += int(center)
    out = _pd.DataFrame(index = dt_array, columns=['pearson_r'], dtype=float)
    for dtt in dt_array:
        tst = other.copy()
        tst.data.index += _np.timedelta64(int(dtt), dt[1])
        corr = ts.correlate_to(tst)
        minshape = _np.array([tst.data.shape[0], ts.data.shape[0]]).min()
        if corr._data.shape[0] < (min_good_ratio * minshape):
            out.loc[dtt] = _np.nan
        else:
            out.loc[dtt] = corr.pearson_r[0]
        # if not out:
        #     out = corr
        #     import pdb
        #     pdb.set_trace()
        #     return out
        #     out.data.columns = [dtt]
        # else:
        #     out.data[dtt] = corr.data  # [self._data_column]

    # if normalize:
    #     out = TimeSeries_2D(out.data.apply(lambda line: line / line.max(), axis=1))
    # else:
    #     out = TimeSeries_2D(out.data)

    # def aaa(line):
    #     if (_np.isnan(line)).sum() > ((~_np.isnan(line)).sum() * 0.3):
    #         return _np.nan
    #     col_t = cols.copy()
    #     lt = line[~ _np.isnan(line)]
    #
    #     col_t = col_t[~ _np.isnan(line)]
    #     argmax = lt.argmax()
    #     realMax = col_t[argmax]
    #     return realMax
    #
    # cols = out.data.columns
    # dt_max = _np.apply_along_axis(aaa, 1, out.data.values)

    # cols = out.data.columns
    # dt_max = _np.apply_along_axis(lambda arr: cols[arr.argmax()], 1, out.data.values)
    # dt_max = TimeSeries(_pd.DataFrame(dt_max, index=out.data.index))
    # if dt[1] == 'm':
    #     ylt = 'min.'
    # else:
    #     ylt = dt[1]
    # dt_max._y_label = 'Time lag (%s)' % ylt
    return out

class Rolling(_pd.core.window.Rolling):
    def __init__(self, ts, window, min_valid_data_pts = None, min_good_ratio=None,
                 verbose=True, center = True, use_old= False,
                 **kwargs):
        """I wrote this because it was originally not included in pandas for time_windows. However, some of the
        functions below are now available in pandas and i updated the functions to use the pandas versions. The old
        functions should still be available though the argument use_old. For more documentation see also the doc string
        of pd.Dataframe.rolling().
        """
        if min_good_ratio and not use_old:
            raise ValueError('min_good_ration only works in combination with use_old, set min_valid_data_pts instead')
        elif min_valid_data_pts and use_old:
            raise ValueError('min_valid_data_pts only works if use_old is False, set min_good_ration instead.')

        self._use_old = use_old
        self.ts = ts
        if ts.data.columns.shape[0] == 1:
            self._data_column = ts.data.columns[0]
        else:
            txt = 'please make sure the timeseries has only one collumn'
            raise ValueError(txt)

        if self._use_old:
            window = _np.timedelta64(window[0], window[1])
            window = int(window / _np.timedelta64(int(ts._data_period), 's'))
            min_periods = int(window * min_good_ratio)


            if verbose:
                print('Each window contains %s data points of which at least %s are not nan.' % (window,
                                                                                                 min_periods))
            super().__init__(ts.data[self._data_column],
                             window,
                             min_periods=min_periods,
                             center = center,
                             **kwargs)
        else:
            super().__init__(ts.data[self._data_column],
                             _pd.to_timedelta(window[0], window[1]),
                             min_periods=min_valid_data_pts,
                             # center=center,
                             **kwargs)

    def corr(self, other, *args, **kwargs):
        raise AttributeError('Sorry this function is currently not working')
        if other.data.columns.shape[0] == 1:
            other_column = other.data.columns[0]
        else:
            txt = 'please make sure the timeseries has only one collumn'
            raise ValueError(txt)
        # other_column = 'Bs_G_Dry_1um_Neph3W_1'
        if self._use_old:
            other = other.align_to(self.ts)

        other = other.data[other_column]
        corr_res = super().corr(other, *args, **kwargs)
        corr_res_ts = TimeSeries(_pd.DataFrame(corr_res))
        corr_res_ts._data_period = self.ts._data_period
        return corr_res_ts

    def corr_timelag(self, other, dt=(5, 'm'), no_of_steps=10, center=0, direction = None, normalize=True, **kwargs):
        """
        Parameters
        ----------
        dt: tuple
            first arg of tuple can be int or array-like of dtype int. Second arg is unit. if array-like no_of... is ignored
        direction: bool or string
            if direction is set the center parameter will be ignored
            p for positive, n for negative"""


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
            dt_array = _np.arange(0, dt[0] * no_of_steps, dt[0])

            if direction:
                if direction == 'p':
                    pass
                elif direction == 'n':
                    dt_array *= -1
                else:
                    txt = 'Direction has to be None, "p", or "n" ... it is %s'%direction
                    raise ValueError(txt)

            else:
                dt_array += int(center) - int(no_of_steps * dt[0] / 2)

        # if center:
        #     dt_array += int(center)
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
            return realMax.astype(_np.float)

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

class WrappedPlot(list):
    def __init__(self, axes_list):
        for a in axes_list:
            self.append(a)


def plot_wrapped(ts,periods = 1, frequency = 'h', ylabel = 'auto', max_wraps = 10, ylim = None, ax = None, twin_x = None, skip_first = 0,  **plot_kwargs):
    """frequency: http://docs.scipy.org/doc/numpy/reference/arrays.datetime.html#datetime-units

    if ax is set, all other parameters will be ignored
    ylim: set to False if you don't want """
    if 'cb_kwargs' in plot_kwargs.keys():
        cb_kwargs = plot_kwargs.pop('cb_kwargs')
    else:
        cb_kwargs = False

    ax_blank = False
    if _np.any(ax):
        if hasattr(ax, '_periods'):
            periods = ax._periods
            frequency =  ax._frequency
            ylabel = ax._ylabel
            max_wraps = ax._max_wraps
            ylim_old = ax._ylim
            col_no = _np.array([len(a.get_lines()) for a in ax]).max()
        else:
            a = ax
            f = a[0].get_figure()
            ax_blank = True
            ax = None

    if periods >1:
        raise ValueError('Sorry periods larger one is not working ... consider fixing it?!?')
    start, end = ts.get_timespan()
    length = end - start

    if frequency == 'Y':
        start = _np.datetime64(str(start.year))
        end = _np.datetime64(str(end.year))
        periods_no = int((end - start + 1) / _np.timedelta64(1, 'Y'))
        autofmt_xdate = True
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
    if _np.any(ax):
        a = ax
        f = a[0].get_figure()
    else:
        if not ax_blank:
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
        else:
            col_no = 0 #not sure what that is good for ... but is needed ;-)
    bbox_props = dict(boxstyle="round,pad=0.3", fc=[1,1,1,0.4], ec="black", lw=0.5 * _plt.rcParams['axes.linewidth'])

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
                at.xaxis.set_major_locator(_DayLocator(interval=4))
                at.xaxis.set_major_formatter(df)

            elif frequency == 'Y':
                def month_formatter(x, pos):
                    month_str = ['XXX', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
                    dt = _np.datetime64('0') + _np.timedelta64(int(x), 'D')
                    dt = _pd.to_datetime(dt)
                    out = dt.month
                    return month_str[out]

                df = _FuncFormatter(month_formatter)
                at.xaxis.set_major_locator(_MonthLocator(interval=1))
                at.xaxis.set_major_formatter(df)
                text = text[0].split('-')[0]
            elif frequency == 'D':
                text = start_t
            else:
                text = 'not set'
            at.text(txtpos[0], txtpos[1], text, transform=at.transAxes, bbox=bbox_props)


        if tst:
            tst.data.index  = tst.data.index - _pd.to_datetime(start_t)
            tst.data.index = _pd.to_datetime(tst.data.index + _np.datetime64('1900'))
            if 'TimeSeries' in (tst.__class__.__bases__[0].__name__ , type(tst).__name__):
                if twin_x:
                    plt_out = tst.plot(ax=at, autofmt_xdate=autofmt_xdate,
                                       # color=_plt_tools.color_cycle[col_no],
                                       **plot_kwargs)
                else:
                    plt_out = tst.plot(ax=at, autofmt_xdate = autofmt_xdate,
                                       # color = _plt_tools.color_cycle[col_no],
                                       **plot_kwargs)
                    # plt_out = tst.plot(ax=at)
            elif 'TimeSeries_2D' in (tst.__class__.__bases__[0].__name__, type(tst).__name__):
                plt_out = tst.plot(ax=at, autofmt_xdate=autofmt_xdate, color=_plt_tools.color_cycle[col_no], cb_kwargs = False, **plot_kwargs)
                plt_out[2].set_clim(ylim)
            else:
                txt = 'Time series or parent of Time series must be of type "TimeSeries" or "TimeSeries_2D. It is %s and %s ,respectively.'%(tst.__class__.__bases__[0].__name__, type(tst).__name__)
                raise TypeError(txt)

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
    _plt_tools.set_shared_label(a,ylabel, axis = 'y')
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
    sampling_period: int
        This is the period that the data is roughly sampled at in Seconds. None will cause some operations to fail, e.g. align and merg!
    """

    def __init__(self, data, sampling_period='auto', info=None):
        # if not type(data).__name__ == 'DataFrame':
        #     raise TypeError('Data has to be of type DataFrame. It currently is of type: %s'%(type(data).__name__))

        self.data = data
        self.info = info
        self.statistics = _statistics.Statistics(self)

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

        self.data_structure = data_structure.DataStructure(self)
        if sampling_period == 'auto':
            self._data_period = self.data_structure.estimate_sampling_period()
        else:
            self._data_period = sampling_period


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
                out = other.data.truediv(self_t.data.iloc[:,0], axis = 0)
                out = 1/out

            elif _np.all(self.data.columns == other.data.columns):
                out = self.data.divide(other.data)
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
        elif _np.all(self.data.columns == other.data.columns):
            out = self.data.add(other.data)
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

        elif _np.all(self.data.columns == other.data.columns):
            out = self.data.subtract(other.data)

        else:
            txt = 'at least one of the dataframes have to have one column only'
            raise ValueError(txt)

        ts = TimeSeries(out)
        ts._data_period = self._data_period
        return ts

    def __mul__(self,other):
        self = self.copy()
        other = other.copy()

        if not _np.array_equal(self.data.index, other.data.index):
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
        elif _np.all(self.data.columns == other.data.columns):
            out = self.data.multiply(other.data)
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



    def convert2verticalprofile(self, altitude_column = 'Altitude', resolution = None, return_std = False):
        """Convertes the time series into a vertical profile based on a column containing altitude
        information. In its simplest form it replaces the index with the altitude column. If resolution
        is set the data will be binned into altitude bins.
        Arguments
        ---------
        altitude_column: str ['Altitude']
            column label which contains the altitude information
        resolution: float or tuple of len == 3
            if float: altitude resolution in the same unit as data in the altitude column
            if tuple: (resolution, lower_limit, upper_limit)
        return_std: bool
            If to return an additional vertical profile which contains the
            standard deviation rather than the mean of all values within a 
            altitude bin.
        """

        ts_tmp = self.copy()
    #     hk_tmp.data['Time'] = hk_tmp.data.index
    #     if alt_label:
    #         label = alt_label
    #     else:
    #         label = 'Altitude'
    #     if alt_timeseries:
    #         alt_timeseries = alt_timeseries.align_to(ts_tmp)
    #         _pandas_tools.ensure_column_exists(alt_timeseries.data, 'Altitude', col_alt=alt_label)
    #         ts_tmp.data.index = alt_timeseries.data['Altitude']
    #     else:
        _pandas_tools.ensure_column_exists(ts_tmp.data, altitude_column)
        ts_tmp.data['DateTime'] = ts_tmp.data.index
        ts_tmp.data.index = ts_tmp.data[altitude_column]

        if resolution:
            ts_tmp.data.sort_index(inplace=True)
            if type(resolution) == tuple:
                start = resolution[1]
                end = resolution[2]
                resolution = resolution[0]
            else:
                start = _np.floor(ts_tmp.data[altitude_column].min())
                end = _np.ceil(ts_tmp.data[altitude_column].max())
            vertical_bin_edges = _np.arange(start, end + 1, resolution)
            vertical_bin_edges = _np.array([vertical_bin_edges[0:-1], vertical_bin_edges[1:]]).transpose()
            index = _np.apply_along_axis(lambda x: x.sum(), 1, vertical_bin_edges) / 2.
            df = _pd.DataFrame(_np.zeros((vertical_bin_edges.shape[0], ts_tmp.data.shape[1])),
                              index=index, columns=ts_tmp.data.columns)
            if return_std:
                dfstd = df.copy()

            for l in vertical_bin_edges:
                where = _np.where(_np.logical_and(ts_tmp.data.index > l[0], ts_tmp.data.index < l[1]))[0]
                mean = ts_tmp.data.iloc[where, :].mean()
                df.loc[(l[0] + l[1]) / 2] = mean
                if return_std:
                    std = ts_tmp.data.iloc[where, :].std()
                    dfstd.loc[(l[0] + l[1]) / 2] = std

        else:
            df = ts_tmp.data

        out = atmPy.general.vertical_profile.VerticalProfile(df)
        out._x_label = self._y_label

        if return_std:
            out_std = atmPy.general.vertical_profile.VerticalProfile(dfstd)
            out_std._x_label = self._y_label
            return out, out_std
        else:
            return out


    def _del_all_columns_but(self, keep, inplace = False):
        """as it says, deletes all columns but ...
        Parameters:
        -----------
        keep: string or array-like
            column name(s) to keep
        """
        if inplace:
            ts = self
        else:
            ts = self.copy()
        all_keys = list(ts.data.keys())

        if isinstance(keep, str):
            keep = [keep]

        for k in keep:
            try:
                all_keys.remove(k)
            except ValueError:
                pass

        ts.data = ts.data.drop(labels=all_keys, axis=1)
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

    def average_time(self, window, std = False, envelope = False, verbose = False):
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

        if window[1] == 'm':
            window_rs = (window[0],'min')
        else:
            window_rs = window
        resample = ts.data.resample(window_rs,
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
        # ts._data_period = _pd.Timedelta(*window) / _pd.Timedelta(1,'s')
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

    corr_timelag = corr_timelag

    correlate_rolling = correlate_rolling

    merge = merge

    copy = _deepcopy

    # rollingR = Rolling

    def rolling(self, window,
                # data_column=False,
                #  correlant_column=False,
                use_old=False,
                min_valid_data_pts=None,
                min_good_ratio=None,
                verbose=True):
        """This will hopefully work one day. Currently as of (2018-08) the rolling and then corr will create
        unreasonable values, e.g. the result changes dramatically when with the length of the time series even if the
        length exceeds the window width by alot. Pleas use the rolling correlation instead .... it might be possible to
        bring that function into a similar structure as the rolling().corr() ... some time"""

        return Rolling(self, window,
                       # data_column=data_column,
                       # correlant_column=correlant_column,
                       use_old=use_old,
                       min_valid_data_pts=min_valid_data_pts,
                       min_good_ratio=min_good_ratio,
                       verbose=verbose)

    def remove_artefacts(self, which, sigma=1.5, window=3, verbose=False, inplace = False):
        """Removes artifacts by testing if the std is above a certain threshold"""
        data = self.data[which]
        roll = data.rolling(window=window, center=True)
        roll_std = roll.std()
        threshold = roll_std.mean() * sigma
        spike_pos = roll_std > threshold
        no_removed = spike_pos.sum()
        data_ft = data.copy()
        data_ft[spike_pos] = _np.nan
        if verbose:
            print("Number of points removed: {}".format(no_removed))
        if inplace:
            self.data[which] = data_ft
        return data_ft

    def plot(self, columns2plot=None, ax = None, label = None, autofmt_xdate = True, times_of_interest = None, plot_engine = 'pandas', **kwargs):
        """Plot each parameter separately versus time
        Arguments
        ---------
        columns2plot: column names or list of column names
            if only some of the columns ought to be plotted. Only effects plot_engine:matplotlib
        times_of_interest: dict or list of dicts
            excepted keys for each dict:
                datetime: e.g. '2017-04-02 23:50:00'
                annotate: (str, float)
                    This is the text and the yposition
                annotate_kwargs: dict of annotation kwargs
                vline_kwargs: dict of vline kwargs
        plot_engine: str, ('pandas', ['matplotlib'])
            Depending on the pandas version plotting through pandas (its still matplotlib of course) can resultin errors
            or be beneficial ... decide for yourself.
            Problems encountered 2019-11-25:
                I had an dataframe with nan everywhere but the very end and plotted with matplotlib
                -> it came up with a few timestamps that were not even present in the dataframe!!
        kwargs: keyword argurments passed to matplotlib plot function e.g.:
        picker: float
            in addition to the normal behaviour this will add text with the position to the plot and append the
            timeseries attribute plot_click_position. Matplotlib has to be in some kind of interactive backend.


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
        if 'picker' in kwargs.keys():
            self.plot_click_positions = []
            self.plot_events = []

            def onclick(event):
                self.plot_events.append(event)
                x = event.mouseevent.xdata
                xdate = _pd.Timestamp(_dates.num2date(x))
                y = event.mouseevent.ydata
                self.plot_click_positions.append([xdate, y])
                ax.text(x, y, 'x = {:02d}:{:02d}:{:02d}\ny = {:.0f}'.format(xdate.hour, xdate.minute, xdate.second, y))

            f.canvas.mpl_connect('pick_event', onclick)

        if plot_engine == 'pandas':
            self.data.plot(ax = ax, **kwargs)

        elif plot_engine == 'matplotlib':
            did_plot = False  # had to implement that since matploglib cept crashing when all where nan
            for e,k in enumerate(self.data.keys()):
                if not isinstance(columns2plot, type(None)):
                    if not isinstance(columns2plot, list):
                        columns2plot = [columns2plot,]
                    # if isinstance(columns2plot, str):
                    if k not in columns2plot:
                        continue
                    # else:
                    #     if e not in columns2plot:
                    #         continue

                if not label:
                    label_t = k
                else:
                    label_t = label

                if _np.all(_np.isnan(self.data[k].values)):
                    continue

                ax.plot(self.data.index.values, self.data[k].values, label = label_t, **kwargs)

                if self._time_format == 'timedelta':
                    formatter = _FuncFormatter(timeTicks)
                    ax.xaxis.set_major_formatter(formatter)

                did_plot = True

            if did_plot:
                if autofmt_xdate:
                    f.autofmt_xdate()

        ax.set_xlabel(self._x_label)
        ax.set_ylabel(self._y_label)




        if times_of_interest:
            add_time_of_interest2axes(ax, times_of_interest)
            # if type(times_of_interest) == dict:
            #     times_of_interest = [times_of_interest]
            # for toi in times_of_interest:
            #     ts = toi.pop('datetime')
            #     if 'color' not in toi.keys():
            #         toi['color'] = 'black'
            #     try:
            #         annotate = toi.pop('annotate')
            #         annotate_kwargs = toi.pop('annotate_kwargs')
            #     except:
            #         annotate = None
            #
            #     if 'vline_kwargs' not in toi.keys():
            #         toi['vline_kwargs'] ={}
            #
            #     ax.vlines(ts, -300, 40000, **toi['vline_kwargs'])
            #
            #     if annotate:
            #
            #         if 'bbox' not in annotate_kwargs:
            #             annotate_kwargs['bbox'] = dict(boxstyle="round,pad=0.3", fc=[1, 1, 1, 0.8], ec=toi['color'], lw=1)
            #         if 'ha' not in annotate_kwargs:
            #             annotate_kwargs['ha'] = 'center'
            #         # pos_y = annotate_kwargs.pop('pos_y')
            #         ax.annotate(annotate[0], (ts, annotate[1]), **annotate_kwargs)

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
        try:
            housek._start_time = housek.data.index[0]
        except IndexError:
            pass

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

    def save_csv(self, fname):
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
        super().__init__(*args)

    def plot(self, xaxis = 0, ax = None, autofmt_xdate = True, times_of_interest = None, cb_kwargs = None, pc_kwargs = None,  **kwargs):
        """

        Args:
            xaxis:
            ax:
            autofmt_xdate:
            times_of_interest:
            cb_kwargs:
            pc_kwargs:
            **kwargs:

        Returns:

        Examples:
            cmap = calipso.get_cmap(norm = 'log',log_min=-2.5, reverse=True)
            toi = {'datetime': cali.path.get_closest2location(loc_oliktok).index[0],
                   'linestyle': '--',
                   'color': 'black',
                   'annotate': 'Min. dist. to Oliktok ({} km)'.format(int(closest.distance[0])),
                   'annotate_kwargs': {'pos_y': 7000}}

            f,a,pc,cb = cali.total_attenuated_backscattering.plot(times_of_interest = [toi],
                                                                 pc_kwargs={'cmap': cmap,
                                                                            'vmin': 0,
                                                                            'vmax': 2
            #                                                                 'norm': plt.Normalize(vmin=0,vmax=1.)
                                                                           },
                                                                 cb_kwargs=True)
        """
        # if 'cb_kwargs' in kwargs.keys():
        #     cb_kwargs = kwargs['cb_kwargs']
        # if 'pc_kwargs' in kwargs.keys():
        #     pc_kwargs = pc_kwargs
        f, a, pc, cb = _pandas_tools.plot_dataframe_meshgrid(self.data, xaxis=xaxis, ax=ax, pc_kwargs=pc_kwargs, cb_kwargs=cb_kwargs)
        if autofmt_xdate:
            f.autofmt_xdate()

        if times_of_interest:
            add_time_of_interest2axes(a,times_of_interest)

        return f, a, pc, cb

    def _del_all_columns_but(self, keep, inplace = False):
        deled = super()._del_all_columns_but(keep, inplace = inplace)
        out = TimeSeries(deled.data)
        out._data_period = deled._data_period
        return out


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


def get_altitude(ts,
                 standard_atmosphere = False,
                 pressure_ref = None,
                 altitude_ref = None,
                 laps_rate = 0.0065,
                 col_name_press = 'Barometric_pressure',
                 col_name_temp = 'Temperature'):
    """Calculates the altitude from the measured barometric pressure
    Arguments
    ---------
    temperature: bool or array-like, optional
        False: temperature according to international standard is assumed.
        arraylike: actually measured temperature in Kelvin.

    Returns
    -------
    returns altitude and adds it to this instance
    """
    if standard_atmosphere:
        alt, tmp = _atm_std.standard_atmosphere(ts.data.loc[:, 'Barometric_pressure'], quantity='pressure')
    else:
        bf = _atmosphere.Barometric_Formula(pressure=ts.data.loc[:, col_name_press],
                                            pressure_ref= pressure_ref,
                                            alt_ref= altitude_ref,
                                            temp=ts.data.loc[:, col_name_temp],
                                            laps_rate=laps_rate)
        alt = bf.altitude

    ts.data['Altitude'] = alt
    return alt