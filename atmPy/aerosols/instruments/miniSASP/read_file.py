import numpy as _np
import pandas as _pd
import os as _os
import warnings as _warnings
from scipy import stats as _stats
from atmPy.aerosols.instruments.miniSASP import _miniSASP
import datetime as _datetime


def _simplefill(series):
    """Very simple function to fill missing values. Should only be used for
    values which basically do not change like month and day.
    Will most likely give strange results when the day does change.
    Returns: nothing everything happens inplace"""

    series.values[0] = series.dropna().values[0]
    series.values[-1] = series.dropna().values[-1]
    series.fillna(method='ffill', inplace=True)
    return


def _extrapolate(x, y):
    """This is a very simple extrapolation.
    Takes: two series of the same pandas DataTable. x is most likely the index of the DataTable
    assumption:
        - x is are the sampling points while y is the dataset which is incomplete and needs to be extrapolated
        - the relation ship is very close to linear
    proceedure:
        - takes the fist two points of y and performs a linear fit. This fit is then used to calculate y at the very first
          x value
        - similar at the end just with the last two points.
    returns: nothing. everthing happens inplace
    """

    xAtYnotNan = x.values[~_np.isnan(y.values)][:2]
    YnotNan = y.values[~_np.isnan(y.values)][:2]
    slope, intercept, r_value, p_value, slope_std_error = _stats.linregress(xAtYnotNan, YnotNan)

    fkt = lambda x: intercept + (slope * x)
    y.values[0] = fkt(x.values[0])

    xAtYnotNan = x.values[~_np.isnan(y.values)][-2:]
    YnotNan = y.values[~_np.isnan(y.values)][-2:]
    slope, intercept, r_value, p_value, slope_std_error = _stats.linregress(xAtYnotNan, YnotNan)

    fkt = lambda x: intercept + (slope * x)
    y.values[-1] = fkt(x.values[-1])

    return

def read_csv(fname, version = 'current', year = 2015, verbose=False):
    """Creates a single ULR instance from one file or a list of files.

    Arguments
    ---------
    fname: string or list
    version: str
        0.1: files till 2016-07-18 ... this includes Svalbard data
        current: files since 2016-07-18
    """
    if type(fname).__name__ == 'list':
        first = True
        for file in fname:
            if _os.path.split(file)[-1][0] != 'r':
                continue
            if verbose:
                print(file)
            # ulrt = miniSASP(file, verbose=verbose)
            ulrt = _process(file, version = version, year = year, verbose = verbose)
            if first:
                ulr = ulrt
                first = False
            else:
                ulr.data = _pd.concat((ulr.data, ulrt.data))
    else:
        # ulr = miniSASP(fname, verbose=verbose)
        ulr = _process(fname, version=version, year = year, verbose=verbose)

    return ulr

def _process(fname,
             version = 'current',
             year = 2015,
             create_timestamp = True,
             remove_data_withoutGPS = True,
             sort = True,
             remove_unused_cols = True,
             classify = True,
             verbose = False):
    df = _read_file(fname)
    df = _recover_negative_values(df)
    df = _norm2integration_time(df)
    if create_timestamp:
        df = _create_timestamp(df, version = version, year = year, verbose = verbose)
    if sort:
        df = df.sort_index()
    if remove_data_withoutGPS:
        df = _remove_data_withoutGPS(df)
    if remove_unused_cols:
        df = _remove_unused_columns(df)
    if classify:
        df = _miniSASP.MiniSASP(df)
    return df

def _remove_unused_columns(df):
    dropthis = ['var1', 'var2', 'var3', 'time', 'GPSHr', 'MonthDay', 'YearMonth',
                'Month', 'Day', 'HKSeconds', 'Modeflag', 'PhotoOffArr', 'Year',
                'GPSReadSeconds', 'GateLgArr', 'GateShArr', 'caseflag', 'homep', 'MicroUsed', ]
    for dt in dropthis:
        df.drop(dt, axis = 1, inplace = True)
    return df

def _remove_data_withoutGPS(df, day='08', month='01'):
    """ Removes data from before the GPS is fully initiallized. At that time the Date should be the 8th of January.
    This is an arbitray value, which might change

    Arguments
    ---------
    day (optional): string of 2 digit integer
    month (optional): string of 2 digit integer
    """
    df = df[((df.Day != day) & (df.Month != month))]
    return df

def _create_timestamp(df, version= 'current', year = 2015, millisscale = 10, verbose = False):
    if verbose:
        print(('================\n'
               'create timestamp\n'
               '-------\n'))
    df.Seconds *= (millisscale/1000.)
    df.index = df.Seconds

    df.YearMonth = _np.floor(df.MonthDay / 100.)
    _simplefill(df.YearMonth)

    df.Day = df.MonthDay - df.YearMonth*100
    _simplefill(df.Day)
    if version == 'current': #since 2016-07-18
        df.Year = '20' + df.YearMonth.astype(int).apply(lambda x: "{:.2}".format(str(x)))
        year_tmp = df.Year
        df.Month = df.YearMonth.astype(int).apply(lambda x: "{:2}".format(str(x)[2:]))
    elif version == '0.1': #older than 2016-07-18
        if verbose:
            print('version 0.1')
        year_tmp = str(year)
        df.Month = df.YearMonth.astype(int).apply(lambda x: "{:0>2}".format(str(x)))
    df.Day = df.Day.astype(int).apply(lambda x: '{0:0>2}'.format(x))
    # GPSHr_P_3 = df.GPSHr.copy()
    # Month_P_1 = df.Month.copy()
    # Day_P_1 = df.Day.copy()

    GPSunique = df.GPSHr.dropna().unique()
    for e,i in enumerate(GPSunique):
        where = _np.where(df.GPSHr == i)[0][1:]
        df.GPSHr.values[where] = _np.nan

    # GPSHr_P_1 = df.GPSHr.copy()

    # extrapolate and interpolate the time
    _extrapolate(df.index, df.GPSHr)
    df['GPSHr_p'] = df.GPSHr.copy()
    df.GPSHr.interpolate(method='index', inplace=True)
    # return df
    df.GPSHr.dropna(inplace=True)
    # GPSHr_P_2 = df.GPSHr.copy()
    # def GPSHr2timestr(x):
    #     h = x
    #     m = 60 * (x % 1)
    #     s = round(60 * ((60 * (x % 1)) % 1))
    #     if s >= 60.:
    #         s = 0.
    #         m += 1.
    #     if m >= 60.:
    #         m = 0.
    #         h += 1
    #     time_str = '%02i:%02i:%09.6f' % (h, m, s)
    #     return time_str
    def GPSHr2timestr(x):
        # looks complicated ... the only way I could make it work though
        pdt = _pd.Timedelta(x, 'h')
        pdt = pdt.to_pytimedelta()
        pdt = pdt + _datetime.datetime(2000, 1, 1)
        pdt_str = '{0:%H:%M:%S.%f}'.format(pdt)
        return pdt_str
    # return df
    df.GPSHr = df.GPSHr.apply(GPSHr2timestr)
    # return df
    ###### DateTime!!
    dateTime = year_tmp + '-' +  df.Month + '-' + df.Day +' ' + df.GPSHr
    # return df
    df['Time_new'] = _pd.to_datetime(dateTime, format="%Y-%m-%d %H:%M:%S.%f")
    # return df
    df.index = _pd.Series(_pd.to_datetime(dateTime, format="%Y-%m-%d %H:%M:%S.%f"), name='Time_UTC')

    df = df[_pd.notnull(df.index)]  # gets rid of NaT
    return df

def _norm2integration_time(df):
    columns = ['PhotoA', 'PhotoB', 'PhotoC', 'PhotoD', 'PhotoAsh', 'PhotoBsh', 'PhotoCsh', 'PhotoDsh']
    for col in columns:
        if 'sh' in col:
            df[col] = df[col] / df.GateShArr
        else:
            df[col] = df[col] / df.GateLgArr
    return df

def _recover_negative_values(df):
    columns = [df.PhotoA, df.PhotoB, df.PhotoC, df.PhotoD, df.PhotoAsh, df.PhotoBsh, df.PhotoCsh, df.PhotoDsh]
    do_warn = 0
    for col in columns:
        # print(col)
        series = col.values
        where = _np.where(series > 2 ** 16)
        series[where] = (series[where] * 2 ** 16) / 2 ** 16
        if where[0].shape[0] > 0:
            do_warn = where[0].shape[0]
            # print(do_warn)

    if _np.any(do_warn):
        #     print(do_warn)
        _warnings.warn("""This has to be checked!!! Dont know if i implemented it correctly! Arduino negatives become very large positive (unsigned longs)
       ;  recover using deliberate overflow""")

    return df

def _read_file(fname):
    df = _pd.read_csv(fname,
                      encoding="ISO-8859-1",
                      skiprows=16,
                      header=None,
                      error_bad_lines=False,
                      warn_bad_lines=False
                      )

    #### set column labels
    collabels = ['PhotoAsh',
                 'PhotoBsh',
                 'PhotoCsh',
                 'PhotoDsh',
                 'PhotoA',
                 'PhotoB',
                 'PhotoC',
                 'PhotoD',
                 'Seconds',
                 'caseflag',
                 'var1',
                 'var2',
                 'var3']

    df.columns = collabels

    #### Drop all lines which lead to errors
    # df = df.convert_objects(convert_numeric='force') # deprecated
    df = df.apply(_pd.to_numeric, errors = 'coerce')
    df = df.dropna(subset=['Seconds'])
    # self.data = df
    # return
    df = df.astype(float)
    # self.data = df
    # return

    #### add extra columns
    df['time'] = _np.nan
    df['azimuth'] = _np.nan
    df['homep'] = _np.nan
    df['MicroUsed'] = _np.nan
    df['lat'] = _np.nan
    df['lon'] = _np.nan
    df['Te'] = _np.nan
    df['GPSHr'] = _np.nan
    df['MonthDay'] = _np.nan
    df['YearMonth'] = _np.nan
    df['Year'] = _np.nan
    df['Month'] = _np.nan
    df['Day'] = _np.nan
    df['GPSReadSeconds'] = _np.nan
    df['HKSeconds'] = _np.nan
    df['Yaw'] = _np.nan
    df['Pitch'] = _np.nan
    df['Roll'] = _np.nan
    df['BaromPr'] = _np.nan
    df['BaromTe'] = _np.nan
    df['Modeflag'] = _np.nan
    df['GateLgArr'] = _np.nan
    df['GateShArr'] = _np.nan
    df['PhotoOffArr'] = _np.nan

    ##### Case 0
    case = _np.where(df.caseflag.values == 0)
    df.azimuth.values[case] = df.var1.values[case]
    df.homep.values[case] = df.var2.values[case]
    df.MicroUsed.values[case] = df.var3.values[case]

    ##### Case 1
    case = _np.where(df.caseflag.values == 1)
    df.lat.values[case] = df.var1.values[case]
    df.lon.values[case] = df.var2.values[case]
    df.Te.values[case] = df.var3.values[case]

    ##### Case 2
    case = _np.where(df.caseflag.values == 2)
    df.GPSHr.values[case] = df.var1.values[case]
    df.MonthDay.values[case] = df.var2.values[case]
    df.GPSReadSeconds.values[case] = df.var3.values[case].astype(float) / 100.
    df.HKSeconds.values[case] = df.Seconds.values[case]

    ##### Case 3
    case = _np.where(df.caseflag.values == 3)
    df.Yaw.values[case] = df.var1.values[case]
    df.Pitch.values[case] = df.var2.values[case]
    df.Roll.values[case] = df.var3.values[case]

    ##### Case 4
    case = _np.where(df.caseflag.values == 4)
    df.BaromPr.values[case] = df.var1.values[case]
    df.BaromTe.values[case] = df.var2.values[case]
    df.Modeflag.values[case] = df.var3.values[case].astype(float) + 0.5

    ##### Case 5
    case = _np.where(df.caseflag.values == 5)
    df.GateLgArr.values[case] = df.var1.values[case]
    df.GateShArr.values[case] = df.var2.values[case]
    df.PhotoOffArr.values[case] = df.var3.values[case]
    _simplefill(df.GateLgArr)
    _simplefill(df.GateShArr)
    return df