import numpy as np
import pandas as pd
import os
import warnings


def _simplefill(series):
    """Very simple function to fill missing values. Should only be used for
    values which basically do not change like month and day.
    Will most likely give strange results when the day does change.
    Returns: nothing everything happens inplace"""

    series.values[0] = series.dropna().values[0]
    series.values[-1] = series.dropna().values[-1]
    series.fillna(method='ffill', inplace=True)
    return

def read_csv(fname, version = 'current', verbose=False):
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
            if os.path.split(file)[-1][0] != 'r':
                continue
            if verbose:
                print(file)
            # ulrt = miniSASP(file, verbose=verbose)
            ulrt = _process(file, version = version, verbose = verbose)
            if first:
                ulr = ulrt
                first = False
            else:
                ulr.data = pd.concat((ulr.data, ulrt.data))
    else:
        # ulr = miniSASP(fname, verbose=verbose)
        ulr = _process(fname, version=version, verbose=verbose)
    ulr.data = ulr.data.sort_index()
    return ulr

def _process(fname, version = 'current', verbose = False):
    df = _read_file(fname)
    df = _recover_negative_values(df)
    return df


def _recover_negative_values(df):
    columns = [df.PhotoA, df.PhotoB, df.PhotoC, df.PhotoD, df.PhotoAsh, df.PhotoBsh, df.PhotoCsh, df.PhotoDsh]
    do_warn = 0
    for col in columns:
        print(col)
        series = col.values
        where = np.where(series > 2 ** 16)
        series[where] = (series[where] * 2 ** 16) / 2 ** 16
        if where[0].shape[0] > 0:
            do_warn = where[0].shape[0]
            print(do_warn)

    if np.any(do_warn):
        #     print(do_warn)
        warnings.warn("""This has to be checked!!! Dont know if i implemented it correctly! Arduino negatives become very large positive (unsigned longs)
       ;  recover using deliberate overflow""")

    return df

def _read_file(fname):
    df = pd.read_csv(fname,
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
    df = df.convert_objects(convert_numeric='force')
    df = df.dropna(subset=['Seconds'])
    # self.data = df
    # return
    df = df.astype(float)
    # self.data = df
    # return

    #### add extra columns
    df['time'] = np.nan
    df['azimuth'] = np.nan
    df['homep'] = np.nan
    df['MicroUsed'] = np.nan
    df['lat'] = np.nan
    df['lon'] = np.nan
    df['Te'] = np.nan
    df['GPSHr'] = np.nan
    df['MonthDay'] = np.nan
    df['Month'] = np.nan
    df['Day'] = np.nan
    df['GPSReadSeconds'] = np.nan
    df['HKSeconds'] = np.nan
    df['Yaw'] = np.nan
    df['Pitch'] = np.nan
    df['Roll'] = np.nan
    df['BaromPr'] = np.nan
    df['BaromTe'] = np.nan
    df['Modeflag'] = np.nan
    df['GateLgArr'] = np.nan
    df['GateShArr'] = np.nan
    df['PhotoOffArr'] = np.nan

    ##### Case 0
    case = np.where(df.caseflag.values == 0)
    df.azimuth.values[case] = df.var1.values[case]
    df.homep.values[case] = df.var2.values[case]
    df.MicroUsed.values[case] = df.var3.values[case]

    ##### Case 1
    case = np.where(df.caseflag.values == 1)
    df.lat.values[case] = df.var1.values[case]
    df.lon.values[case] = df.var2.values[case]
    df.Te.values[case] = df.var3.values[case]

    ##### Case 2
    case = np.where(df.caseflag.values == 2)
    df.GPSHr.values[case] = df.var1.values[case]
    df.MonthDay.values[case] = df.var2.values[case]
    df.GPSReadSeconds.values[case] = df.var3.values[case].astype(float) / 100.
    df.HKSeconds.values[case] = df.Seconds.values[case]

    ##### Case 3
    case = np.where(df.caseflag.values == 3)
    df.Yaw.values[case] = df.var1.values[case]
    df.Pitch.values[case] = df.var2.values[case]
    df.Roll.values[case] = df.var3.values[case]

    ##### Case 4
    case = np.where(df.caseflag.values == 4)
    df.BaromPr.values[case] = df.var1.values[case]
    df.BaromTe.values[case] = df.var2.values[case]
    df.Modeflag.values[case] = df.var3.values[case].astype(float) + 0.5

    ##### Case 5
    case = np.where(df.caseflag.values == 5)
    df.GateLgArr.values[case] = df.var1.values[case]
    df.GateShArr.values[case] = df.var2.values[case]
    df.PhotoOffArr.values[case] = df.var3.values[case]
    _simplefill(df.GateLgArr)
    _simplefill(df.GateShArr)
    return df