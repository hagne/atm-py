from netCDF4 import Dataset as _Dataset
from netCDF4 import num2date as _num2date
import numpy as _np
import pandas as _pd
from .general import timeseries as _timeseries
from atmPy.aerosols.instruments.miniSASP import _miniSASP
from atmPy.general import vertical_profile as _vertical_profile
from atmPy.aerosols.instruments.POPS import housekeeping as _pops_hk
import atmPy.aerosols.size_distribution.sizedistribution as _sizedistribution
import xarray as _xr
import pathlib as pl

# import warnings as _warnings

importable_types = {#########
                    ### Time series
                    'TimeSeries':           {'call': _timeseries.TimeSeries, 'category': 'timeseries'},
                    'TimeSeries_2D':        {'call': _timeseries.TimeSeries_2D, 'category': 'timeseries'},
                    'TimeSeries_3D':        {'call': _timeseries.TimeSeries_3D, 'category': 'timeseries'},
                    'Sun_Intensities_TS':   {'call': _miniSASP.Sun_Intensities_TS, 'category': 'timeseries'},
                    'POPSHouseKeeping':     {'call': _pops_hk.POPSHouseKeeping, 'category': 'timeseries'},
                    ##########
                    ### Vertical profiles
                    'VerticalProfile':      {'call': _vertical_profile.VerticalProfile, 'category': 'verticalprofile'}
                    }



def files_between(root: pl.Path, start: _pd.Timestamp|str, end: _pd.Timestamp|str, 
                  globpattern: str = '*{date:%Y%m%d}*', 
                  input_directory_structure: str = "yearly", 
                  verbose = True,
                  return_date_instances = False
                  ):
    """ Generator that yields all files between start and end dates (inclusive) in the given root
    directory. Dates are determined from the filename. How the date is represented in the filename 
    needs to be specified in the globpattern kwarg.

    Parameters
    ----------
    root : pl.Path
        Root directory containing year subdirectories with files.
    start : pd.Timestamp|str
        Start date.
    end : pd.Timestamp|str
        End date.
    globpattern : str, optional
        Glob pattern to match files. How to contain the information how the date is formatted in the file names, e.g.  
        standard: tbl_20250513.nc -> '*{date:%Y%m%d}*', 
        y2k bug and julian date: tbl_25123.nc '*{date:%Y%m%d}*'
    Yields
    -------
    pl.Path
        Paths to files between start and end dates.
    """ 
    start = _pd.Timestamp(start)
    end = _pd.Timestamp(end)
    assert('{date:' in globpattern), f"globpattern ({globpattern}) has to define how timestamps (dates) are formated in the filename. E.g. '*{{date:%Y%m%d}}*' "
    assert(end > start), f'End must come after start! (end: {end}, start{start})'
    root = pl.Path(root)
    day = _pd.to_datetime(start.date())
    while day <= end:
        if input_directory_structure == 'yearly':
            year_dir = root / f"{day.year}"
        else:
            year_dir = root
        if not year_dir.exists():
            if verbose:
                print(f'The directory {year_dir} does not exists, be more carefull with you start, end, day arguments.')
        # yield from year_dir.glob(globpattern.format(date = d))
        if return_date_instances:
            yield from ((day, path) for path in year_dir.glob(globpattern.format(date = day)))
        else:
            yield from (path for path in year_dir.glob(globpattern.format(date = day)))
        day += _pd.to_timedelta(1, 'D')


def open_atmpy(fname):
    sd = _xr.open_dataset(fname)
    data_type = sd._atmPy.to_pandas().loc['type']
    if data_type in ['SizeDist_LS', 'SizeDist_TS', 'SizeDist']:
        out = _sizedistribution.open_netcdf(fname)
    else:
        raise ValueError('This file opener is not capable yet to open the data type {}. Try the other options!!'.format(data_type))

    return out

def open_netCDF(fname, data_type = None, error_unknown_type = True, verbose = False):
    """

    Parameters
    ----------
    fname
    data_type: string
        You can overwrite the type of the file. This is handy in case you get an error that the type is unkown

    Returns
    -------

    """

    ni = _Dataset(fname, 'r')

    # for very old files which do not set this attribute yet
    if not data_type:
        try:
            data_type = ni.getncattr('_type')
        except AttributeError:
            try:
                data_type = ni.getncattr('_ts_type')
            except AttributeError:
                txt = 'File has no attribute "_type". You can set kwarg data_type if you know what the type is supposed to be. E.g. data_type = TimeSeries'
                if error_unknown_type:
                    raise TypeError(txt)
                # _warnings.warn(txt)
                # data_type = 'TimeSeries'
                # print('Warning do not seam to be working ... hier a printout of the warning: %s'%txt)

        # also for older file types
        if not data_type:
            data_type = ni.getncattr('_ts_type')
            print('pos3', data_type)

    if verbose:
        print('Type is %s.'%data_type)
    # test if type among known types
    if data_type in importable_types.keys():
        category = importable_types[data_type]['category']

    else:
        txt = 'Type %s is unknown, programming required.' % data_type
        raise TypeError(txt)

    if category == 'timeseries':
        # load time
        time_var = ni.variables['time']
        # time_var.units
        ts_time = _num2date(time_var[:], time_var.units)
        index = _pd.DatetimeIndex(ts_time)

    elif category == 'verticalprofile':
        # load altitude
        alt_var = ni.variables['altitude']
        index = alt_var[:]


    # load  data
    var_data = ni.variables['data']
    ts_data = _pd.DataFrame(var_data[:], index=index)

    # load column names
    var_data_col = ni.variables['data_columns']
    ts_data = _pd.DataFrame(var_data[:], index=index,
                           columns=var_data_col[:])

    # test which type of timeseries (1D, 2D, 3D)


    # create time series
    # if ts_type == 'Sun_Intensities_TS':
    #     ts_out = _miniSASP.Sun_Intensities_TS(ts_data)
    # if type in importable_types.keys():
    ts_out = importable_types[data_type]['call'](ts_data)
    #
    # elif ts_type == 'TimeSeries':
    #     ts_out = TimeSeries(ts_data)
    # elif ts_type == 'TimeSeries_2D':
    #     ts_out = TimeSeries_2D(ts_data)
    # elif ts_type == 'TimeSeries_3D':
    #     ts_out = TimeSeries_3D(ts_data)
    # else:
    #     txt = 'Type %s is unknown, programming required.' % type
    #     raise TypeError(txt)

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