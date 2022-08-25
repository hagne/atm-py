# -*- coding: utf-8 -*-
"""
@author: Hagen Telg
"""

import numpy as np
import pandas as pd

from atmPy.general import timeseries
from atmPy.tools import time_tools


def _drop_some_columns(data):
    data.drop('Clock', axis=1, inplace=True)
    data.drop('Year', axis=1, inplace=True)
    data.drop('Month', axis=1, inplace=True)
    data.drop('Day', axis=1, inplace=True)
    data.drop('Hours', axis=1, inplace=True)
    data.drop('Minutes', axis=1, inplace=True)
    data.drop('Seconds', axis=1, inplace=True)


def _read_file(fname):
    picof = open(fname, 'r')
    header = picof.readline()
    picof.close()

    header = header.split(' ')
    header_cleaned = []

    for head in header:
        bla = head.replace('<', '').replace('>', '')
        where = bla.find('[')
        if where != -1:
            bla = bla[:where]
        header_cleaned.append(bla)

    data = pd.read_csv(fname,
                       names=header_cleaned,
                       sep=' ',
                       skiprows=1,
                       header=0)

    data.drop(range(20), inplace=True)  # dropping the first x lines, since the time is often dwrong

    time_series = data.Year.astype(str) + '-' + data.Month.apply(lambda x: '%02i' % x) + '-' + data.Day.apply(
        lambda x: '%02i' % x) + ' ' + data.Hours.apply(lambda x: '%02i' % x) + ':' + data.Minutes.apply(
        lambda x: '%02i' % x) + ':' + data.Seconds.apply(lambda x: '%05.2f' % x)
    data.index = pd.Series(pd.to_datetime(time_series, format=time_tools.get_time_formate()))

    _drop_some_columns(data)

    # convert from rad to deg
    data.Lat.values[:] = np.rad2deg(data.Lat.values)
    data.Lon.values[:] = np.rad2deg(data.Lon.values)

    data['Altitude'] = data['Height']
    data = data.drop('Height', axis=1)

    data.sort_index(inplace=True)

    return timeseries.TimeSeries(data, {'original header': header})


def read_csv(fname):
    """ reads in a piccolo log file or list of log files and returns a housekeeping instance
    """
    picco = None
    if type(fname).__name__ == 'list':
        first = True
        for file in fname:
            if '.log' not in file:
                print('%s is not a piccolo log file ... skipped' % file)
                continue
            print('%s ... processed' % file)
            picco_t = _read_file(file)
            if first:
                picco = picco_t
                first = False
            else:
                picco.data = pd.concat((picco.data, picco_t.data))

    else:
        picco = _read_file(fname)

    return picco
#
# class AutoPilot(object):
# def __init__(self, data, info):
# self.data = data
# self.info = info
