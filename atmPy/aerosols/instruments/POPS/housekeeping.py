# -*- coding: utf-8 -*-
"""
@author: Hagen Telg
"""
import datetime

import pandas as pd
# import os
# import pylab as plt
# from atmPy.tools import conversion_tools as ct
from atmPy.general import atmosphere_standards as atm_std, timeseries


def _read_housekeeping(fname):
    """Reads housekeeping file (test_data_folder; csv-format) returns a pandas data frame instance."""
    print('reading %s'%fname)
    try:
        df = pd.read_csv(fname, error_bad_lines=False)
    except ValueError:
        return False
#    data = df.values
#    dateString = test_data_folder.split('_')[0]
    dt = datetime.datetime.strptime('19700101', "%Y%m%d") - datetime.datetime.strptime('19040101', "%Y%m%d") 
    dts = dt.total_seconds()
    # todo: (low) what is that delta t for, looks fishi (Hagen)
    dtsPlus = datetime.timedelta(hours=0).total_seconds()
    # Time_s = data[:,0]
    # data = data[:,1:]
    df.index = pd.Series(pd.to_datetime(df.Time_s-dts-dtsPlus, unit = 's'), name = 'Time_UTC')
    # if 'P_Baro' in df.keys():
    #     df['barometric_pressure'] = df.P_Baro
    #     df.drop('P_Baro', 1, inplace=True)
    #     df['altitude'] = ct.p2h(df.barometric_pressure)
    return POPSHouseKeeping(df)


def read_csv(fname, verbose = False):
    """
    Parameters
    ----------
    fname: string or list of strings.

    Returns
    -------
    TimeSeries instance
    """
    # test_data_folder = os.listdir()
    # test_data_folder = '20150419_000_POPS_HK.csv'

    houseKeeping_file_endings = ['HK.csv', 'HK.txt']

    first = True

    if type(fname).__name__ == 'list':
        for file in fname:
            for i in houseKeeping_file_endings:
                if i in file:
                    is_hk = True
                    break
                else:
                    is_hk = False
            if verbose and not is_hk:
                print('%s is not a housekeeping file ... continue')

            if is_hk:
                hktmp = _read_housekeeping(file)
                if not hktmp:
                    print('%s is empty ... next one' % file)
                elif first:
                    data = hktmp.data.copy()
                    first = False
                    hk = POPSHouseKeeping(data)
                    # continue
                else:
                    data = pd.concat((data, hktmp.data))
                    hk = POPSHouseKeeping(data)
        if first:
            txt = """Either the prvided list of names is empty, the files are empty, or none of the file names end on
the required ending (*HK.csv)"""
            raise ValueError(txt)
    else:
        hk = _read_housekeeping(fname)
    hk.data = hk.data.dropna(how='all')  # this is necessary to avoid errors in further processing

    if 'P_Baro' in hk.data.keys():
        hk.data['Barometric_pressure'] = hk.data.P_Baro
        hk.data.drop('P_Baro', 1, inplace=True)
        # hk.data['Altitude'] = ct.p2h(hk.data.barometric_pressure)

    return hk


class POPSHouseKeeping(timeseries.TimeSeries):
    def get_altitude(self, temperature=False):
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
        alt, tmp = atm_std.standard_atmosphere(self.data.Barometric_pressure, quantity='pressure')
        self.data['Altitude'] = alt
        return alt

# todo: (low) this has never been actually implemented
# def read_housekeeping_allInFolder(concatWithOther = False, other = False, skip=[]):
# """Read all housekeeping files in current folder and concatinates them.
#     Output: pandas dataFrame instance
#     Parameters
#         concatWithOther: bool, if you want to concat the created data with an older set given by other
#         other: dataframe to concat the generated one
#         skip: list of file which you want to exclude"""
#
#     files = os.listdir('./')
#     if concatWithOther:
#         counter = True
#         hkdf = other.copy()
#     else:
#         counter = False
#     for e,i in enumerate(files):
#         if 'HK.csv' in i:
#             if i in skip:
#                 continue
#             hkdf_tmp = read_housekeeping(i)
#             if not counter:
#                 hkdf = hkdf_tmp
#             else:
#                 hkdf = pd.concat([hkdf,hkdf_tmp])
#             counter = True
#     return hkdf

