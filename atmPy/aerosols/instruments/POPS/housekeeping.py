# -*- coding: utf-8 -*-
"""
@author: Hagen Telg
"""
import datetime

import pandas as pd
import numpy as _np
import os
# import pylab as plt
# from atmPy.tools import conversion_tools as ct
from atmPy.general import timeseries
from atmPy.atmosphere import standards as atm_std
import atmPy.aerosols.size_distribution.sizedistribution as atmsd
import pathlib as pl



def read_file(path,
              version = 'BBB_02',
              pattern = 'HK',
              skip_histogram = False,
              size_bins = None,
              # calibration_file = None,
              ignore_colums = [],  #['Flow_Rate_ccps', 'LED_P_MON', 'AI_4', 'AI_5', 'AI_7', 'AI_8', 'AI_9', 'AI_10', 'AI_11', 'LED_P_Mon_Therm', 'AO_Flow', 'AO_LaserPower', 'No_Pts', 'ValidParts', 'writeTime', 'currMax'],
              verbose = False,
              ):
    """
    Parameters
    ----------
    path: string or list of strings.
        This can either be a file name, a list of filenames or a folder.
    pattern: str
        if folder is given than this is the pattern housekeeping files will be identified by.
    version: string ['BBB_01']
        BBB_02: Hendix version, not sure since when. At least since 2022-06, but 
                most likely way earlier ...
        BBB_01: Beagle bone (original)
        sbRio: sbRio
    size_bins: int or pathlib.Path
        Path to a file containing the bin edges (EDGES not CENTERS!!). Structure: currently one value per line.
    verbose: bool
    Returns
    -------
    TimeSeries instance
    """
    # test_data_folder = os.listdir()
    # test_data_folder = '20150419_000_POPS_HK.csv'


    def read_sbRio(fname, skip_histogram = False, verbose=False):
        """Reads housekeeping file (test_data_folder; csv-format) returns a pandas data frame instance.
        """
        if verbose:
            print('reading %s' % fname)
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
        df.index = pd.Series(pd.to_datetime(df.Time_s - dts - dtsPlus, unit='s'), name='Time_UTC')
        # if 'P_Baro' in df.keys():
        #     df['barometric_pressure'] = df.P_Baro
        #     df.drop('P_Baro', 1, inplace=True)
        #     df['altitude'] = ct.p2h(df.barometric_pressure)
        return POPSHouseKeeping(df)

    def read_BBB(fname, skip_histogram = False, verbose = False):
        if verbose:
            print(f'read pops house keeping bbb file: {fname}')
        col_names = pd.read_csv(fname, sep=',', nrows=1, header=None,
                                #             index_col=1,
                                #             usecols=np.arange()
                                ).values[0][:-1].astype(str)
        col_names = _np.char.strip(col_names)

        if skip_histogram:
            usecols = list(range(27))
        else:
            usecols = None
        data = pd.read_csv(fname, sep=',', skiprows=1, header=None, usecols = usecols
                           #             index_col=1,
                           #             usecols=np.arange()
                           )

        data_hk = data.iloc[:, :27]
        data_hk.columns = col_names
        data_hk.index = pd.to_datetime(data_hk['DateTime'], unit='s')
        data_hk.drop('DateTime', axis=1, inplace=True)
        #     hk = atmPy.general.timeseries.TimeSeries(data_hk, sampling_period = 1)

        hk = POPSHouseKeeping(data_hk, sampling_period=1)
        hk.data['Barometric_pressure'] = hk.data['P']
        return hk
    
    def read_BBB_02(fname, skip_histogram = False, verbose = False):
        if verbose:
            print(f'read pops house keeping file: {fname}')

        if skip_histogram:
            usecols = list(range(27))
        else:
            usecols = None
        data = pd.read_csv(fname, sep=',', 
                           usecols = usecols
                           )

        data.columns = [col.strip() for col in data.columns]
        data.index = pd.to_datetime(data['DateTime'], unit='s')
        data.drop('DateTime', axis=1, inplace=True)

        hk = POPSHouseKeeping(data, sampling_period=1)
        hk.data['Barometric_pressure'] = hk.data['P']
        return hk
    
    dist = f'Extraction of the sizedistribution is currently not implemented for the file_version {version}'
    
    #### assign version
    if version == 'sbRio':
        read = read_sbRio
    elif version == 'BBB_01':
        read = read_BBB
    elif version == 'BBB_02':
        read = read_BBB_02
    else:
        raise ValueError('Housekeeping version {} is unknown!'.format(version))

    #### workplan
    path = pl.Path(path)
    if path.is_dir():
        file_paths = sorted(list(path.glob('*{}*'.format(pattern))))
    elif path.is_file():
        file_paths = [path]
    elif type(path) == list:
        file_paths = path
    else:
        raise TypeError('fname is of unknown type: {}'.format(type(path).__name__))

    file_paths.sort()
    
    #### read files
    hk_data = []
    for file in file_paths:

        hktmp = read(file, skip_histogram=skip_histogram, verbose=verbose)
        hk_data.append(hktmp.data)
        
    data = pd.concat(hk_data)
    

        
            
    #### generate POPSHouseKeeping instance and condition data
    hk = POPSHouseKeeping(data)
    hk.data = hk.data.dropna(how='all')  # this is necessary to avoid errors in further processing

    if ('P_Baro' in hk.data.keys()) or ('P_Ambient' in hk.data.keys()):
        if 'P_Baro' in hk.data.keys():
            hk.data['Barometric_pressure'] = hk.data.P_Baro
            hk.data.drop('P_Baro', 1, inplace=True)
        if 'P_Ambient' in hk.data.keys():
            hk.data['Barometric_pressure'] = hk.data.P_Ambient
            hk.data.drop('P_Ambient', 1, inplace=True)
            # try:
                # hk.data['Altitude'] = ct.p2h(hk.data.barometric_pressure)

    if ignore_colums:
        hk.data = hk.data.drop(ignore_colums, axis=1)
    
    #### separate housekeeping and sizedistribution
    if version == 'BBB_02':
        data = hk.data
        hist_cols  = [col for col in data.columns if (col[0] == 'b' and col[1:].isnumeric())]
        dist = data.loc[:,hist_cols]
        data.drop(hist_cols, axis=1, inplace = True)
        
        
        #### read size bin file
        fn = pl.Path(size_bins)
        with open(fn, 'r') as rein:
            lines = rein.readlines()
        bins = _np.array([float(l) for l in lines])
        
        #### generate size distribution timeseries instance
        dist = atmsd.SizeDist_TS(dist, bins, 'numberConcentration')
        
        dist.housekeeping = hk
    return {'housekeeping': hk, 'sizedistribution': dist}


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
        alt, tmp = atm_std.standard_atmosphere(self.data.loc[:,'Barometric_pressure'].copy(), quantity='pressure')
        self.data['Altitude'] = alt
        return alt


