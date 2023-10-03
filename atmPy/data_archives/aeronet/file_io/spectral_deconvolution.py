#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 09:58:35 2022

@author: hagen
"""

import pandas as pd
import numpy as np
import atmPy.aerosols.size_distribution.diameter_binning as atmdb
import atmPy.aerosols.size_distribution.sizedistribution as atmsd
import atmPy.general.timeseries as atmts

def read_file(path2file, verbose = False):
    """
    So far it works only for version 2. There is a version 3 out now... some 
    programming will be needed

    Parameters
    ----------
    path2file : TYPE
        DESCRIPTION.

    Returns
    -------
    aero_inv : TYPE
        DESCRIPTION.

    """
    #### read the header
    header = ''
    with open(path2file) as rein:
        for i in range(3):
            header += rein.readline()
    
    # if 'Version 2' in header:
    #     version = 2
    #     skiprows = 3
    #     date_column = "Date(dd-mm-yyyy)"
    #     day_of_year = 'Julian_Day'
    #     time_column = "Time(hh:mm:ss)"
        
    if 'Version 3' in header:
        version = 3
        skiprows = 6
        day_of_year = 'Day_of_Year'
        if 'SDA Version 4.1' in header:
            date_column = 'Date_(dd:mm:yyyy)'
            time_column = "Time_(hh:mm:ss)"
            sda_version = '4.1'
        else:
            assert(False), f'unknown SDA version: {header}'
    else:
        assert(False), f'unkonwn aeronet version: {header}'
        
    if verbose:
        print(f'retrieval version: {version}')
        
    if version in [3,]:
        if verbose:
            print('Re-reading header')
        header = ''
        with open(path2file) as rein:
            for i in range(skiprows):
                header += rein.readline()
        
            
    #### read the data
    df = pd.read_csv(path2file, skiprows = skiprows)
            
    #### create timestamp
    df.index = df.apply(lambda row: pd.to_datetime(f'{row[date_column]} {row[time_column]}', format='%d:%m:%Y %H:%M:%S'), axis = 1)
    df = df.drop([date_column,time_column, day_of_year], axis = 1)
    df.index.name = 'datetime'
    
    #### parse the data and add to AeronetInversion instance
    ds = df.to_xarray()
    # replace invalid with nan
    ds = ds.where(ds != -999.0, np.nan)
    aero_inv = AeronetAODInversion(ds)
    aero_inv.header = header
    aero_inv.retrieval_version = version
    
    dist = extract_sizedistribution(df)
    # if dist:
    aero_inv.sizedistribution = dist
    
    
    return aero_inv

# def extract_singlescatteringalbedo(df, version):
#     """
#     Extract the single scattering albedo for all aerosols (there is no 
#     seperation into fine and coarse).

#     Parameters
#     ----------
#     df : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None.


#     """
#     if version == 2:
#         ssa_txt = 'SSA'
#     elif version == 3:
#         ssa_txt = 'Single_Scattering_Albedo'
    
#     # def 
#     ssa = df.loc[:,[i for i in df.columns if ssa_txt in i]]
#     ssa.columns = [''.join([e for e in i if e.isnumeric()]) for i in ssa.columns]
#     ssa.columns.name = 'channel (nm)'
    
#     return atmts.TimeSeries(ssa)
    
    
    
    
def extract_sizedistribution(df):   
    #### get the size distribution data
    cols = df.columns
    cols = [i for i in cols if i.replace('.','').isnumeric()]
    dist = df.loc[:, cols]
    if len(cols) == 0:
        return False
    
    # create bins for atmpy
    bins, _ = atmdb.bincenters2binsANDnames(np.array([float(i) for i in cols]))
    bins*=2 #radius to diameter
    bins*=1e3 # um to nm
    
    #### create sizedistribution instance
    #### todo: there is a scaling error since AERONET uses 'dVdlnDp' and I use 'dVdlogDp'
    dist_ts  = atmsd.SizeDist_TS(dist, bins, 'dVdlogDp', 
                                  # fill_data_gaps_with=np.nan, 
                                  ignore_data_gap_error=True,
                                  ) 
    
    return dist_ts

class AeronetAODInversion(object):
    def __init__(self, data):
        self.data = data