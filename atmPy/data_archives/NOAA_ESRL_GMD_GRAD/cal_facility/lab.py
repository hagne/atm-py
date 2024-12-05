#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 14:00:39 2022

@author: hagen
"""
import xarray as xr
import pandas as pd
import numpy as np

def read_mfrsr_cal(path2file, verbose = False):
    ds = xr.Dataset()
    
    #### read the header
    with open(path2file, 'r', encoding='ISO-8859-1') as rein:
        maxit = 50
        header = []
        for i in range(maxit):
            # print(i)
            line = rein.readline()
            # try:
            #     line = rein.readline()
            # except UnicodeDecodeError:
            #     continue
            if 'Caculated_results_lables' in line:
                skiprowsstats = i + 1
                header_txt = '\n'.join([h.strip() for h in header])
            elif "Tabulated_Response_Data_Lables" in line:
                break
            else:
                header.append(line)
        assert(i != range(maxit)[-1]), 'End of header not reached, increase maxit?'

    skiprows = i+1
    
    
    
    ds.attrs['header'] = header_txt
    
    #### read the channel statistics
    df = pd.read_csv(path2file, 
                skiprows=skiprowsstats,
                nrows = skiprows-skiprowsstats-3,
                # encoding='unicode_escape', 
                encoding = 'ISO-8859-1',
                sep = '\t')
    
    df.index = df.Band
    df.index.name = 'channel'
    df.drop('Band', axis = 1, inplace=True)
    df = df.loc[:,[c for c in df.columns if not 'Unnamed' in c]]
    if 'Open' in df.index:
        df.drop('Open', inplace = True)
    elif 'Thermopile' in df.index:
        df.drop('Thermopile', inplace = True)
    elif 'Thermop' in df.index: # for some reason the name of this row varies
        df.drop('Thermop', inplace = True)
    else:
        assert(False), 'what else?'
     
    #### rename stats channels to nominal channels
    wl_nominal = np.array([415,500,1625, 615, 670, 870, 940])
    df.rename({col: wl_nominal[abs(wl_nominal - int(col)).argmin()] for col in df.index}, axis = 0, inplace = True)   
    # df.index = df.index.astype(int)
    df.columns.name = 'stats'
    ds['statistics'] = df
    
    #### read the data
    if verbose:
        print(f'skiprows before data: {skiprows}')
    df = pd.read_csv(path2file, 
                skiprows=skiprows,
                # encoding='unicode_escape', 
                encoding = 'ISO-8859-1',
                # sep = '\t',
                    sep=r"\s+")
    # return None, df
    if verbose:
        print('responds data ... inital reading')
        print(df)
    df = df[1:].copy()
    
    df.index = df.WAVEL
    
    df.drop('WAVEL', axis=1, inplace=True)
    
    
    df.index.name = 'wavelength'
    df.index = df.index.astype(np.float32)
    # return df
    # return None, df
    # if 'responds_open' in df.columns:
    #     # open responds, not sure what exactly this is?
    #     ds['responds_open'] = df['0']
    #     ds['responds_open_error'] = df['0ERR']
    # elif 'Thermop' in df.columns:
    #     ds['Thermop']
    
    # else:
    #     assert(False), 'what else?'

    # These are the thermopile columns which are always have a different name!
    df.drop([df.columns[0], df.columns[7]], axis = 1, inplace = True)
    df.columns = [c.replace('.1', '') for c in df.columns]
    
    if 0:
        for c in df.columns:
            if c in ['0', '0ERR', 'Therm', 'ThermERR', 'Thermop', 'ThermopERR', '2000', '2000ERR', 'Thermopile', 'ThermopileERR']:
                df.drop([c], axis = 1, inplace = True)
    # if '0' in df.columns:
    #     df.drop(['0','0ERR'], axis = 1, inplace = True)
    # elif 'Thermop' in df.columns:
    #     df.drop(['Thermop','ThermopERR'], axis = 1, inplace = True)
    # elif '2000' in df.columns:
    #     df.drop(['2000','2000ERR'], axis = 1, inplace = True)  
    
    df.drop([col for col in df.columns if 'Unnamed' in col], axis = 1, inplace = True)
    
    # return ds
    #### rename the channel wavelength to the nominal wavelength by finding the closest
    df.rename({col: f"{wl_nominal[abs(wl_nominal - int(col.strip('ERR'))).argmin()]}ERR" for col in df.columns if 'ERR' in col}, axis = 1, inplace = True)
    df.rename({col: str(wl_nominal[abs(wl_nominal - int(col)).argmin()]) for col in df.columns if 'ERR' not in col}, axis = 1, inplace = True)
    # return df
    #### responds
    df_res = df.loc[:,[c for c in df.columns if not 'ERR' in c]]
    df_res.columns = [int(c) for c in df_res.columns]
    df_res.columns.name = 'channel'
    # when a scan section ends with a non-zero value convolution with e.g. a RTM irradiance spectum will lead to problems! Therefore -> 
    df_res[df_res == 0] = np.nan # this can potenoally lead to nans where we don't want them .... consider interpolating afterwards ... as there is no extrapolating this should not lead to a reemergance of the above problem.
    # return ds, df_res
    if verbose:
        print('df_res')
        print(df_res)
    ds['responds'] = df_res
    
    #### respnds error
    # return ds
    df_err = df.loc[:,[c for c in df.columns if 'ERR' in c]]
    df_err.columns = [int(c.replace('ERR','')) for c in df_err.columns]
    df_err.columns.name = 'channel'
    ds['responds_error'] = df_res
    
    return ds