#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 14:00:39 2022

@author: hagen
"""
import xarray as xr
import pandas as pd
import numpy as np

def read_mfrsr_cal(path2file):
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
    df.drop('Open', inplace = True)
    df.index = df.index.astype(int)
    df.columns.name = 'stats'
    ds['statistics'] = df
    
    
    #### read the data
    df = pd.read_csv(path2file, 
                skiprows=skiprows,
                # encoding='unicode_escape', 
                encoding = 'ISO-8859-1',
                sep = '\t')
    
    df = df[1:].copy()
    
    df.index = df.WAVEL
    
    df.drop('WAVEL', axis=1, inplace=True)
    
    
    df.index.name = 'wavelength'
    df.index = df.index.astype(np.float32)
    
    # open responds, not sure what exactly this is?
    ds['responds_open'] = df['0']
    ds['responds_open_error'] = df['0']
    df.drop(['0','0ERR'], axis = 1, inplace = True)
    
    # responds
    df_res = df.loc[:,[c for c in df.columns if not 'ERR' in c]]
    df_res.columns = [int(c) for c in df_res.columns]
    df_res.columns.name = 'channel'
    # when a scan section ends with a non-zero value convolution with e.g. a RTM irradiance spectum will lead to problems! Therefore -> 
    df_res[df_res == 0] = np.nan # this can potenoally lead to nans where we don't want them .... consider interpolating afterwards ... as there is no extrapolating this should not lead to a reemergance of the above problem.
    ds['responds'] = df_res
    
    # respnds error
    df_err = df.loc[:,[c for c in df.columns if 'ERR' in c]]
    df_err.columns = [int(c.replace('ERR','')) for c in df_err.columns]
    df_err.columns.name = 'channel'
    ds['responds_error'] = df_res
    
    return ds