#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:52:06 2025

@author: htelg
"""

import xarray as xr
import pandas as pd
import numpy as np

def read_T640(p2f):
    """
    Read T640 data as generate within NFAN

    Parameters
    ----------
    p2f : TYPE
        typical locations of the data are '/nfs/iftp/aftp/aerosol/bos/T640'.

    Returns
    -------
    ds : xarray.Dataset
        DESCRIPTION.

    """
    df = pd.read_csv(p2f, 
                     # skiprows=1,
                     header = [0,1]
                    )
    # df
    
    
    
    df.index = pd.to_datetime(df.iloc[:,0])
    
    df.index.name = 'datetime'
    
    df = df.drop('Date String (YYYY-MM-DD hh:mm:ss) UTC', axis=1)
    
    
    var_dict = dict(F1 = 'instrument_flag',
                    P = 'pressure',
                    Q1 = 'flow_sample',
                    Q2 = 'flow_bypass',
                    T1 = 'temperature_sample',
                    T2 = 'temperature_ambient',
                    T3 = 'temperature_conditioner',
                    T4 = 'temperature_led',
                    T5 = 'temperature_box',
                    U1 = 'rh',
                    ZSPAN = 'span_deviation',
                    X0 = 'pm10',
                    X1 = 'pm1',
                    X2 = 'pm2_5',
                   )
    
    ds = xr.Dataset()
    
    for lab, col in df.transpose().iterrows():
        lab_abb = lab[1].split('_')[0]
        lab_ds = var_dict[lab_abb]
        if lab_abb == 'F1':
            ds[lab_ds] = col
        else:
            ds[lab_ds] = col.astype(np.float32)
        ds[lab_ds].attrs['long_name'] = lab[0]

    ds.attrs = dict(instrument = 'T640')
    return ds

def open_nfan(fn = '/nfs/iftp/aftp/aerosol/bos/2020/ESRL-GMD-AEROSOL_v1.0_HOUR_BOS_s20200101T000000Z_e20210101T000000Z_c20210214T053935Z.nc'):
    '''
    Open a standard nfan file. You might have to install h5py for this to work?

    Parameters
    ----------
    fn : TYPE, optional
        path to file.

    Returns
    -------
    dsall : TYPE
        DESCRIPTION.

    '''
    
    ## time, pressure and temp
    
    ds = xr.open_dataset(fn, group='data')
    
    ## light_scattering
    
    dss = xr.open_dataset(fn, group='data/light_scattering')
    dss = dss.assign_coords({'time': ds.time})#xr.open_dataset(fn, group='data').time})
    if 'wavelength' in dss.variables:
        dss = dss.rename({'wavelength': 'wl_neph'})
    
    # particle_concentration
    
    dspc = xr.open_dataset(fn, group='data/particle_concentration')
    dspc = dspc.assign_coords({'time': ds.time})
    
    if 'sample_pressure' in dspc.variables:
        dspc = dspc.drop_vars(['sample_pressure']) #They interfer with same variable in scattering, also they both were empty in my test file
    if 'sample_temperature' in dspc.variables:
        dspc = dspc.drop_vars(['sample_temperature']) #They interfer with same variable in scattering, also they both were empty in my test file
    ## Absorbtion
    
    dsa = xr.open_dataset(fn, group='data/light_absorption')
    dsa = dsa.assign_coords({'time': ds.time})
    if 'wavelength' in dsa.variables:
        dsa = dsa.rename({'wavelength': 'wl_clap'})
    
    ## merge
    
    dsall = xr.merge([dss, dsa, dspc])
    
    return dsall