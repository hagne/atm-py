#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:52:06 2025

@author: htelg
"""

import xarray as xr

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
    dss = dss.rename({'wavelength': 'wl_neph'})
    
    # particle_concentration
    
    dspc = xr.open_dataset(fn, group='data/particle_concentration')
    dspc = dspc.assign_coords({'time': ds.time})
    dspc = dspc.drop_vars(['sample_pressure','sample_temperature']) #They interfer with same variable in scattering, also they both were empty in my test file
    
    ## Absorbtion
    
    dsa = xr.open_dataset(fn, group='data/light_absorption')
    dsa = dsa.assign_coords({'time': ds.time})
    dsa = dsa.rename({'wavelength': 'wl_clap'})
    
    ## merge
    
    dsall = xr.merge([dss, dsa, dspc])
    
    return dsall