#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 10:50:28 2021

@author: hagen
"""
import numpy as _np
import pandas as _pd

import atmPy.aerosols.size_distribution.sizedistribution as _sd
import atmPy.aerosols.size_distribution.diameter_binning as _db

def uhsas2sizedist(df):
    """
    Creates size distribution time series instance from uhsas data (as 
    returned by the read_file function)

    Parameters
    ----------
    df : pandas.DataFrame
        as put out by the read_file function.

    Returns
    -------
    dist : TYPE
        DESCRIPTION.

    """
    ## make bins (based on whats mentioned in the header)
    bins = _np.linspace(40,1000, 99)

    ## the size distribution data
    data = df.iloc[:,:-1].copy()
    data.columns = bins

    ### to my knowledge the uhsas can not measure below ~70 nm
    data_trunc = data.loc[:,69:]

    # make the size distribution
    bined, _ = _db.bincenters2binsANDnames(data_trunc.columns.values)
    dist = _sd.SizeDist_TS(data_trunc, bined, 'numberConcentration')
    return dist

def read_file(fn, guess_product = True):
    """
    Reads a file with the ICARTT format.
    https://www-air.larc.nasa.gov/missions/intexna/DataManagement_plan.htm

    Parameters
    ----------
    fn : TYPE
        DESCRIPTION.
    guess_product : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    """
    # get header
    with open(fn, 'r') as rein:
        fl = rein.readline()
        header_length = int(fl.split(',')[0])
        header_lines = [fl,]
        for i in range(header_length-1):
            header_lines.append(rein.readline())
    
    # get data
    ## get date
    date_sp = [i.strip() for i in header_lines[6].split(',')]
    date = _pd.to_datetime(f'{date_sp[0]}-{date_sp[1]}-{date_sp[2]}')
    
    names = [i.strip() for i in header_lines[-1].split(',')]
    df = _pd.read_csv(fn, skiprows = header_length, names=names)
    df.index = df.apply(lambda row: date + _pd.to_timedelta(row.name, unit = 's'), axis = 1)
    df.index.name = 'datetime'
    df[df == -999999.00000] = _np.nan
    df = df.drop('UTC_datetime', axis = 1)
    
    if header_lines[3].strip().lower() == 'uhsas':
        return uhsas2sizedist(df)
    
    return header_lines,df