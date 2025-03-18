#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 14:00:39 2022

@author: hagen
"""
import xarray as xr
import pandas as pd
import numpy as np
import io

def read_mfrsr_cal(path2file, verbose = False):
    """
    This reads the Spectral Responds tests, which are files that usually end on .SPR, e.g. V0648_09721Cd.SPR

    Parameters
    ----------
    path2file : TYPE
        DESCRIPTION.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    """
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

def read_mfrsr_cal_cos(p2f):
    """
    Reads cosine corrections commonly used in GRAD.

    Parameters
    ----------
    p2f : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    with open(p2f,'r') as rein:
        # pass
        lines = rein.readlines()
    
    verbose = True
    header = []
    header_done = False
    data1 = []
    data1_done = False
    data2 = []
    for e,l in enumerate(lines):
        if l.strip() == 'Tabulated_Response_Data_Lables:':
            if verbose:
                print('Done reading header')
            header_done = True
            continue
        elif l.strip() == 'Labels:':
            if verbose:
                print('Done reading data1')
            data1_done = True
            continue
    
        lc = l.strip()
        if len(lc) == 0:
            if verbose:
                print(f'Line {e} is empty')
            continue
        if not header_done:
            header.append(lc)
        elif not data1_done:
            data1.append(lc)
        else:
            data2.append(lc)
    
    def parse_data(data, verbose = True):
        data.pop(1) #remove the units
        data.pop(1) #remove some text
    
        data = '\n'.join(data)
        df = pd.read_csv(io.StringIO(data),sep = '\t' )
        df.index = df.Angle
        # test if NS (north-south) or EW (...)
        if 'NS' in df.columns[2]:
            scan_direction = 'NS'
        elif 'EW' in df.columns[2]:
            scan_direction = 'EW'
        else:
            ValueError('Charlse probably parsed this file differently ... get in there and fix it')
        
        
    
        # get & format only relevant colums, consider getting the other ones at some point too?
        cols = [col for col in df.columns if f'{scan_direction}_Corr' in col]
        df_cc = df.loc[:, cols]
        df_cc = df_cc.dropna(axis = 1) #there are placeholder collums for the UV MFRSRs, this gets rid of them ... or the other way aroudn
        cols = [col.strip(f'{scan_direction}_Corr_') for col in df_cc.columns]
        cols = [col.strip('nm') for col in cols]
        df_cc.columns = cols
    
        # 90 and -90 degrees are singularities ==> remove
        for i in [-90, 90]:
            if (idx:= i) in df.index:
                df_cc.drop(idx, inplace=True)
    
        # other stuff
        df_cc.columns.name = 'channel'
        return dict(data = df_cc, scan_direction = scan_direction)
    
    data1dict = parse_data(data1)
    data2dict = parse_data(data2)
    
    ds = xr.Dataset()
    
    scand = data1dict['scan_direction']
    ds[f'broadband_{scand}'] = data1dict['data'].Thermopile
    ds[f'spectral_{scand}'] = data1dict['data'].drop('Thermopile', axis = 1)
    
    scand = data2dict['scan_direction']
    ds[f'broadband_{scand}'] = data2dict['data'].Thermopile
    ds[f'spectral_{scand}'] = data2dict['data'].drop('Thermopile', axis = 1)
    
    ds.attrs['header'] = '\n'.join(header)
    return ds

def read_mfrsr_cal_responsivity(p2f):
    xls = pd.ExcelFile(p2f)
    
    df = xls.parse(sheet_name='D.ABS', 
                   # header = 37,
                  )
    dfhead = df.iloc[:32, :2]
    df = df.iloc[39:, :2]
    df.columns = ['channel', 'responsivity']
    df.index = df.channel
    df.drop('channel', inplace = True, axis = 1)
    assert(df.index[0] == 'Thermopile'), 'Format is inconsistant!!! Kick Charles in the ...!'
    
    ds = xr.Dataset(df)
    
    header = {row[dfhead.columns[0]]: row[dfhead.columns[1]] for _,row in dfhead.iterrows()}
    header['Serial_number'] = header.pop('Instrument')
    ds.attrs = header
    
    # lamp calibration
    df = xls.parse(sheet_name='licor_1030L_cal_data', header = 7)
    assert(df.columns[0] == 'nm'), 'arrrg, files are inconsistant'
    df_cal_lamp = df.drop('mw/m^2-nm', axis = 1)
    df_cal_lamp.rename({'nm': 'wavelength', 'W/m^2-nm': 'intensity'}, axis=1, inplace=True)
    df_cal_lamp.index = df_cal_lamp.wavelength
    df_cal_lamp.drop('wavelength', axis = 1, inplace=True)
    
    ds['lamp_data'] = df_cal_lamp.intensity
    
    # read the spectral intensity
    # just in case I want to put it back in at some point. This is the same data as in the SPR files
    if 0:
        df = xls.parse(sheet_name='SPR_Data', 
                       header = 47, 
                      )
        df = df.iloc[1:] #in my case there was a line with useless text
        
        assert(df.iloc[0,0] == 395), 'This should be the topleft corner of the data section, showing the lowest wavelength. Unless does not keep this the same this is expected to be 395)'
        
        df.columns = df.columns.astype(str)
        
        # for some reason charles often does not name the first collumn correctly, arrrrrrg
        if df.columns[1] == '940':
            df.rename({'940': 'broadband'}, axis = 1, inplace=True)
            df.rename({'940.1': '940'}, axis = 1, inplace=True)
        
        df.rename({'WAVEL': 'wavelength'}, axis = 1, inplace = True)
        
        df.drop([c for c in df.columns if 'ERR' in c], axis = 1, inplace = True)
        
        df.index = df.wavelength
        
        df_spc = df.drop('wavelength', axis = 1)

    # change the channel wavelength to nominal wavelent and add the channel wavelenth as a variable
    wl_nominal = [415, 500, 1625,670,870,940]
    wl_nominal = np.array(wl_nominal)
    ch_nominal = []
    channels_orig = ds.channel.values
    for c in channels_orig:
        if isinstance(c, str):
            ch_nominal.append(c)
        else:
            nc = wl_nominal[abs(c - wl_nominal).argmin()]
            ch_nominal.append(nc)
    
    ds = ds.assign_coords(channel = ch_nominal)
    
    ds['channel_wavelength'] = xr.DataArray(channels_orig, dims=('channel',))

    return ds