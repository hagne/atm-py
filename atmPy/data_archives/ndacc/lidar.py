"""
Module to read lidar data in the NDACC data base.

Download data from here:
    https://lidar.jpl.nasa.gov/ndacc/data/general.php

Not sure if there is an API
"""

from pyhdf.SD import SD, SDC
import pandas as pd
import xarray as xr
import numpy as np
# Replace 'your_file.hdf' with your actual HDF file path

def read_jpl_hdf(file_path):
    """
    Read files at leased for MLO jpl lidar files. There is a standard for the hdf lidar files, so this might work for hdf lidar files within NDACC.
    This is a very simple read only including very little metadata of what is included in the hdf files. Progamming required if you want more.
    """
    def extract_hdf_data(hdf_file):
        dataset_name = 'ALTITUDE'
        alt = hdf_file.select(dataset_name)
        
        dataset_name = 'DATETIME'
        dt = hdf_file.select(dataset_name)
        dtres = pd.to_datetime('2000-01-01') + pd.to_timedelta(dt[:], 'd')
        dt.endaccess()
        
        ds = xr.Dataset()
        var2get = ['AEROSOL.BACKSCATTER.RATIO_BACKSCATTER', 'AEROSOL.BACKSCATTER.COEFFICIENT_DERIVED']
        for dataset_name in var2get:
            dataset = hdf_file.select(dataset_name)
            ds[dataset_name.replace('.','_')] = xr.DataArray(dataset[:], coords = {'altitude': alt[:],})
            dataset.endaccess()
        
        ds =ds.expand_dims(datetime = dtres)
        
        alt.endaccess()
    
        return ds
        
    # Open the HDF file in read mode
    if isinstance(file_path, pl.Path):
        file_path = file_path.as_posix()
    hdf_file = SD(file_path, SDC.READ)
    ds = extract_hdf_data(hdf_file)
    hdf_file.end()
    return ds


def read_NOAA_ames(file_path):
    # Open and parse the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    columns = ['backscatter_ratio', 'backscatter_ratio_erro', 'density']
    # Initialize lists to store data
    altitude = []
    backscatter_ratio = []
    backscatter_ratio_error = []
    density_log = []
    
    # Parse the main data block
    lastline = ''
    lineiter = iter(lines)
    for line in lineiter:  # Adjust line number as needed
        if line.strip() == '0' and lastline.strip() == '0':
           break
        lastline = line
    
    # while 1:
    # header = next(lineiter).split()
    datablocks = []
    # thisistheend = False
    while 1:
        # print('next')
        try:
            header = next(lineiter).split()
        except StopIteration:
            # thisistheend = True
            break
        data = []
        for i in range(int(header[1])):
            line = next(lineiter)
            dline = line.split()
            dline = [int(i) for i in dline]
            data.append(dline)
    
    
        data = np.array(data)
        df = pd.DataFrame(data[:,1:], columns=columns, index = data[:,0])
        df.index.name = 'altitude'
    
        dst = df.to_xarray()
        # header = header.split()
        dt = pd.to_datetime(f'{header[2]}-{int(header[3]):02d}-{int(header[4]):02d} {int(header[5]):02d}:{int(header[6]):02d}:00')
        dst = dst.expand_dims(datetime = [dt])
        datablocks.append(dst)
    
    ds = xr.concat(datablocks, dim = 'datetime')
    return ds