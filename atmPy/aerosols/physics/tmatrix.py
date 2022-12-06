#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 21:06:16 2022

@author: hagen
"""

import pytmatrix.tmatrix
from pytmatrix import tmatrix
import pytmatrix.scatter
import multiprocessing as mp
import numpy as np
import xarray as xr


def runthis(return_dict, d, ndgs, ddelt):
    # print('.', end = '', flush=True)
    return_dict[d] = np.nan
    do = d
    scale = 1e-12
    wl = 500 * 1e-3 * scale # the  * 1e-3 is to get to the right unit, the scale to improve the speed, no idea why?!?
    # d = .500 #
    ratio = 1#3.
    # d = 2822* scale #limit with ndgs = 5
    # d = 500 * scale
    d = d * 1e-3 * scale
    scatterer = tmatrix.Scatterer(radius=d/2, 
                                  wavelength=wl, 
                                  m=complex(1.5, 0), 
                                  axis_ratio=ratio,
                                  ndgs = ndgs,
                                  ddelt = ddelt, 
                                  # orient = pytmatrix.orientation.orient_single,

                                 )
    scat = pytmatrix.scatter.sca_xsect(scatterer, 
                                       # h_pol=True,
                                       )
    res = scat / scale**2
    # results.append(res)
    # results['test'] = res
    # print(res)
    return_dict[do] = res
    print('|', end = '', flush=True)
    return res

def diameter2ds(diameters):
    data = np.zeros(diameters.shape)
    data[:] = np.nan
    
    ds = xr.Dataset()
    ds['scatt_cross_scect'] = xr.DataArray(data.copy(), coords={'diameter': diameters})
    ds['ddelt'] = xr.DataArray(data.copy(), coords={'diameter': diameters})
    ds['ndgs'] = xr.DataArray(data.copy(), coords={'diameter': diameters})
    return ds

# if __name__ == "__main__":
    
def playground(ds = None):
    # if 0:
    #     runthis()
    # elif 1:
    try:
        mp.set_start_method('spawn')
    except RuntimeError as er:
        if er.args[0] !='context has already been set':
            raise
            
    manager = mp.Manager()
    return_dict = manager.dict()
    
    if isinstance(ds, type(None)):
        diameters = np.array([500,
                  # 2822, 
                  # 500,
                  ])
        # diameters = np.logspace(np.log10(3050), np.log10(15000), 100)  # converged for all numbers, many at  
        ds = diameter2ds(diameters)

        
        
    ndgsmax = 14
    ndgs_list = np.linspace(2, ndgsmax, ndgsmax-1) 
    ddeltlist =  [1e-3, 1e-2, 1e-1, 1, 1e1]
    
    # just for testing
    # test = False
    # if test:
    #     diameters = [3050]
    #     ndgsmax = 20
    #     ndgs_list = np.linspace(2, ndgsmax, ndgsmax-1) 
    #     ddeltlist = [0.1]
        
    # for d in diameters:
    for d in ds.diameter:
        d = float(d)
        # print(f'{d}', end = '\t')
        dsel = ds.sel(diameter = [d])
        if not np.isnan(float(dsel.ddelt[0])):
            continue
        # ddeltlist = [1e-3,1e-4,1e-5,]
        for ddelt in ddeltlist:
            # res_list = []
            for ndgs in ndgs_list:
                # ndgs = i + 1
                # print(ndgs, end = ' ')
                # print('>', end = '')
                process = mp.Process(target = runthis, 
                                      args = (return_dict,
                                              d,
                                              ndgs,
                                              ddelt,
                                              ),
                                    )
                process.start()
                process.join()
                

                if np.isnan(return_dict[d]):
                    continue
                else:
                    break    
            if not np.isnan(return_dict[d]):
                break
                            
        # print(f'{ddelt}', end = ' ')                    
        # print(f'{ndgs}', end = ' ')
        ds.ddelt.loc[dict(diameter = d)] = np.log10(ddelt)
        ds.ndgs.loc[dict(diameter = d)] = ndgs    
        
        if not np.isnan(return_dict[d]):
            # print(f'{return_dict[d]}',
            #       # end = '',
            #       )
            ds.scatt_cross_scect.loc[dict(diameter = d)] = return_dict[d]
        else:
            # print('Fail!!',
            #       # end = '',
            #       )
        
            
            ds.ddelt.loc[dict(diameter = d)] = np.log10(ddelt)
            ds.ndgs.loc[dict(diameter = d)] = ndgs
    # if test:
    #     mean = np.mean(res_list)
    #     std = np.std(res_list)
    #     print(f'mean: {mean}')
    #     print(f'std: {std}')
    #     print(f'rel. std {(1 - (mean - std)/mean) *100}')
        
    # print('')
    # print(return_dict)
    return ds
        # print(process.exitcode)
    
    # else:
    #     process.run()
    
    # results
            
        