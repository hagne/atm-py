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


def round_sig(x, sig=2):
    out = np.zeros(x.shape)
    for e,xi in enumerate(x):
        out[e] = round(xi, sig-(np.floor(np.log10(abs(xi)))).astype(int)-1)
    return out

def get_samplingpoints(dmin, dmax, scale = 'log', round2 = 0, iteration = 0):
    bla = lambda n: n-1 + n
    reslist = [3]
    for i in range(iteration):
        reslist.append(bla(reslist[-1]))
    if scale == 'log':       
        d = np.logspace(np.log10(dmin), np.log10(dmax), reslist[-1])
        d = round_sig(d, round2)
    elif scale == 'lin':
        d = np.round(np.linspace(dmin, dmax, reslist[-1]), round2)
    else:
        assert(False)
        
    d = np.unique(d)
    return d

def calculate_optical_properties(return_dict, d,
                                    pshape,
                                    n_real,
                                    n_imag, 
                                    ndgs, ddelt, verbose = False):
    """perform t-matrix approximation for particular set of parameters"""
    
    if verbose:
        print('.', end = '', flush=True)
    return_dict[d] = np.nan
    do = d
    scale = 1e-12
    wl = 500 * 1e-3 * scale # the  * 1e-3 is to get to the right unit, the scale to improve the speed, no idea why?!?
    # d = .500 #
    ratio = pshape
    # d = 2822* scale #limit with ndgs = 5
    # d = 500 * scale
    d = d * 1e-3 * scale
    scatterer = tmatrix.Scatterer(radius=d/2, 
                                  wavelength=wl, 
                                  m=complex(n_real, n_imag), 
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
    if verbose:
        print('|', end = '', flush=True)
    return res

def deprecated_generate_lut(diameters):
    data = np.zeros(diameters.shape)
    data[:] = np.nan
    
    ds = xr.Dataset()
    ds['scatt_cross_scect'] = xr.DataArray(data.copy(), coords={'diameter': diameters})
    ds['status'] = xr.DataArray(data.copy(), coords={'diameter': diameters})
    # ds['ndgs'] = xr.DataArray(data.copy(), coords={'diameter': diameters})
    return ds

# if __name__ == "__main__":
    
def deprecated_work_lut(ds = None):
    

    try:
        mp.set_start_method('spawn')
    except RuntimeError as er:
        if er.args[0] !='context has already been set':
            raise
            
    manager = mp.Manager()
    return_dict = manager.dict()
    
    if isinstance(ds, type(None)):
        diameters = np.array([500, ])
        ds = generate_lut(diameters)

    ndgsmax = 14
    ndgs_list = np.linspace(2, ndgsmax, ndgsmax-1) 
    ddeltlist =  [1e-3, 1e-2, 1e-1, 1, 1e1]
    
        
    for d in ds.diameter:
        d = float(d)
        print(f'{d}', end = '\t')
        dsel = ds.sel(diameter = [d])
        if not np.isnan(float(dsel.ddelt[0])):
            continue
        for ddelt in ddeltlist:
            for ndgs in ndgs_list:
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
                            
        print(f'{ddelt}', end = ' ')                    
        print(f'{ndgs}', end = ' ')
        status = int((abs(np.log10(ddelt)) * 100) + ndgs)
        ds.status.loc[dict(diameter = d)] = status
        # ds.ndgs.loc[dict(diameter = d)] = ndgs    
        
        if not np.isnan(return_dict[d]):
            print(f'{return_dict[d]}',
                  # end = '',
                  )
            ds.scatt_cross_scect.loc[dict(diameter = d)] = return_dict[d]
        else:
            print('Fail!!',
                  # end = '',
                  )
        
            
            ds.ddelt.loc[dict(diameter = d)] = np.log10(ddelt)
            ds.ndgs.loc[dict(diameter = d)] = ndgs

    return ds




def optimize_optical_properties(parent_dict, d, pshape, n_real, n_imag,  verbose = False):
    """Recognizes if t-matrix approximation fails and addopts parameters until
    calculations are successfull, or all attempts fail."""
    
    manager = mp.Manager()
    return_dict = manager.dict()
    ndgsmax = 14
    ndgs_list = np.linspace(2, ndgsmax, ndgsmax-1) 
    ddeltlist =  [1e-3, 1e-2, 1e-1, 1, 1e1]
    
    d = float(d)
    if verbose:
        print(f'{d}', end = '\t')
    # dsel = ds.sel(diameter = [d])
    # if not np.isnan(float(dsel.ddelt[0])):
    #     continue
    for ddelt in ddeltlist:
        for ndgs in ndgs_list:
            process = mp.Process(target = calculate_optical_properties, 
                                  args = (return_dict,
                                          d,
                                          pshape,
                                          n_real,
                                          n_imag,
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
    if verbose:                    
        print(f'{ddelt}', end = ' ')                    
        print(f'{ndgs}', end = ' ')
    status = int((abs(np.log10(ddelt)) * 100) + ndgs)   
    # parent_dict[d] = return_dict[d]
    out = {'scatt': return_dict[d],
           'diameter': d,
           'pshape': pshape,
           'n_real': n_real,
           'n_imag': n_imag, 
           'status': status}
    parent_dict.append(out)
    return return_dict[d]
    # if not np.isnan(return_dict[d]):
    #     print(f'{return_dict[d]}',
    #           # end = '',
    #           )
    #     # ds.scatt_cross_scect.loc[dict(diameter = d)] = return_dict[d]
    # else:
    #     print('Fail!!',
    #           # end = '',
    #           )
    
        
    #     ds.ddelt.loc[dict(diameter = d)] = np.log10(ddelt)
    #     ds.ndgs.loc[dict(diameter = d)] = ndgs
    

           
