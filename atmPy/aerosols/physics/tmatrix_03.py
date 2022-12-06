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
import pathlib as pl


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
                                  orient = pytmatrix.orientation.orient_averaged_fixed,
                                  or_pdf = pytmatrix.orientation.uniform_pdf(),
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
    
def generate_lut(diameters, pshapes, n_real, n_imag):
    data = np.zeros([diameters.shape[0], pshapes.shape[0], n_real.shape[0], n_imag.shape[0]])
    data[:] = np.nan
    
    ds = xr.Dataset()
    da = xr.DataArray(data.copy(), coords={'diameter': diameters, 'pshape': pshapes, 'n_real': n_real, 'n_imag': n_imag})
    ds['scatt_cross_scect'] = da
    ds['status'] =            da.copy()
    return ds
    
if __name__ == "__main__":
    #### settings
    p2fld = pl.Path('/home/grad/htelg/projects/uncertainty_paper/lut03')
    
    p2out = 'lut03_{d_shape}_{ps_shape}_{nr_shape}_{ni_shape}_ch{chunk}.nc'

    no_cpu = 30
    chunksize = no_cpu * 20 #save at the end of every chunk, 
    iterations_d = 0
    
    assert(pl.isdir())
    #### ---- limites
    
    #### diamters
    # according to the paper the expeced diameter varies by pm 20% for 2sigma -> 3 sigma are 30%, to make sure it works lets do 40%
    
    aps_min = 750 * 0.6
    aps_max = 15000 * 1.4
    
    diameters = get_samplingpoints(aps_min, aps_max, scale='log', round2=3, iteration=iterations_d)
    
    #### refractive index real
    
    # for coarse:
    # The assumption of quartz and feldspar as the predominant constituents of dust particles would suggest a refractive
    # 216index between 1.5 and 1.6. However, experiments on ambient particles suggest larger values of up to 1.67 (Eidhammer
    # 217et al., 2008; Ishida et al., 1991). In addition, it has been shown that n of coarse mode particles has a significant
    # 218imaginary part of up to 0.01 (Eidhammer et al., 2008). Varying n of coarse mode particles from 1.5 to 1.67i0.01
    
    # 0.05/2 * 3
    
    sigma3 = 0.075
    mean = 1.55
    n_real = get_samplingpoints(mean - sigma3, mean + sigma3, scale='log', round2=3, iteration=1)
      
    # For accu!!!!!!!
    # Instead of
    # 210considering the case where each σ calculation is using an individual n based on the particular ACSM measurement
    # 211we assume a constant value of n = 1.5, which we find for the mean value of all recorded values in 2012. For the
    # 212uncertainty we use δn = 0.028, which are two standard deviations from the mean.
    
    # That means 3 stds are 0.042
    
    sigma3 = 0.042
    mean = 1.5
    
    n_real = get_samplingpoints(mean - sigma3, mean + sigma3, scale='log', round2=2, iteration=0)
    
    #### refractive index imaginary
    
    n_imag = get_samplingpoints(0.001, 0.01, scale='log', round2=1, iteration=1)
    n_imag[0] = 0
    
    ### particle shape
    
    spmin = -5
    spmax = 5

    pshapes = get_samplingpoints(spmin, spmax, scale='lin', round2=3, iteration=1)
    
    for e,psi in enumerate(pshapes):
        if psi <0:
            pshapes[e] = abs(psi)+1
        elif psi>= 0:
            pshapes[e] = 1/(psi + 1)
        else:
            assert(False)
            
    pshapes = round_sig(pshapes, 2)
    pshapes
    
    #### combine to empty LUT
    
    print(f'diameters: {diameters.shape[0]}')
    print(f'pshapes: {pshapes.shape[0]}')
    print(f'n_real: {n_real.shape[0]}')
    print(f'n_imag: {n_imag.shape[0]}')
    no_combinations = n_real.shape[0] * n_imag.shape[0] * diameters.shape[0] * pshapes.shape[0]
    print(f'all combos: {no_combinations}')
    
    ds = generate_lut(diameters, pshapes, n_real, n_imag)
    
    #### ---- run it
    
    # p2out = f'lut_{diameters.shape[0]}_{pshapes.shape[0]}_{n_real.shape[0]}_{n_imag.shape[0]}' + '_ch{}.nc'

    #### Merge with previeous dataset
    # ds = xr.merge([ds_old, ds])
    
    ds_stack = ds.stack(dim = ['diameter', 'pshape', 'n_real', 'n_imag'])
    print(ds_stack.status.shape)
    ds_stack = ds_stack.status[ds_stack.status.isnull()]
    print(ds_stack.shape)
    
    #### Split into chunk. At the end of these chunks the data is saved
    i = 0
    while True:
        start = i * chunksize
        if start >= ds_stack.shape[0]:
            break
        print('{',end = '')
        i += 1
        end = i * chunksize
        ds_stack_chunk = ds_stack[start: end]
        
        #### split up into setsthis is the set that is parallized
        e = 0
        while True:
            start = e * no_cpu
            if start >= ds_stack_chunk.shape[0]:
                break
            print('.',end = '')
            e += 1
            end = e * no_cpu
            ds_stack_chunk_sel = ds_stack_chunk[start:end]
    
            # assert(False)
            manager = mp.Manager()
            res_list = manager.list()
            subproslist = []
            for params in ds_stack_chunk_sel.dim:
                # break
    
                d = float(params.diameter)
                pshape = float(params.pshape)
                n_r = float(params.n_real)
                n_i = float(params.n_imag)
                # if not np.isnan(ds.loc[dict(diameter = d, pshape = pshape, n_real = n_r, n_imag = n_i)].status):
                #     continue
                assert(np.isnan(ds.loc[dict(diameter = d, pshape = pshape, n_real = n_r, n_imag = n_i)].status)), 'I would have thought that all status != nan ar removed ?'
                process = mp.Process(target = optimize_optical_properties, 
                                      args = (res_list,
                                              d,
                                              pshape,
                                              n_r,
                                              n_i,
                                              ),
                                    )
                process.start()
                # process.join()
                subproslist.append(process)
    
            [p.join() for p in subproslist]
    
            #### assign results to the dataset
            results = list(res_list)
            while 1:
                res = results.pop(0)
                scatt = res.pop('scatt')
                status = res.pop('status')
                ds.scatt_cross_scect.loc[res] = scatt
                ds.status.loc[res] = status
                if len(results) == 0:
                    break
        p2o_t = p2out.format(d_shape = diameters.shape[0],
                               ps_shape = pshapes.shape[0],
                               nr_shape = n_real.shape[0],
                               ni_shape = n_imag.shape[0],
                               chunk = i)
        ds.to_netcdf(p2fld.joinpath(p2o_t))
    print('done')
           
