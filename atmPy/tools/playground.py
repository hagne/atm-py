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

def runthis(return_dict, d, ndgs, ddelt):
    print('.', end = '', flush=True)
    return_dict[d] = np.nan
    do = d
    scale = 1e-12
    wl = 500 * scale
    # d = .500 #
    ratio = 3.
    # d = 2822* scale #limit with ndgs = 5
    # d = 500 * scale
    d = d * scale
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


if __name__ == "__main__":
    if 0:
        runthis()
    elif 1:
        mp.set_start_method('spawn')
        manager = mp.Manager()
        return_dict = manager.dict()
        diameters = [4000,
                  # 2822, 
                  # 500,
                  ]
        diameters = np.logspace(np.log10(3050), np.log10(15000), 100)  # converged for all numbers, many at  
        ndgsmax = 14
        ndgs_list = np.linspace(2, ndgsmax, ndgsmax-1) 
        ddeltlist =  [1e-3, 1e-2, 1e-1, 1, 1e1]
        
        # just for testing
        test = False
        if test:
            diameters = [3050]
            ndgsmax = 20
            ndgs_list = np.linspace(2, ndgsmax, ndgsmax-1) 
            ddeltlist = [0.1]
            
        for d in diameters:
            print(f'\n{d}', end = ' ')
            
            # ddeltlist = [1e-3,1e-4,1e-5,]
            for ddelt in ddeltlist:
                print(f'\n\t{ddelt}', end = ' ')
                res_list = []
                for ndgs in ndgs_list:
                    # ndgs = i + 1
                    print(ndgs, end = ' ')
                    process = mp.Process(target = runthis, 
                                          args = (return_dict,
                                                  d,
                                                  ndgs,
                                                  ddelt,
                                                  ),
                                        )
                    process.start()
                    process.join()
                    
                    if test:
                        print(f'\t{return_dict[d]}')
                        if not np.isnan(return_dict[d]):
                            res_list.append(return_dict[d])
                    elif np.isnan(return_dict[d]):
                        pass
                    else:
                        print(f'\t{return_dict[d]}', end = ' ')
                        break
                    
                if not np.isnan(return_dict[d]):
                    break
        if test:
            mean = np.mean(res_list)
            std = np.std(res_list)
            print(f'mean: {mean}')
            print(f'std: {std}')
            print(f'rel. std {(1 - (mean - std)/mean) *100}')
            
        print('')
        print(return_dict)
            # print(process.exitcode)
        
        # else:
        #     process.run()
        
        # results
                
        
        
        