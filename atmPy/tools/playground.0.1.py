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

def runthis(d, return_dict):

    return_dict[d] = np.nan
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
                                  ndgs = 2,
                                  # orient = pytmatrix.orientation.orient_single,

                                 )
    scat = pytmatrix.scatter.sca_xsect(scatterer, 
                                       # h_pol=True,
                                       )
    res = scat / scale**2
    # results.append(res)
    # results['test'] = res
    # print(res)
    return_dict[d/scale] = res
    return res


if __name__ == "__main__":
    if 0:
        runthis()
    elif 1:
        mp.set_start_method('spawn')
        manager = mp.Manager()
        return_dict = manager.dict()
        for d in [5000, 500]:
            process = mp.Process(target = runthis, 
                                  args = (d,return_dict),
                                )
            
            process.start()
            process.join()
        print(return_dict)
            # print(process.exitcode)
        
        # else:
        #     process.run()
        
        # results
                
        
        
        