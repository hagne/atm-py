# -*- coding: utf-8 -*-
import xarray as _xr
from atmPy.data_archives.goes_satellites import satlab

def open_file(fname, verbose = False):
    if verbose:
        print(f'opening {fname}')
    ds = _xr.open_dataset(fname)
    if 'OR_ABI-L2-MCMIPC-M6' in ds.dataset_name:
        if verbose:
            print(f'data porduct:  OR_ABI-L2-MCMIPC-M6')
        out = satlab.OR_ABI_L2_MCMIPC(ds)
    elif 'OR_ABI-L2-AODC-M6' in ds.dataset_name:
        if verbose:
            print(f'data porduct:  OR_ABI-L2-AODC-M6')
        out = satlab.OR_ABI_L2_AODC_M6(ds)
    return out