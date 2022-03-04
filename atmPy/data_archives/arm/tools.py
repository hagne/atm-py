import pandas as _pd
import os as _os
from pathlib import Path
import numpy as _np

def path2info(path, verbose = False):
    path = Path(path)
    if path.is_dir():
        path = list(path.iterdir())[0]

#suffix
    suffix = path.suffix
    suffixlist = ['.nc', '.cdf']
    if suffix not in suffixlist:
        raise ValueError('Suffix {} is unknown. Choose from {}'.format(suffix, suffixlist) )

    out = {}
    out['suffix'] = suffix

#site, product, facility
    site_product_facility = path.name.split('.')[0]
    idx_facil = _np.argmax(
        [i.isupper() for i in site_product_facility])  # index where the uppercase letter is, e.g. 9 for sgpaosapsE13
    idx_site = 3
    site, product, facility = site_product_facility[:idx_site], site_product_facility[idx_site:idx_facil], site_product_facility[
                                                                                                 idx_facil:]
    out['product'] = product
    out['facility'] = facility
    out['site'] = site

#timestamp
    out['timestamp'] = _pd.to_datetime(' '.join(path.name.split('.')[2:4]))

#qc_level
    out['qc_level'] = path.name.split('.')[1]

    if verbose:
        print('suffix: {}'.format(suffix))
        print('site: {}'.format(site))
        print('product: {}'.format(product))
        print('facility: {}'.format(facility))
        print('qc_level: {}'.format(out['qc_level']))
        # print(': {}'.format())
        # print(': {}'.format())
        # print(': {}'.format())
        # print(': {}'.format())
        # print(': {}'.format())
    return out


def path2filelist(path = '/Volumes/HTelg_4TB_Backup/arm_data/SGP/aosapsE13/', suffix = '.nc',
                  product = None, site = None, facility = None, window = None, qc_level = None,
                  # start_time = '2016-11-15', end_time = '2016-12-02'
                  verbose = False,
                  raise_error_when_empty = True,
                  ):
    """
    Takes a path and returns a list of files that meet various criteria defined by the kwargs
    Parameters
    ----------
    path: str
        Path to file or folder. If folder
    suffix: str ['.nc']
        file suffix
    product: str [None]
        which arm product to look for
    site: str ['sgp']
        which site ... SGP, NSA
    facility: str
        facility within the site, e.g. 'E13'
    window: tuple [None]
        e.g. '2016-11-15'

    Returns
    -------

    """

    path = Path(path)

    # start_time = _pd.to_datetime(start_time)
    # end_time = _pd.to_datetime(end_time)
    if not isinstance(window, type(None)):
        window = [_pd.to_datetime(t) for t in window]

    if path.is_dir():
        files = list(path.iterdir())
    else:
        files = [path]

    tfiles = []
    for file in files:
        info = path2info(file)
        tpath = type('TPath', (object,), {})
        tpath.path = file
        tpath.info = info
        tfiles.append(tpath)

    if verbose:
        print('Total files in folder: {}'.format(len(tfiles)))

    # only nc files
    # files = [file for file in files if file.suffix == suffix]
    tfiles = [file for file in tfiles if file.info['suffix'] == suffix]
    if verbose:
        print('Files with correct suffix: {}'.format(len(tfiles)))

    # select site
    # files = [file for file in files if site in file.name.split('.')[0]]
    if isinstance(site, type(None)):
        site = tfiles[0].info['site']
    tfiles = [file for file in tfiles if file.info['site'] == site]
    if verbose:
        print('Files with correct site: {}'.format(len(tfiles)))

    # select product
    # files = [file for file in files if product in file.name.split('.')[0]]
    if isinstance(product, type(None)):
        product = tfiles[0].info['product']
    tfiles = [file for file in tfiles if file.info['product'] == product]
    if verbose:
        print('Files with correct product: {}'.format(len(tfiles)))

    # select facility
    # files = [file for file in files if facility in file.name.split('.')[0]]
    if isinstance(facility, type(None)):
        facility = tfiles[0].info['facility']
    tfiles = [file for file in tfiles if file.info['facility'] == facility]
    if verbose:
        print('Files with correct facility: {}'.format(len(tfiles)))

    # select qc_level
    # files = [file for file in files if facility in file.name.split('.')[0]]
    if isinstance(qc_level, type(None)):
        qc_level = tfiles[0].info['qc_level']
    tfiles = [file for file in tfiles if file.info['qc_level'] == qc_level]
    if verbose:
        print('Files with correct qc_level: {}'.format(len(tfiles)))

    # files between start and end time
    if not isinstance(window, type(None)):
        # files = [file for file in files if start_time < _pd.to_datetime(' '.join(file.name.split('.')[2:4])) < end_time]
        tfiles = [file for file in tfiles if window[0] <= file.info['timestamp'] <= window[1]]
        if verbose:
            print('Files in window: {}'.format(len(tfiles)))

    # check if info is the same for all files
    tocheck = lambda tf: [tf.info[k] for k in ['suffix', 'product', 'site', 'facility', 'qc_level']]
    for tf in tfiles[1:]:
        assert(tocheck(tf) == tocheck(tfiles[0]))

    files = [file.path for file in tfiles]


    if len(files) == 0:
        txt = 'No files found that meet requirements!'
        if raise_error_when_empty:
            raise ValueError(txt)
        if verbose:
            print(txt)

    if verbose:
        print('Files to be opened:')
        for file in files:
            print('\t {}'.format(file))
    return files

#todo: this can be simplified with an internal function (from netCDF4 import num2date)
def _get_time(file_obj):
    bt = file_obj.variables['base_time']
    toff = file_obj.variables['time_offset']
    time = _pd.to_datetime(0) + _pd.to_timedelta(bt[:].flatten()[0], unit = 's') + _pd.to_timedelta(toff[:], unit = 's')
    return time


def is_in_time_window(f,time_window, verbose = False):
    out = True
    if time_window:
        fnt = _os.path.split(f)[-1].split('.')
        ts = fnt[-3]
        file_start_data = _pd.to_datetime(ts)
        start_time = _pd.to_datetime(time_window[0])
        end_time = _pd.to_datetime(time_window[1])
        dt_start = file_start_data - start_time
        dt_end = file_start_data - end_time
        out = file_start_data

        if dt_start.total_seconds() < -86399:
            if verbose:
                print('outside (before) the time window ... skip')
            out = False
        elif dt_end.total_seconds() > 86399:
            if verbose:
                print('outside (after) the time window ... skip')
            out = False
    return out