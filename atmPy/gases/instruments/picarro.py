import types
import pathlib as pl
import pandas as pd


def open_mobile(p2f=pl.Path('/home/ccg/mund/dev/mobileLabSites/bos/ML_dataManager/output/bos_picarro_60sec_avg.csv')):

    df = pd.read_csv(p2f)
    
    # set index
    df.index = pd.to_datetime(df.dt,)
    df.index.name = 'time'
    df.drop('dt', axis=1, inplace=True)
    
    dt = (df.index[1:] - df.index[:-1]).median()
    
    dtunits = ['days', 'hours', 'minutes', 'seconds']
    
    for dfu in dtunits:
        # print(getattr(dt, dfu))
        dtf = dt / pd.to_timedelta(1, dfu)
        if dtf >= 1:
            dt_freq = (dtf, dfu)
            break
    
    dt_freq
    
    ds = df.to_xarray()
    
    for v in ds.variables:
        if v == 'time':
            continue
        # print(v)
        unit = v.split('_')[-1]
        ds[v].attrs['unit'] = unit
        if v.split('_')[0] == 'avg':
            ds[v].attrs['avg'] = dt_freq
            
        ds = ds.rename({v: v.replace('avg', '').replace(unit, '').replace('_', '')})
        
    return Picarro(ds)


def open_prelim(path, verbose = False):
    p2f = pl.Path(path)

    if p2f.is_file():
        p2fs = [p2f]
    
    elif p2f.is_dir():
        p2fs = p2f.glob('*')
    
    elif isinstance(p2f, (types.GeneratorType, list, tuple)):
        p2fs = p2f
    
    else:
        assert(False), f'Type {type(p2fs)} not implemented'
    
    
    
    df = pd.DataFrame()
    for p2f in p2fs:
        if verbose:
            print(p2f.as_posix())
        dft = pd.read_csv(p2f)
        df = pd.concat([df,dft])
    
    df.index = pd.to_datetime(df.TMSTAMP)
    df.sort_index(inplace=True)
    df.index.name = 'time'
    coords = df.loc[:,['Latitude', 'Longitude', 'Elevation']].iloc[0]
    df = df.drop(['TMSTAMP','Frac_DOY', 'Frac_time', 'DOY', 'Year', 'Hour', 'Min','Latitude', 'Longitude', 'Elevation'], axis = 1)
    
    ds = df.to_xarray()
    
    ds.attrs['coordinates'] = '; '.join([f'{i[0]}:{i[1]}' for i in zip(coords.index, coords.values)])
    return Picarro(ds)


class Picarro(object):
    def __init__(self, dataset):
        self.dataset = dataset