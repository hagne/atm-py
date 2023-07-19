import numpy as _np
import pandas as _pd
# import os as _os
import atmPy.general.timeseries as _timeseries
import atmPy.aerosols.physics.column_optical_properties as _column_optical_properties
import atmPy.general.measurement_site as _measurement_site
import atmPy.atmosphere.sounding as atmsound
import pathlib as _pl
# import warnings as _warnings
import xarray as _xr
import atmPy.radiation.observations.spectral_irradiance as atmradobs
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.cal_facility.lab as atmcucf

# from atmPy.general import measurement_site as _measurement_site
# from atmPy.radiation import solar as _solar

class DataSpecs(object):
    def __init__(self, parent):
        self.parent = parent
        self._window = None
        self._cloudscreened = True
        self._path = _pl.Path('/nfs/grad/')
         
    def _reset_data(self):
        """reset all the data so it will be reloaded"""
        for station in self.parent.stations._stations_list:
            station.data._reset()
    
    @property
    def path(self):
        return self._path
    
    @path.setter
    def path(self, value):
        self._reset_data()
        self._path = _pl.Path(value)
    
    @property
    def window(self):
        return self._window
    
    @window.setter
    def window(self,value):
        self._reset_data()
        self._window = _pd.to_datetime(value)
        
    @property
    def cloud_screened(self):
        return self._cloudscreened
    
    @cloud_screened.setter
    def cloud_screened(self, value):
        self._reset_data()
        self._cloudscreened = value

class Surfrad_Data(object):
    def __init__(self, parent):
        self._parent = parent
        self._aod = None
    
    def _reset(self):
        self._aod = None
    
    @property
    def aod(self):
        if isinstance(self._aod, type(None)):
            self._aod = self.load_data(param = 'aod', )
        return self._aod
            
    def load_data(self, param = 'aod'):
        # get_data function
        
        # check requried parameters
        assert(not isinstance(self._parent.parent_network.data_specs.window, type(None))), "set data_specs"

        window = self._parent.parent_network.data_specs.window

        # assert(_np.unique([w.year for w in window]).shape[0] == 1), 'Data range spannes over more than one year ... programming required'

        # # get path
        # p2f = self._parent.parent_network.data_specs.path.joinpath(param)
        p2f = self._parent.parent_network.data_specs.path
        if self._parent.abb.lower() == 'bnd':
            siteabb = 'bon'
        else:
            siteabb = self._parent.abb.lower()

        # p2f = p2f.joinpath(siteabb)

        # p2f = p2f.joinpath(str(window[0].year))

        # assert(p2f.is_dir()), f'not an existing path: {p2f}'

        # load data
        if param == 'aod':
            mfrsr = open_path(path = p2f, site = siteabb, window = window, local2UTC=True, 
                              cloud_sceened = self._parent.parent_network.data_specs.cloud_screened,)
        else:
            assert(False), f'the param:{param} is not set up yet ... programming required'
        return mfrsr
    
    
_locations = [{'name': 'Bondville',
              'state' :'IL',
              'abbreviation': ['BND', 'bon'],
              'lon': -88.37309,
              'lat': 40.05192,
              'alt' :230,
              'timezone': -6},
              {'name': 'Sioux Falls',
              'state': 'SD',
              'abbreviation': ['SXF', 'sxf'],
              'lon': -96.62328,
              'lat': 43.73403,
              'alt': 473,
              'timezone': -6},
              {'name': 'Table Mountain',
              'state': 'CO',
              'abbreviation': ['TBL', 'tbl'],
              'lon': -105.23680,
              'lat': 40.12498,
              'alt': 1689,
              'timezone': -7},
              {'name': 'Desert Rock',
              'state': 'NV',
              'abbreviation': ['DRA', 'dra'],
              'lon': -116.01947,
              'lat': 36.62373,
              'alt': 1007,
              'timezone': -8,
              'type': 'permanent'},
              {'name': 'Fort Peck',
              'state': 'MT',
              'abbreviation': ['FPK', 'fpk', 'FPE'],
              'lon': -105.10170,
              'lat': 48.30783,
              'alt': 634,
              'timezone': -7,
              'type': 'permanent'},
              {'name': 'Goodwin Creek',
              'state': 'MS',
              'abbreviation': ['GWN', 'gwn'],
              'lon': -89.8729,
              'lat': 34.2547,
              'alt': 98,
              'timezone': -6,
              'type': 'permanent'},
              {'name': 'Penn. State Univ.',
              'state': 'PA',
              'abbreviation': ['PSU', 'psu'],
              'lon': -77.93085,
              'lat': 40.72012,
              'alt': 376,
              'timezone': -5,
              'type': 'permanent'},
              # {'name': 'ARM Southern Great Plains Facility',
              # 'state': 'OK',
              # 'abbreviation': ['SGP', 'sgp'],
              # 'lon': -97.48525,
              # 'lat': 36.60406,
              # 'alt': 314,
              # 'timezone': 6,
              # 'type': 'permanent'},
              # {'name': '',
              #  'state': '',
              #  'abbriviations': ['', ''],
              #  'lon': -,
              #  'lat': ,
              #  'alt': ,
              #  'timezone': ,
              #  'type': 'permanent'}
              ]

_channel_labels = _np.array([415, 500, 614, 673, 870, 1625])

network = _measurement_site.Network(_locations)
network.name = 'surfrad'
network.data_specs = DataSpecs(network)
for station in network.stations._stations_list:
    station.data = Surfrad_Data(station)

#todo: remove the following dictionary ... its not used anymore
_deprecated_col_label_trans_dict = {'OD413': 415,
                         'OD414': 415,
                         'OD415': 415,
                         'OD416': 415,
                         'OD417': 415,
                         'AOD417':415,
                         'OD495': 500,
                         'OD496': 500,
                         'OD497': 500,
                         'OD499': 500,
                         'OD500': 500,
                         'OD501': 500,
                         'OD502': 500,
                         'OD503': 500,
                         'OD504': 500,
                         'OD505': 500,
                         'OD609': 614,
                         'OD612': 614,
                         'OD614': 614,
                         'OD615': 614,
                         'OD616': 614,
                         'AOD616':614,
                         'OD664': 673,
                         'OD670': 673,
                         'OD671': 673,
                         'OD672': 673,
                         'OD673': 673,
                         'OD674': 673,
                         'OD676': 673,
                         'OD861': 870,
                         'OD868': 870,
                         'OD869': 870,
                         'AOD869':870,
                         'OD870': 870,
                         'OD1623': 1625,
                         'OD1624': 1625,
                         ### sometime AOD was used instead of OD =>
                         'AOD414': 415,
                         'AOD415': 415,
                         'AOD497': 500,
                         'AOD501': 500,
                         'AOD609': 614,
                         'AOD615': 614,
                         'AOD664': 673,
                         'AOD673': 673,
                         'AOD861': 870,
                         'AOD870': 870,
                         }


def _path2files(path2base_folder = '/nfs/grad/surfrad/aod/', site = 'bon', window = ('2020-07-31 00:00:00', '2020-10-30 23:00:00'), product = 'aod'):
    path2base_folder = _pl.Path(path2base_folder)

# first loading all files takes too long, so get site first
    path2aodatsite = path2base_folder.joinpath(site)
    if product == 'albedo':
        path2aodatsite = path2aodatsite.joinpath(product)
        product = 'alb'
    assert(path2aodatsite.is_dir()), f'{path2aodatsite} not a directory'
    
    df = _pd.DataFrame([p for p in path2aodatsite.glob(f'**/*.{product}') if p.is_file()], columns=['path'])
    df.index = df.apply(lambda row: _pd.to_datetime(''.join([i for i in row.path.name if i.isdigit()])), axis = 1)
    # df.sort_index(inplace=True)
    df = df.sort_index().truncate(window[0], window[1])
    return df.path.values


def _read_header(fname):
    """Read the header of file in folder and reterns a dict with relevant data"""
    header_size = 5
    with fname.open() as myfile:
        head = [next(myfile) for x in range(header_size)]

    out = {}
    # header size
    out['header_size'] = header_size
    # site
    out['site'] = head[0].split()[0]
    # channels
    channels = head[2].split()
    out['channels'] = channels[:channels.index('channel')]

    # date
    out['date'] = _pd.to_datetime(head[1].split()[0])
    #     return head
    return out

def _read_data(fname, UTC = False, header=None):
    """Reads the file takes care of the timestamp and returns a Dataframe
    """
    if not header:
        header = _read_header(fname)

    # dateparse = lambda x: _pd.datetime.strptime(x, "%d:%m:%Y %H:%M:%S")
    df = _pd.read_csv(fname, skiprows=header['header_size'],
                     delim_whitespace=True,
                     #                      na_values=['N/A'],
                     #                   parse_dates={'times': [0, 1]},
                     #                   date_parser=dateparse
                     )

    datetimestr = '{0:0>4}{1:0>2}{2:0>2}'.format(header['date'].year, header['date'].month, header['date'].day)+ df.ltime.apply \
        (lambda x: '{0:0>4}'.format(x)) + 'UTC'  # '+0000'
    df.index = _pd.to_datetime(datetimestr, format="%Y%m%d%H%M%Z")
    if UTC:
        try:
            timezone = [l for l in _locations if header['site'] in l['name']][0]['timezone']
        except IndexError:
            try:
                timezone = [l for l in _locations if header['site'] in l['abbreviation']][0]['timezone']
            except IndexError:
                raise ValueError('Site name {} not found in _locations (neither in name not in abbreviation)'.format(header['site']))
        df.index += _pd.to_timedelta(-1 * timezone, 'h')
        df.index.name = 'Time (UTC)'
    else:
        df.index.name = 'Time (local)'
    # for col in df.columns:
    #     if col in ['ltime', '0=good', 'p_mb', 'Ang_exp']:
    #         continue
    #     elif 'E' in col:
    #         continue
    #     elif col not in _col_label_trans_dict.keys():
    #         print('not in transdict: {}'.format(col))

    trans_dict = {}
    for key in df.columns:
        if not (('OD' in key) or ('E' in key)):
            continue
        wl = int(key.replace('AOD', '').replace('OD', '').replace('E', ''))
        col_lab = _channel_labels[abs(_channel_labels - wl).argmin()]
        #     print('{} -> {}'.format(key, col_lab))
        #     print('{}\t -> {}\t -> {}\t -> {}'.format(key, wl, col_lab, key.replace(str(wl), str(col_lab))))
        if 'E' in key:
            col_lab = '{}E'.format(col_lab)
        trans_dict[key] = col_lab


    df.rename(columns=trans_dict, inplace=True)
    return df

def _read_aod_files(files, verbose, UTC = False, cloud_sceened = True):

    if len(files) == 0:
        raise ValueError('no Files to open')

    if verbose:
        print('Reading files:')
    data_list = []
    wl_match_list = []
    header_first = _read_header(files[0])
    for fname in files:
        if verbose:
            print('\t{}'.format(fname), end=' ... ')
        header = _read_header(fname)

        # make sure that all the headers are identical
        if header_first['site'] != header['site']:
            try:
                site = [site for site in _locations if header_first['site'] in site['name']][0]
            except IndexError:
                try:
                    site = [site for site in _locations if header_first['site'] in site['abbreviation']][0]
                except IndexError:
                    raise ValueError(
                        'Site name {} not found in _locations (neither in name not in abbreviation)'.format(
                            header['site']))

            if  header['site'] in site['abbreviation']:
                header['site'] = header_first['site']
            elif header['site'] in site['name']:
                header_first['site'] = header['site']
                # _warnings.warn('The site name changed from {} to {}! Since its the same site we march on.'.format(header_first['site'], header['site']))
            else:
                raise ValueError('The site name changed from {} to {}!'.format(header_first['site'], header['site']))
        # read the data
        data = _read_data(fname, UTC = UTC, header=header)
        data_list.append(data)

        # matching table that gives the exact wavelength as a function of time (identical for
        # cols = [col for col in data.columns if str(col).isnumeric()]
        cols = [int(str(col).replace('AOD', '').replace('OD', '')) for col in data.columns if
                str(str(col).replace('AOD', '').replace('OD', '')).isnumeric()]

        wls = header['channels']
        try:
            wl_match_list.append(_pd.DataFrame([wls] * data.shape[0], columns=cols, index=data.index))
        except:
            pass

        if verbose:
            print('done')

    # concatinate and sort Dataframes and create Timeseries instance
    data = _pd.concat(data_list, sort=False)
    data.sort_index(inplace=True)
    data[data == -999.0] = _np.nan
    data[data == -9.999] = _np.nan
    data = _timeseries.TimeSeries(data, sampling_period=1 * 60)
    if cloud_sceened:
        data.data[data.data['0=good'] == 1] = _np.nan
    out = {'data': data}
    out['header_first'] = header_first

    # concatenate wavelength match dataframe
    wl_match = _pd.concat(wl_match_list,sort=False)
    wl_match.sort_index(inplace=True)
    wl_match = wl_match.astype(float)
    wl_match = _timeseries.TimeSeries(wl_match, sampling_period=60)
    out['wavelength_match'] = wl_match

    if verbose:
        print('done')

    return out

class Surfrad_AOD(_column_optical_properties.AOD_AOT_20221216):
    pass
#     def __init__(self, lat, lon, elevation = 0, name = None, name_short = None, timezone = 0):
#         self._aot = None
#         self._aot = None
#         self._sunposition = None
#         self._timezone = timezone
#
#         self.site = _measurement_site.Site(lat, lon, elevation, name=name, abbreviation=name_short)
#
#
#     @property
#     def sun_position(self):
#         if not self._sunposition:
#             if self._timezone != 0:
#                 date = self.AOD.data.index +  _pd.to_timedelta(-1 * self._timezone, 'h')
#             else:
#                 date = self.AOD.data.index
#             self._sunposition = _solar.get_sun_position(self.site.lat, self.site.lon, date)
#             self._sunposition.index = self.AOD.data.index
#             self._sunposition = _timeseries.TimeSeries(self._sunposition)
#         return self._sunposition
#
#     @property
#     def AOT(self):
#         if not self._aot:
#             if not self._aod:
#                 raise AttributeError('Make sure either AOD or AOT is set.')
#             aot = self.AOD.data.mul(self.sun_position.data.airmass, axis='rows')
#             aot.columns.name = 'AOT@wavelength(nm)'
#             aot = _timeseries.TimeSeries(aot)
#             self._aot = aot
#         return self._aot
#
#     @ AOT.setter
#     def AOT(self,value):
#         self._aot = value
#
#     @property
#     def AOD(self):
#         if not self._aod:
#             if not self._aot:
#                 raise AttributeError('Make sure either AOD or AOT is set.')
#             aod = self.AOT.data.dif(self.sun_position.data.airmass, axis='rows')
#             aod.columns.name = 'AOD@wavelength(nm)'
#             aod = _timeseries.TimeSeries(aod)
#             self._aod = aod
#         return self._aod
#
#     @ AOD.setter
#     def AOD(self,value):
#         self._aod = value

def read_albedo(path2file, path2readme = '/nfs/grad/Inst/MFR/README.alb', verbose = True):
    """
    Read a path or list of paths to albedo files

    Parameters
    ----------
    path2file : str, pathlib.Path, list
        Where is the data: /nfs/grad/Inst/MFR/SURFRAD/tbl/albedo/2021/.
    path2readme : TYPE, optional
        This readme needs to be read to proved the column names. The default is '/nfs/grad/Inst/MFR/README.alb'.
    verbose : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    p2r = _pl.Path(path2readme)
    if p2r.name.split('.')[-1] == 'preliminary':
        filetype = 'preliminary'
    elif p2r.name.split('.')[-1] == 'alb':
        filetype = 'albedo'

    def read_readme(p2r):
        with open(p2r, 'r') as rein:
            lines = rein.readlines()

        clines = [l for l in lines if l[0].isnumeric()]
        colnames = [l.split(' - ')[1] for l in clines]
        colnames = [c.strip() for c in colnames]
        return colnames



    def read_data(p2f, colnames):
        df = _pd.read_csv(p2f, delim_whitespace=True,
                    names = colnames,
                   )


        df.index = df.apply(lambda row: _pd.to_datetime(str(int(row.YYYYMMDDhhmmss)), format = '%Y%m%d%H%M%S'), axis = 1)
        df.index.name = 'datetime'

        # remove datetime collumns
        df = df.iloc[:,5:].copy()

        # get the albedo and channels
        alb_cols = [c for c in df.columns if 'Albedo' in c]
        alb_cols = [c for c in alb_cols if not 'Broadband' in c]
        alb_df = df.loc[:,alb_cols]
        alb_df[alb_df == -9.999] = _np.nan
        # return {"alb_df": alb_df}
        if filetype == 'preliminary':
            wl_df = _pd.Series({e:int(c.split(':')[1].strip().split()[0]) for e,c in enumerate(alb_df.columns)}).transpose()
        elif filetype == 'albedo':
            wl_df = _pd.Series({e:int(c.split(':')[1].strip().replace('615nm/','').replace('nm','')) for e,c in enumerate(alb_df.columns)}).transpose()
        wl_df.index.name = 'channel'
        alb_df.columns = wl_df.index

        # df without albedo and columnames change

        dfwoa = df.drop(alb_cols, axis=1)
        dfwoa.rename({'NDVI "Normailized Difference Vegetation Index"': 'NDVI'}, axis = 1, inplace=True)
        dfwoa.columns = [c.replace(' ', '_') for c in dfwoa.columns]

        # change dtypes
        dfwoa = dfwoa.astype(_np.float32)
        alb_df = alb_df.astype(_np.float32)

        # make dataset
        ds = _xr.Dataset(dfwoa)

        # add albedo and wavelength
        ds['Albedo'] = alb_df

        ds['Wavelength_nominal'] = wl_df
        
        ds = ds.rename_vars({v:v.replace('/','_') for v in ds.variables if '/' in v})

        return ds
    

    if isinstance(path2file, (list, _np.ndarray, tuple)):
        pass
    elif _pl.Path(path2file).is_dir():
        assert(False), 'just a little bit of programming needed to make this work!!'
    elif isinstance(path2file, (str, _pl.Path)):
        path2file = [path2file,]
    else:
        TypeError(f"Don't know what to do with type: {type(path2file).__name__}")
        
    colnames = read_readme(p2r) 
    ds = read_data(path2file[0], colnames)
    # return ds
    ds = _xr.concat([read_data(fn, colnames) for fn in path2file], dim = 'datetime')
    ds['Wavelength_nominal'] = ds.Wavelength_nominal[0].drop('datetime')
    return ds

        
def get_mfrsr_filter_responds(serial_no, path2folder = '/nfs/grad/Calibration_facilities/cucf/Surfrad/working'):
    p2fld_filter = _pl.Path(path2folder)
    p2f_filter_resp = list(p2fld_filter.glob(f'*{serial_no:04d}*SPR.txt'))
    assert(len(p2f_filter_resp) == 1), f'there should only be one responds function!! {len(p2f_filter_resp)} found'
    p2f_filter_resp = p2f_filter_resp[0]

    filter_resp = atmcucf.read_mfrsr_cal(p2f_filter_resp)
    return filter_resp


def read_ozon(p2f):
    """
    This reads the ozon files ({site}_ozone.dat) in /home/grad/surfrad/aod/ozone/

    Parameters
    ----------
    p2f : TYPE
        DESCRIPTION.

    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    """
    df = _pd.read_csv(p2f, names=['year', 'month','day','epoch', 'ozone'])
    df.index = df.apply(lambda row: _pd.to_datetime(f'{row.year:04d}{row.month:02d}{row.day:02d}'), axis = 1)
    df = df.drop(['year', 'month', 'day', 'epoch'], axis=1)
    df.sort_index(inplace=True)
    df.index.name = 'datetime'
    ds = df.to_xarray()
    ds.ozone.attrs['unit'] = 'dobson'
    ds.attrs['info'] = "Ozone absorption is accessed from a file that is written by a Perl script that gets daily ozone over the station being processed from NASA TOMS OMI or OMPS.(from documentation in aod_analysis.f)"
    return ds

def read_sounding(p2f):
    """
    This reads the interpolated functions that John's code uses and which are
    stored here: /nfs/grad/surfrad/sounding/

    Parameters
    ----------
    p2f : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # outlist = []
    # ds_sonde = xr.Dataset()
    dt_sonde = _pd.to_datetime(p2f.name, format = '%Y%m%d_%H.int')
    
    with open(p2f, 'r') as rein:
        lines = rein.readlines()

    header = lines[0]
    lines = lines[1:]
    
    sounding_no = int(header.split()[0])
    
    chsize = 39
    sitechuncks = [lines[i * chsize: (i+1) * chsize] for i in range(sounding_no)]
    
    for e,sc in enumerate(sitechuncks):
        # break
        # print(sc[0].strip().split(' '))
        site = ' '.join(sc[0].strip().split(' ')[:2]).strip()
        # print(site)
        # tp_1 = sc.copy()
        sc = sc[1:]
    
        data = [[float(v) for v in l.strip().split()] for l in sc]
    
        columns = ['pressure', 'height', 'temperature', 'dewpoint', 'uwind', 'vwind']
    
        data = _pd.DataFrame(data, columns=columns)
    
        data[data == -999] = _np.nan
    
        data.dropna(inplace=True)
    
        # data.loc[:,'temperature']= data.temperature + 273.16
        # tc = data.temperature - 273.16
        # data['esat'] = 6.1078 * np.exp((17.2693882*tc)/(tc+237.3))
        # data['water_massfrac'] = 621.97*data.esat/(data.pressure-data.esat)
    
        # tpw = spint.simps(data[::-1].water_massfrac, data[::-1].pressure)/980 #precipitable water in cm, no idea what that 980 is for?
    
        # outlist.append({'site': site, 'tpw': tpw})
        
        
        # data.index = data.pressure
        ds = data.to_xarray()
        ds = ds.expand_dims({'site': [site,]})
        
        if e == 0:
            ds_sonde = ds
        else:
            ds_sonde = _xr.concat([ds_sonde, ds], 'site')
        # if site == 'Goodwin':
        #     break
    
    # tpwdf = pd.DataFrame(outlist)
    ds_sonde.attrs['datetime'] = dt_sonde
    
    return atmsound.BalloonSounding(ds_sonde)

def read_ccc(p2f, verbose = False):
    """
    Reads surfrads ccc (cosine corrected) files typically located somewhere 
    around here:
    /nfs/grad/Inst/MFR/SURFRAD/dra/mfrsr/ccc/2019/dra20190728_0669.ccc

    Parameters
    ----------
    p2f : TYPE
        DESCRIPTION.

    Returns
    -------
    ds : TYPE
        DESCRIPTION.

    """
    with open(p2f, 'r') as rein:
        for i in range(2):
            l = rein.readline()

    filetype = 'ccc'
    instrument_type = 'MFRSR'
    fieldno = len(l.split())
    
    if verbose:
        print(f'fieldno: {fieldno}')
    
    if fieldno == 36:
        filetype = 'tu'
    elif fieldno == 18:
        filetype = 'tu'
        instrument_type = 'MFR'
    
    if verbose:
        print(f'filetype: {filetype}')
        print(f'instrument_type: {instrument_type}')
        
    #### parse name
    site_abb = p2f.name[:3]
    serial_no = int(p2f.name.split('_')[-1].split('.')[0])

    site = [s for s in network.stations._stations_list if s.abb == site_abb.upper()][0]
    
    cols = list(range(fieldno))
    cols[0] = 'datetime'

    if filetype == 'tu':
        skiprows = 1
        timezoneadjust = site.time_zone['diff2UTC_of_standard_time']
    else:
        skiprows = None
        timezoneadjust = 0
    df = _pd.read_csv(p2f,
                names =cols,
                delim_whitespace = True,
                skiprows = skiprows)

    index =  df.apply(lambda row: _pd.to_datetime('1900-01-01') + _pd.to_timedelta(row.datetime - 1, 'days'), axis = 1) # the -1 is bacause the year starts with day 1 not zero
    index -= _pd.to_timedelta(timezoneadjust, 'h')
    df.index = index
    df.sort_index(inplace = True)

    df.index.name = 'datetime'

    ds = _xr.Dataset()

    # return df
    if filetype == 'ccc':
        ds['solar_elevation'] = df[1]
        ds['instrument_temperature'] = df[23]
        ds['instrument_power'] = df[24]

    elif filetype == 'tu':
        ds['temp_sensor'] = df[1]
        
        ds['temp_housing_2'] = df.iloc[:,-5]
        da = _xr.DataArray(df[2])
        da = da.where(da != -9999.0, _np.nan)
        da = da.where(da != -9998.0, _np.nan)
        ds['broadband'] = da#df[2]
        no = 7
        if instrument_type == 'MFRSR':
            # ds['temp_housing_2'] = df.iloc[:,-5]
            start_col = 2
        elif instrument_type == 'MFR':
            ds['temp_housing'] = df[10]
            # ds['temp_housing_2'] = df.iloc[:,-5]
            start_col = 2
        else:
            assert(False), 'nenenenene'
            
        sel = df.iloc[:,start_col: start_col + no].copy()

        #### TODO: assign a QF invalid due to ...
        sel[sel == -9999.0] = _np.nan
        sel[sel == -9998.0] = _np.nan

        sel.columns = _np.array(range(len(sel.columns))) + 1
        sel.columns.name = 'channel'
        var_name = 'all_time'
        ds[var_name] = sel
        ds[var_name].attrs['unit'] = 'mV'
        # ds[var_name].attrs['corrections'] = 'cosine' It not corrected!!!!
    else:
        assert(False), 'not possible!'
    groupnames = ['global_horizontal_irradiation', 'diffuse_horizontal_irradiation', 'direct_normal_irradiation']

    if filetype == 'tu':
        start_col = 2+7+2
    else:
        start_col = 2
    for block in range(3):
        if instrument_type == 'MFR':
            break # there is only the alltime collumn in the MFR?
        # block = 0
        no = 7
        sel = df.iloc[:,start_col + (block * no): start_col + (block * no) + no].copy()

        #### TODO: assign a QF invalid due to ...
        sel[sel == -9999.0] = _np.nan
        sel[sel == -9998.0] = _np.nan

        sel.columns = _np.array(range(len(sel.columns))) + 1
        sel.columns.name = 'channel'

        ds[groupnames[block]] = sel
        ds[groupnames[block]].attrs['unit'] = 'mV'
        ds[groupnames[block]].attrs['corrections'] = 'cosine'

    if filetype == 'tu':
        sel = df.iloc[:, -4:].copy()
        ds['airmass'] = sel.iloc[0]
        ds['solar_azimuth'] = sel.iloc[1]
        ds['solar_zenith'] = sel.iloc[2]
        ds['solar_elevation'] = sel.iloc[3]
        # sel.columns = ['airmass', 'azimuth', 'zenith', 'elevation']
        # sel.columns.name = 'sun_pos_params'
        # ds['sun'] = sel

    ds.solar_elevation.attrs['unit'] = 'radian'
    ds.attrs['info'] = 'Cosine corrected SURFRAD MFRSR measurments.'
    
    ds.attrs['site_latitude'] = site.lat
    ds.attrs['site_longitude'] = site.lon
    ds.attrs['site_elevation'] =  site.alt
    ds.attrs['site_name'] = site.name
    ds.attrs['site'] =site_abb
    ds.attrs['instrument_type'] = instrument_type
    ds.attrs['serial_no'] = serial_no
    
    # name channels and clean up
    fresp = get_mfrsr_filter_responds(serial_no)
    ds = ds.drop_sel({'channel':1}) #this is the broadband channel
    ds = ds.assign_coords(channel = fresp.channel) 
    ds['channel_center'] = fresp.statistics.sel(stats = 'CENT').to_pandas()
    # return ds
    out = atmradobs.CombinedGlobalDiffuseDirect(ds)
    out.filter_functions = fresp
    return out

def open_path(path = '/nfs/grad/',
              product = 'aod',
              site = 'bon',
              window = ('2017-01-01', '2017-01-02'),
              cloud_sceened = False,
              local2UTC = False,
              perform_header_test = False,
              verbose = False,
              fill_gaps= False,
              keep_original_data = False,
              test = False):
    """
    

    Parameters
    ----------
    path : TYPE, optional
        Path to data. If directory this is considered the base directory 
        (without site). If this is a file name, only this file will be opened.
        The default is '/nfs/grad/surfrad/aod/'.
    site : TYPE, optional
        DESCRIPTION. The default is 'bon'.
    window : TYPE, optional
        DESCRIPTION. The default is ('2017-01-01', '2017-01-02').
    cloud_sceened : TYPE, optional
        DESCRIPTION. The default is False.
    local2UTC : TYPE, optional
        DESCRIPTION. The default is False.
    perform_header_test : TYPE, optional
        DESCRIPTION. The default is False.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.
    fill_gaps : TYPE, optional
        DESCRIPTION. The default is False.
    keep_original_data : TYPE, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    saod : TYPE
        DESCRIPTION.

    """
    
    path = _pl.Path(path)
        
    if path.is_file():
        files = [path]
        site = path.name[:3]
    else:
        if product == 'aod':
            path = path.joinpath('surfrad/aod/')
        elif product == 'albedo':
            path = path.joinpath('Inst/MFR/SURFRAD/')
        files = _path2files(path, site, window, product)#, perform_header_test, verbose)
    if test:
        return files
    
    if files[0].name.split('.')[-1] == 'alb':
        if verbose:
            print('File type: albedo')
        assert(_np.all([f.name.split('.')[-1] == 'alb' for f in files])), 'Not all files are of the same type (albedo).'
        
        ds = read_albedo(files)
        return ds
        
    elif files[0].name.split('.')[-1] == 'aod':
        if verbose:
            print('File type: AOD')
        assert(_np.all([f.name.split('.')[-1] == 'aod' for f in files])), 'Not all files are of the same type (AOD).'
        file_content = _read_aod_files(files, verbose, UTC=local2UTC, cloud_sceened=cloud_sceened)   
        data = file_content['data']
        wl_match = file_content['wavelength_match']
        header_first = file_content['header_first']
        if fill_gaps:
            if verbose:
                print('filling gaps', end=' ... ')
            data.data_structure.fill_gaps_with(what=_np.nan, inplace=True)
            wl_match.data_structure.fill_gaps_with(what = _np.nan, inplace = True)
            if verbose:
                print('done')
    
        # add Site class to surfrad_aod
        if site:
            if len([loc for loc in _locations if site in loc['abbreviation']]) == 0:
                raise ValueError('The site {} has not been set up yet. Add relevant data to the location dictionary'.format(site))
        try:
            site = [l for l in _locations if header_first['site'] in l['name']][0]
        except IndexError:
            try:
                site = [l for l in _locations if header_first['site'] in l['abbreviation']][0]
            except IndexError:
                raise ValueError('Looks like the site you trying to open ({}) is not set up correctly yet in "location"'.format(header_first['site']))
    
        lon = site['lon']
        lat = site['lat']
        alt = site['alt']
        timezone = site['timezone']
        site_name = site['name']
        abb = site['abbreviation'][0]
        # saod.site = _measurement_site.Site(lat, lon, alt, name=site_name, abbreviation=abb)
    
        # generate Surfrad_aod and add AOD to class
        saod = Surfrad_AOD(lat = lat, lon = lon, elevation = alt, name=site_name, name_short=abb, timezone = timezone, site_info = site)
        if keep_original_data:
            saod.original_data = data
    
        if local2UTC:
            saod._timezone = 0
        else:
            saod._timezone = timezone
        ## select columns that show AOD
        data_aod = data.drop(_channel_labels, inverse=True)
    
        ## rename columns
        data_aod.data.columns.name = 'AOD@wavelength(nm)'
        data_aod.data.sort_index(axis = 1, inplace=True)
    
        ## add the resulting Timeseries to the class
        saod.AOD = data_aod
    
        #Angstrom exponent
        saod.ang_exp = data.drop('Ang_exp', inverse=True)
    
        # wavelength matching table
        saod.wavelength_matching_table = wl_match
    
        # errors
        saod.AODerrors = data.drop([col for col in data.data.columns if 'E' in str(col)], inverse = True)
        saod.AODerrors.data.columns = [int(col.replace('E', '')) for col in saod.AODerrors.data.columns]
        saod.files_opened = files
        
        saod.cloudmask.cloudmask_nativ = data.drop('0=good', inverse = True)
        saod.aerosolmask.cloudmask_nativ = data.drop('0=good', inverse = True)
        return saod
    else:
        assert(False), 'noenoenoe'
