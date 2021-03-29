import numpy as _np
import pandas as _pd
import os as _os
import atmPy.general.timeseries as _timeseries
import atmPy.aerosols.physics.column_optical_properties as _column_optical_properties
import atmPy.general.measurement_site as _measurement_site
import pathlib as _pl
import warnings as _warnings

# from atmPy.general import measurement_site as _measurement_site
# from atmPy.radiation import solar as _solar

class DataSpecs(object):
    def __init__(self, parent):
        self.parent = parent
        self._window = None
        self._cloudscreened = True
        self.path = _pl.Path('/nfs/grad/surfrad/')
        
    def _reset_data(self):
        """reset all the data so it will be reloaded"""
        for station in self.parent.stations._stations_list:
            station.data._reset()
        
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
        p2f = self._parent.parent_network.data_specs.path.joinpath(param)

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
              'timezone': 8,
              'type': 'permanent'},
              {'name': 'Fort Peck',
              'state': 'MT',
              'abbreviation': ['FPK', 'fpk', 'FPE'],
              'lon': -105.10170,
              'lat': 48.30783,
              'alt': 634,
              'timezone': 7,
              'type': 'permanent'},
              {'name': 'Goodwin Creek',
              'state': 'MS',
              'abbreviation': ['GWN', 'gwn'],
              'lon': -89.8729,
              'lat': 34.2547,
              'alt': 98,
              'timezone': 6,
              'type': 'permanent'},
              {'name': 'Penn. State Univ.',
              'state': 'PA',
              'abbreviation': ['PSU', 'psu'],
              'lon': -77.93085,
              'lat': 40.72012,
              'alt': 376,
              'timezone': 5,
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


def _path2files(path2base_folder = '/nfs/grad/surfrad/aod/', site = 'bon', window = ('2020-07-31 00:00:00', '2020-10-30 23:00:00')):
    path2base_folder = _pl.Path(path2base_folder)

# first loading all files takes too long, so get site first
    path2aodatsite = path2base_folder.joinpath(site)
    assert(path2aodatsite.is_dir()), f'{path2aodatsite} not a directory'
    
    df = _pd.DataFrame([p for p in path2aodatsite.glob('**/*') if p.is_file()], columns=['path'])
    df.index = df.apply(lambda row: _pd.to_datetime(''.join([i for i in row.path.name if i.isdigit()])), axis = 1)
    df.sort_index(inplace=True)
    df = df.truncate(window[0], window[1])
    return df.path.values

def _path2files_deprecated(path2aod, site, window, perform_header_test, verbose):
    
    assert(not isinstance(window, str))

    if not isinstance(window, (tuple, list)):
        try: 
            window = [w.__str__() for w in window]
        except:
            raise ValueError(f'window has to be tuple or list, it is {type(window)}. Attempt to convert failed')
    
    if type(path2aod) == str:
        path2aod = _pl.Path(path2aod)

    elif type(path2aod) == list:
        for path in path2aod:
            assert (type(path).__name__ == 'PosixPath')
            assert (path.is_file())
        paths = path2aod

    elif type(path2aod).__name__ == 'PosixPath':
        pass

    else:
        raise ValueError('{} is nknown type for path2surfrad (str, list, PosixPath)'.format(path2aod))

    if type(path2aod).__name__ == 'PosixPath':
        if not (path2aod.is_file() or path2aod.is_dir()):
            raise ValueError('File or directory not found: {}'.format(path2aod))
        if path2aod.is_file():
            paths = [path2aod]
        elif path2aod.is_dir():
            paths = list(path2aod.glob('*'))
            keep_going = True
            while keep_going:
                keep_going = False
                for path in paths:
                    if path.is_file():
                        continue
                    elif path.is_dir():
                        dropped = paths.pop(paths.index(path))
                        paths += list(dropped.glob('*'))
                        keep_going = True
        else:
            raise ValueError
    files = paths

    if site:
        files = [f for f in files if site in f.name]
        if verbose:
            print('{} files match site specifications.'.format(len(files)))
    # select time window
    if window:
        start, end = window
        files = [f for f in files if (start.replace('-', '') <= f.name.split('_')[1].split('.')[0] and end.replace('-', '') > f.name.split('_')[1].split('.')[0])]
        if verbose:
            print('{} of remaining files are in the selected time window.'.format(len(files)))

    return files#, folder

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

def _read_files(files, verbose, UTC = False, cloud_sceened = True):

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

class Surfrad_AOD(_column_optical_properties.AOD_AOT):
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



def open_path(path = '/nfs/grad/surfrad/aod/',
              site = 'bon',
              window = ('2017-01-01', '2017-01-02'),
              cloud_sceened = True,
              local2UTC = False,
              perform_header_test = False,
              verbose = False,
              fill_gaps= False,
              keep_original_data = False):

    if site:
        if len([loc for loc in _locations if site in loc['abbreviation']]) == 0:
            raise ValueError('The site {} has not been set up yet. Add relevant data to the location dictionary'.format(site))

    files = _path2files(path, site, window)#, perform_header_test, verbose)
    file_content = _read_files(files, verbose, UTC=local2UTC, cloud_sceened=cloud_sceened)
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
    return saod