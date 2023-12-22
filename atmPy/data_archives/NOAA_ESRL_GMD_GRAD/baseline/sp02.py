#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 13:25:59 2023

@author: hagen
"""
import numpy as np
import magic
import xarray as xr
import pandas as pd
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.baseline.baseline as atmbl
import warnings
import pathlib as pl
import types

class BaselineDatabase(object):
    def __init__(self):
        self.line_table = pd.DataFrame(columns=['site','install', 'line', 'instrument_id', 'comment'])
        self.instrument_table = pd.DataFrame(columns = ['instrument_id','type_id', 'sn', 'config'])

    def add2line_table(self, site,install_datetime, line_id, instrument_id, comment = ''):#20200205
        install_datetime = pd.to_datetime(install_datetime)
        new_line_table_entry = pd.DataFrame({'site':site,'install': install_datetime, 'line': line_id, 'instrument_id': instrument_id, 'comment': comment}, index = [instrument_id])
        # self.line_table = self.line_table.append(new_line_table_entry, ignore_index=True)
        self.line_table = pd.concat([self.line_table, new_line_table_entry], ignore_index=True)
        return
    
    def addnewinstrument(self, instrument_id, type_id, sn, config):
        # self.instrument_table =  self.instrument_table.append({'instrument_id': instrument_id,'type_id':type_id, 'sn': sn, 'config':config_id}, ignore_index=True)
        # new_instrument = pd.DataFrame({'instrument_id': instrument_id,'type_id':type_id, 'sn': sn, 'config':config}, index = [instrument_id])
        new_instrument = pd.DataFrame( [[instrument_id, type_id, sn, config]], columns = ['instrument_id', 'type_id', 'sn', 'config'],index = [instrument_id])
        self.instrument_table = pd.concat([self.instrument_table, new_instrument])#, ignore_index=True)
        return 
    
    def get_instrument(self, site, line, date, when_instrument_not_installed = 'error'):
        """
        

        Parameters
        ----------
        site : TYPE
            DESCRIPTION.
        line : TYPE
            DESCRIPTION.
        date : TYPE
            DESCRIPTION.
        when_instrument_not_installed : string, optional ('error', 'warn', 'silent')
            When silent or warn None is returned. The default is 'error'.

        Returns
        -------
        None.

        """
        fail = False
        lt_site_line = self.line_table[np.logical_and(self.line_table.site == site, self.line_table.line == line)]
        # previous_installs = lt_site_line[lt_site_line.install <= date]
        
        # last install of an instrument at that site on that line
        diff = (lt_site_line.install - date)         
        diff_last = diff[(diff/pd.to_timedelta(1, 'd')) < 0]
        if len(diff_last) == 0:
            txt = f'Instrument not installed at site and line (line:{line}, site: {site}, date: {date})'
            fail = True
        else:
            
            last_install = self.line_table.loc[diff_last.idxmax()]
            diff_last = diff[(diff/pd.to_timedelta(1, 'd')) < 0]
                   
            # has that instrument been installed somewhere else since that install?
            # last install of the particular instrument ... same as above?
            instrument_installs = self.line_table[self.line_table.instrument_id == last_install.instrument_id]
            diff = (instrument_installs.install - date)
            diff_last = diff[(diff/pd.to_timedelta(1, 'd')) < 0]
            last_install_of_inst = self.line_table.loc[diff_last.idxmax()]
            
            # if previous_installs.shape[0] == 0:
            if last_install.name != last_install_of_inst.name:
                txt = (f'Instrument not installed at site and line (line:{line}, site: {site}, date: {date}\n'
                       f'At the given date the instrument was at {last_install_of_inst.site} on line {last_install_of_inst.line} (since {last_install_of_inst.install})'
                       )
                fail = True
            
        if fail:
            if when_instrument_not_installed == 'error':
                raise IndexError(txt)
            elif when_instrument_not_installed == 'warn':
                warnings.warn(txt)
                return None
            elif when_instrument_not_installed == 'silent':
                return None
            else:
                assert(False), f'{when_instrument_not_installed} not an option for when_instrument_not_installed. (error, warn,or silent)'
                
        lt_found = last_install #previous_installs.iloc[-1]

        instrument_found = self.instrument_table[self.instrument_table.instrument_id == lt_found.instrument_id].iloc[0]
        return instrument_found

#### THE DATABASE
database = BaselineDatabase()
database_global = database
#### filter comfigurations
conf_1= {'A': 368, 'B': 1050, 'C': 610, 'D': 778}
conf_2= {'A': 412, 'B': 500, 'C': 675, 'D': 862}

#### Instruments 
database.addnewinstrument(1,1,1032,conf_2)
database.addnewinstrument(2,1,1046,conf_1)
database.addnewinstrument(3,1,1022,conf_2) #typically at SPO

#### instrument linups
installdate = '20131126'
database.add2line_table('mlo', installdate, 121, 2)
database.add2line_table('mlo', installdate, 221, 1)

installdate = '20140104' # something is statring to go wrong on that day!
database.add2line_table('mlo', installdate, 121, 1)
database.add2line_table('mlo', installdate, 221, 2)

installdate = '20141204' 
database.add2line_table('mlo', installdate, 121, 2)
database.add2line_table('mlo', installdate, 221, 1)

installdate = '20151203' 
database.add2line_table('mlo', installdate, 121, 1)
database.add2line_table('mlo', installdate, 221, 2)

installdate = '20161211' 
database.add2line_table('mlo', installdate, 121, 1)
database.add2line_table('mlo', installdate, 221, 2)

installdate = '20171207' 
database.add2line_table('mlo', installdate, 121, 2)
database.add2line_table('mlo', installdate, 221, 1)

#### BRW 2018 lines assigned
installdate = '20180710' # No idea why it was installed so late? 
#uninstalldate = '20191019' # Was never uninstalled in 2018
database.add2line_table('brw', installdate, 121, 1)
database.add2line_table('brw', installdate, 221, 2)
# database.add2line_table('trans', uninstalldate, 0, 1)
# database.add2line_table('trans', uninstalldate, 0, 2)

#### BRW 2019
#installdate = '2019' # Was never uninstalled in 2018
uninstalldate = '20191019' # was taken down a little later, but there where fluctuations in temperature that render the data invalid. 
# database.add2line_table('brw', installdate, 121, 1)
# database.add2line_table('brw', installdate, 221, 2)
database.add2line_table('trans', uninstalldate, 0, 1)
database.add2line_table('trans', uninstalldate, 0, 2)

### MLO 2020

database.add2line_table('mlo', '20200205', 121, 1)
database.add2line_table('mlo', '20200205', 221, 2)

database.add2line_table('mlo', '20200620', 121, 2)
database.add2line_table('mlo', '20200620', 221, 1)

#### BRW 2020 lines assigned
installdate = '20200717' 
uninstalldate = '20201222' 
database.add2line_table('brw', installdate, 121, 1)
database.add2line_table('brw', installdate, 221, 2)
database.add2line_table('trans', uninstalldate, 0, 1)
database.add2line_table('trans', uninstalldate, 0, 2)

#### MLO 2021
installdate = '20210204' 
database.add2line_table('mlo', installdate, 121, 1)
database.add2line_table('mlo', installdate, 221, 2)

### BRW 
installdate = '20210318' 
uninstalldate = '20211008' 
database.add2line_table('brw', installdate, 121, 1)
database.add2line_table('brw', installdate, 221, 2)
database.add2line_table('trans', uninstalldate, 0, 1)
database.add2line_table('trans', uninstalldate, 0, 2)

installdate = '20220101' 
database.add2line_table('mlo', installdate, 121, 1)
database.add2line_table('mlo', installdate, 221, 2)

installdate = '20220309' 
database.add2line_table('mlo', installdate, 121, 3)


def read_file(path2raw, lines = None, 
              # collabels = ['lineID', 'Year', 'DOY', 'HHMM', 'A', 'B', 'C', 'D','temp'],
              collabels = ['lineID', 'Year', 'DOY', 'HHMM'],
              database = None,
              site = None,
              when_instrument_not_installed = 'error',
              verbose = False,
              test = False
              ):
    """
    The particular way I am reading here allows for later implementation of
    reading old data from Longenecker. And also allows to read other raw files

    Parameters
    ----------
    path2raw : str, list, pathlib.Path
        Single or list of path(s) to file(s).
    lines : list, optional
        List of lines to consider (e.g. [121, 221] for sp02 at MLO). The default is None -> all.
    collabels : TYPE, optional
        DESCRIPTION. The default is ['lineID', 'Year', 'DOY', 'HHMM'].
    database : TYPE, optional
        DESCRIPTION. The default is None.
    site : str, optional
        DESCRIPTION. The default is None. If None the site is infered from the
        file path. Set if the path is not the standard path
    when_instrument_not_installed: str ['error', 'skip', 'ignore']:
        What to do if the instrument is not installed according to the database.
        error: ...da!
        skip: will skip processing current file and move to the next one, this is most usefull in operation conditions
        ignore: Continues processing. All instrument related information (including calibrations) will be ignored. this is most usefull when in the process of finding out about installation lines and times.
        

    Returns
    -------
    out_list : TYPE
        DESCRIPTION.

    """ 
    assert(when_instrument_not_installed in ['error', 'skip', 'ignore'])
    if when_instrument_not_installed == 'error':
        womi = 'error'
    elif when_instrument_not_installed in ['skip', 'ignore']:
        womi = 'silent'
    else:
        assert(False), 'haeaeaeaea?'
        
    out = {}
    collabels = np.array(collabels)
    
    if isinstance(database, type(None)):
        database = database_global
    
    #### File list
    if isinstance(path2raw, str):
        path2raw = [pl.Path(path2raw),]
    elif isinstance(path2raw, types.GeneratorType):
        path2raw = list(path2raw)
    
    if not isinstance(path2raw, list):
        path2raw = [path2raw,]
    
    path2raw.sort()
    
    if verbose:
        print(f'opening {len(path2raw)} files.')
    
    assert(np.all(['sp02' in p2f.name for p2f in path2raw])), 'Not all files contains "sp02" which indicates that you are trying to open a non sp02 file. Consider adding "sp02" to your name pattern. A common file name is gradobs.brw-sp02-1.20200614.raw.dat'
    assert(len(path2raw) > 0), "There are no files to open. Adjust you file name pattern.\nHere is a typical file name: gradobs.brw-sp02-1.20200614.raw.dat"
    firstfile = path2raw[0]
    if magic.from_file(firstfile.as_posix()) == 'Hierarchical Data Format (version 5) data':
        ds = xr.open_mfdataset(path2raw) 
        if ('product' not in ds.attrs) or ds.attrs['product'] == 'SP02 raw':
            # return SP02RawData(ds)
            return ds

    else:
        sites = [p2r.name.split('.')[1].split('-')[0] for p2r in path2raw]
        sites = [''.join(char for char in site if not char.isdigit()) for site in sites]
        site = np.unique(sites)
        if len(site) == 1:
            site = site[0]
            site_inst = atmbl.network.stations.find_site(site)
            if verbose:
                print(f'found site: {site_inst.name}')
        else:
            assert(False), f'More than one site found in file names, {site}'
            
        df_all = pd.concat([pd.read_csv(p2r, header=None) for p2r in path2raw])
        
        # df_all = pd.read_csv(path2raw, header=None
        # #             names = False
        #            )
        colsis = df_all.columns.values
        colsis = [int(col) for col in colsis]
        
        # todo: assigne collumn labels accoreding to instrument info
        # if 0:
        colsis[:collabels.shape[0]] = collabels
        df_all.columns = colsis   
        # out['df_all_02'] = df_all.copy()
        # df_all = pd.read_csv(path2raw, names=lines[0]['column_labels'])
    
        # make datetime index
        df_all['HHMM'] = df_all.apply(lambda row: f'{int(row.HHMM):04d}', axis=1)
        df_all.index = df_all.apply(lambda row: pd.to_datetime(f'{int(row.Year)}0101') + pd.to_timedelta(row.DOY - 1 , 'd') + pd.to_timedelta(int(row.HHMM[:2]), 'h') + pd.to_timedelta(int(row.HHMM[2:]), 'm'), axis=1)
        df_all.index.name = 'datetime'
        
        # data_list = []
        # df_inst_temp = pd.DataFrame()
        # df_inst_channels = pd.DataFrame()
        out['df_all'] = df_all.copy()
        # return out
        out_list = []
        # date = df_all.index[0]
        # print(df_all.lineID.unique())
        for lid in df_all.lineID.unique():
            if verbose:
                print(f'line id: {lid}')
            if isinstance(lines, list):
                if lid not in lines:
                    print(f'{lid} not in lines ({lines})')
                    continue
                
            df_lid = df_all[df_all.lineID == lid].copy()
            
            # there was the case that Longenecker must have created an overlab when stiching two days together ... therefor ->
            df_lid = df_lid[~df_lid.index.duplicated()]
            df_lid.sort_index(inplace=True)
            
            
            instrument = database.get_instrument(site, lid, df_all.index[0],
                                                 when_instrument_not_installed = 'warn')
            if isinstance(instrument, type(None)):
                instrument = database.get_instrument(site, lid, df_all.index[-1],
                                                     when_instrument_not_installed = womi)
            if when_instrument_not_installed == 'skip':
                if isinstance(instrument, type(None)):
                    continue
            
            df_lid = df_lid.drop(['lineID', 'Year','DOY', 'HHMM'], axis=1)
            df_lid.columns = ['A', 'B', 'C', 'D', 'temp']
        
            # replace invalid values with nan
            df_lid[df_lid == -99999] = np.nan
            df_lid[df_lid == -7999.0] = np.nan
        
            # seperate photo detector readings from temp 
            df_temp = df_lid.temp
            df_voltages = df_lid.reindex(['A', 'B', 'C', 'D'], axis= 1)
        
            df_voltages.columns.name = 'channel'
        
            # create dataset
            ds = xr.Dataset()
            ds['direct_normal_irradiation'] = df_voltages #this is the raw data
            ds['instrument_temperature'] = df_temp
            if not isinstance(instrument , type(None)):
                ser = pd.Series(instrument.config)
                ser.index.name = 'channel'
                ds['channel_center'] = ser
                sn = instrument['sn']
                ds.attrs['serial_no'] = sn            
                # assign a nominal wavelength as the channel name. currently nominal 
                # wavelength is equal to center wavelength. At some point it would be 
                # nice to have a measured wavelength in addition
                ds['channel_letter'] = ds.channel
                ds = ds.assign_coords(channel = ds.channel_center)  
                
            # ds.attrs['serial_no'] = sn
            ds.attrs['info'] = 'This is raw sp02 data'
            # ds.attrs['site'] = site
            # sn = 'brw'
            ds.attrs['site_latitude'] = site_inst.lat
            ds.attrs['site_longitude'] = site_inst.lon
            ds.attrs['site_elevation'] =  site_inst.alt
            ds.attrs['site_name'] = site_inst.name
            ds.attrs['site'] =site
            ds.attrs['instrument_type'] = 'sp02'
            ds.attrs['line_id'] = lid
            
  
            
            out_list.append(ds)
            

        if test:
            out['out_list'] =out_list
            return out
        else:
            return out_list