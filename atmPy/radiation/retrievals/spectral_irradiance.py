#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 16:30:35 2022

@author: hagen
"""
import numpy as np
import pandas as pd
import xarray as xr
import atmPy.general.measurement_site as atmms
import atmPy.radiation.retrievals.langley_calibration as atmlangcalib
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.cal_facility.lab as atmcal
import atmPy.radiation.rayleigh.lab as atmraylab
import pathlib as pl
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.surfrad.surfrad as atmsrf
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.baseline.baseline as atmbl

import sqlite3
import matplotlib.pyplot as _plt
import copy
import types
from .. import solar as atmsol
import typing

class RadiationDatabase(object):
    def __init__(self, path2db = '/nfs/grad/surfrad/database/surfraddatabase.db', 
                 # parent = None, 
                 verbose = True):
        if not isinstance(path2db, pl.Path):
            path2db = pl.Path(path2db)
        self.path2db = path2db
        self.verbose = verbose
        # self.parent = parent
        
        self._connection = None

    def initiat_new_db(self, si):
        """
        Initiate a new database so tables with the right format are present. This might be out of date since
        
        Arguments
        ==========
        si: SolarIrradiation instance (or inheritances)
        
        """
        
        site_id = self.add_site(si.dataset.site, si.dataset.site_name , 
                     si.dataset.site_elevation, 
                     si.dataset.site_latitude,
                     si.dataset.site_longitude, overwrite=True)
        
        si.add_file_to_database(self)
        
        
        self.add_instrument_type(
            'mfrsr',
            'Multi Filter Rotating Shadowband Radiomenter',
            overwrite=False)
        
        for channel in si.dataset.channel:
            channel = int(channel)
            self.add_channel(channel, 'mfrsr', 'optical_narrowband', overwrite=True)

        return
    
    @property
    def connection(self):
        if isinstance(self._connection, type(None)):
            self._connection = sqlite3.connect(self.path2db)
        # check if connection has been closed ... or not alive for other reasons
        else:
            try:
                self._connection.cursor()
            except sqlite3.ProgrammingError as e:
                self._connection = sqlite3.connect(self.path2db)

        return self._connection

    def add_site(self, abb, name, elevation, latitude, longitude, overwrite = False, closeconnection=True):
        ### add new site
        table_name = 'sites'
        try:
            index = self.connection.execute(f"SELECT max(site_id)+1 FROM {table_name}").fetchall()[0][0]
            if isinstance(index, type(None)):
                index = 1
            exists = pd.read_sql_query(f"SELECT * FROM sites Where abb = '{abb}'", self.connection)
            if exists.shape[0] != 0:  # 'site with that abbriviation already exists'
                index = exists.site_id.iloc[0]
                if overwrite:
                    self.connection.execute(f'delete from sites where site_id={index}').fetchall()
                else:
                    return index
        except Exception as e:
            # Only needed on the very first run, when there is no sites yet
            if e.args[0] == f'no such table: {table_name}':
                index = 1
            else:
                raise
        
        df = pd.DataFrame({'abb': [abb,], 
                           'name': [name,], 
                           'elevation': [elevation,],
                           'latitude': [latitude,],
                           'longitude': [longitude,]},
                          index = [index])
        df.index.name = 'site_id'
        df.to_sql(table_name, self.connection,
                  if_exists = 'append'
                 )
        if closeconnection:
            self.connection.close()
        return index
    
    
    def add_file(self, 
                    path, date=None, date_created=None, site=None, 
                    product_name = None, product_version = None, overwrite = True):
        if date_created in ['NA', 'None']:
            date_created = 'NA'
        else:
            dt = pd.to_datetime(date_created)  #ds.creation_timestamp)
            date_created = f'{dt.year:04d}-{dt.month:02d}-{dt.day:02d}'
        dt = pd.to_datetime(date)  #ds.datetime.values[0])
        date = f'{dt.year:04d}-{dt.month:02d}-{dt.day:02d}'
        site_id = self.connection.execute(f"SELECT site_id FROM sites Where abb='{site}'").fetchall()[0][0]
        table_name = 'files'
        try:
            index = self.connection.execute(f"SELECT max(file_id)+1 FROM {table_name}").fetchall()[0][0]
            exists = pd.read_sql_query(f"SELECT * FROM files Where path = '{path.as_posix()}'", self.connection)
            if exists.shape[0] != 0:
                index = exists.file_id.iloc[0]
                if overwrite:
                    self.connection.execute(f'delete from files where file_id={index}').fetchall()
                else:
                    return index
        except Exception as e:
            # Only needed on the very first run, when there is no sites yet
            if e.args[0] == f'no such table: {table_name}':
                index = 1
            else:
                raise
        
        df = pd.DataFrame({'site_id': [site_id,],
                           'name': [path.name,],
                           'date': [date,],
                           'date_created': [date_created,],
                           'product_name': [product_name,],
                           'product_version': [product_version,],
                           'path': [path.as_posix(),],
                           },
                          index = [index])
        df.index.name = 'file_id'
        df.to_sql(table_name, self.connection, 
                  if_exists='append',
                 )
        return index
    
    
    def add_instrument_type(self, abb, name, overwrite = False):
        ### add new instrument
        table_name = 'instrument_type'
        try:
            index = self.connection.execute(f"SELECT max(instrument_id)+1 FROM {table_name}").fetchall()[0][0]
            exists = pd.read_sql_query(f"SELECT * FROM {table_name} Where abb = '{abb}'", self.connection)
            if exists.shape[0] != 0:
                index = exists.instrument_id.iloc[0]
                if overwrite:
                    self.connection.execute(f'delete from {table_name} where instrument_id={name}').fetchall()
                else:
                    return index
        except Exception as e:
            # Only needed on the very first run, when there is no sites yet
            if e.args[0] == f'no such table: {table_name}':
                index = 1
            else:
                raise
        
        df = pd.DataFrame({'abb': [abb,], 
                           'name': [name,],
                          },
                          index = [index])
        df.index.name = 'instrument_id'
        df.to_sql(table_name, self.connection, 
                  if_exists = 'append',
                 )
        return index
    
    ### add new wavelength channel
    
    def add_channel(self,
                    name,
                    instrument,
                    channel_type,
                    overwrite = False):
    
        instrument_id = self.connection.execute(f"SELECT instrument_id FROM instrument_type Where abb='{instrument}'").fetchall()[0][0]
        table_name = 'channels'
        try:
            index = self.connection.execute(f"SELECT max(channel_id)+1 FROM {table_name}").fetchall()[0][0]
            exists = pd.read_sql_query(f"SELECT * FROM {table_name} Where name = {name}", self.connection)
            if exists.shape[0] != 0:
                index = exists.channel_id.iloc[0]
                if overwrite:
                    self.connection.execute(f'delete from {table_name} where name={name}').fetchall()
                else:
                    return index
        except Exception as e:
            # Only needed on the very first run, when there is no sites yet
            if e.args[0] == f'no such table: {table_name}':
                index = 1
            else:
                raise
        
        df = pd.DataFrame({'instrument_id': [instrument_id,], 
                           'name': [name,],  
                           'channel_type': [channel_type,], 
                          },
                          index = [index])
        df.index.name = 'channel_id'
        df.to_sql(f'{table_name}', self.connection, 
                  if_exists = 'append',
                 )
        return index

    def add_clearsky_params(self, df,
                            overwrite = True):
        conn = self.connection
        dft = pd.DataFrame(index = [0,])
        try:
            file_id = pd.read_sql_query(f"SELECT * FROM files where path='{df.file.iloc[0].as_posix()}'",
                                        conn).file_id.iloc[0]
        except IndexError:
            raise IndexError('File not registerd yet. Do so with SunIrradiation.register_file_in_database.')
        dft['file_id'] = file_id
        site_id = pd.read_sql_query(f"SELECT * FROM sites where abb='{df.site.iloc[0]}'", conn).site_id.iloc[0]
        dft['site_id'] = site_id
        dft['channel_id'] = pd.read_sql_query(f"SELECT * FROM channels where name='{df.channel.iloc[0]}'", conn).channel_id.iloc[0]
        
        df = df.drop(['file', 'site','channel'], axis = 1)
        df = pd.concat([dft, df], axis = 1)
        
        table_name = 'clearsky_params'
        try:
            index = self.connection.execute(f"SELECT max(fit_id)+1 FROM {table_name}").fetchall()[0][0]
            exists = pd.read_sql_query(f'''select * from clearsky_params 
                                              where file_id={df.file_id.iloc[0]}
                                                  and channel_id={df.channel_id.iloc[0]}
                                                  and observation="{df.observation.iloc[0]}"
                                              ''', self.connection)
            if exists.shape[0] != 0:
                index = exists.fit_id.iloc[0]
                if overwrite:
                    self.connection.execute(f'''delete from clearsky_params 
                                              where file_id={df.file_id.iloc[0]}
                                                  and channel_id={df.channel_id.iloc[0]}
                                                  and observation="{df.observation.iloc[0]}"
                                              ''').fetchall()
                else:
                    return index
        except Exception as e:
            # Only needed on the very first run, when there is no sites yet
            if e.args[0] == f'no such table: {table_name}':
                index = 1
            else:
                raise
        
        df.index = [index]
        df.index.name = 'fit_id'
        df.to_sql(f'{table_name}', self.connection, 
                  if_exists = 'append',
                 )
        self.connection.close()
        return index
                


class ClearSky(object):
    def __init__(self,parent):
        self.parent = parent
        self._testresults = None

    @property
    def test_results(self):
        if isinstance(self._testresults, type(None)):
            ds = self.parent.dataset.copy()
            if 'sun_position' not in ds.variables:
                try:
                    df = self.parent.direct_normal_irradiation.sun_position
                    df.columns.name = 'sun_params'
                    ds['sun_position'] = df.drop('ampm', axis = 1)
                    # return df
                except:
                    print('sun_position is not in dataset variables. I tried to get it, but it came back with the following error message')
                    raise
            
            clearskylist = []
            for channel in ds.channel:
                
                dst = xr.Dataset()
                for param in ['global_horizontal', 'diffuse_horizontal', 'direct_normal']:
                    # break
                    if param not in ds.variables:
                        continue
                    # Do a fit as a function of airmass. This is therefore a symmetrical fit for pm and am and should quickly show if there is a problem with am and pm
                    # todo make a test that indicates symmetry problems
                    amlim = 12 # maximum airmass for visualization etc, 
                    airmass = ds.sun_position.sel(sun_params = 'airmass')
                    diffuse = ds[param].sel(channel = channel).dropna('datetime')
                    diffuse = diffuse.where(airmass < amlim, drop = True).where(airmass > 1, drop = True)
                    airmass = airmass.where(~diffuse.isnull())
                    
                    amlim = 10 # maximum airmass to which to consider the clearsky analysis
                    diffusesel = diffuse.where(airmass < amlim, drop = True).where(airmass > 1, drop = True)
                    airmasssel = airmass.where(~diffusesel.isnull(), drop = True)
                
                    degree = 3 # degree of polynomial fit; 2 might also be enough
                    coeff = np.polyfit(airmasssel.values, diffusesel.values, degree)
                    coeffsym = coeff
                    diffpred = np.poly1d(coeff)(airmass)
                    # aa[1].plot(airmass,diffuse - diffpred, label = degree)
                    
                    # Apply the fit from above to the entire day. Take the difference and analyse
                    data = diffuse - diffpred
                    dtfirst = data.datetime.values[0]
                    x = (data.datetime - dtfirst)/pd.to_timedelta(1,'s')
                    
                    coeff = np.polyfit(x, data, 5)
                    datapred = np.poly1d(coeff)(x)
                    dst[f'{param}_clearsky_symmetric'] = xr.DataArray(np.array([diffpred,]), dims = ['channel', 'datetime'], coords = {'channel': [channel], 'datetime': data.datetime})
                    dst[f'{param}_clearsky'] = xr.DataArray(np.array([diffpred + datapred,]), dims = ['channel', 'datetime'], coords = {'channel': [channel], 'datetime': data.datetime})
                    dst[f'{param}_asymmetry_fit'] = xr.DataArray(np.array([datapred,]), dims = ['channel', 'datetime'], coords = {'channel': [channel], 'datetime': data.datetime})
                    dst[f'{param}_asymmetry'] = xr.DataArray(np.array([data,]), dims = ['channel', 'datetime'], coords = {'channel': [channel], 'datetime': data.datetime})
                    dst[f'{param}_asymmetry_fit_params'] = xr.DataArray(np.array([coeff,]), coords = {'channel': [channel], 'polynom_deg': range(6)[::-1]})
                    dst[f'{param}_symmetric_fit_params'] = xr.DataArray(np.array([coeffsym,]), coords = {'channel': [channel], 'polynom_deg': range(4)[::-1]})
                    #### Add some clearsky quality values to it which hopefully help identifying clear sky days?
                    # relative or normalized RMSE (NRMSE)
                    nrmse_sym = np.sqrt(((diffuse-diffpred) ** 2).mean())/diffuse.mean()
                    nrmse_clearsky = np.sqrt(((diffuse-(diffpred + datapred)) ** 2).mean())/diffuse.mean()
                    nrmse_sym = nrmse_sym.expand_dims({'clearsky_quality_params': ['nrmse_sym']})
                    nrmse_clearsky = nrmse_clearsky.expand_dims({'clearsky_quality_params': ['nrmse_clearsky']})
                    
                    clearsky_quality = xr.concat([nrmse_sym, nrmse_clearsky], 'clearsky_quality_params')
                    clearsky_quality = clearsky_quality.drop_vars('sun_params')
                    clearsky_quality = clearsky_quality.expand_dims({"channel":[clearsky_quality.channel]})
                    
                    dst[f'{param}_clearsky_quality'] = clearsky_quality
                    
                clearskylist.append(dst)
            
            clearsky = xr.concat(clearskylist, 'channel')
            clearsky.attrs['site'] = ds.site
            self._testresults = clearsky
        return self._testresults

    def add_clearsky_parameters2database(self, database):
        clearsky = self.test_results
        site = clearsky.site
        dt = pd.Timestamp.now()
        dateoffit = f'{dt.year:04d}-{dt.month:02d}-{dt.day:02d}'
            
        for channel in clearsky.channel:
            channel = int(channel)
            csel = clearsky.sel(channel=channel)

            for param in ['global_horizontal', 'diffuse_horizontal', 'direct_normal']:
                df = csel[f'{param}_symmetric_fit_params'].dropna('polynom_deg').to_dataframe().drop('channel', axis = 1).transpose()
                df.index = [0]
                df.columns.name = None
                df.columns = [f'psym{c}' for c in df.columns]
                psym = df
                
                # format fit rest
                df = csel[f'{param}_asymmetry_fit_params'].dropna('polynom_deg').to_dataframe().drop('channel', axis = 1).transpose()
                df.index = [0]
                df.columns.name = None
                df.columns = [f'prest{c}' for c in df.columns]
                prest = df
                
                # format quality
                df = csel[f'{param}_clearsky_quality'].to_dataframe().drop('channel', axis = 1).transpose()
                df.index = [0]
                df.columns.name = None
                cs_quali = df
                
                df = pd.DataFrame({'file': self.parent.path2file,
                                   'site': site,
                                   'channel': channel,
                                   'dateoffit': dateoffit,
                                   'observation': param,
                                  }, 
                                  index = [0])
                
                df = pd.concat([df,psym, prest, cs_quali], axis = 1)

                database.add_clearsky_params(df)
                
    
class SolarIrradiation(object):
    def __init__(self, dataset, site = None):
        self.dataset = self.unify_variables(dataset)
        self.clearsky = ClearSky(self)
        
        if isinstance(site, type(None)):
            assert('site' in dataset.attrs.keys()), 'If site is None, then the dataset has to have lat,lon,site, site_name, attributes'
            self.site = atmms.Station(lat= dataset.attrs['site_latitude'], 
                                      lon = dataset.attrs['site_longitude'], 
                                      alt = dataset.attrs['site_elevation'], 
                                      name = dataset.attrs['site_name'], 
                                      abbreviation = dataset.attrs['site'],)
        else:
            self.site = site
        
        
        self._sun_position = None
        self._solarspectrum = None
        self.path2solar_spectrum = '/User/htelg/fundgrube/reference_data/solar/spectrum/solar_spectral_irradiance_e490_00a_amo.nc'

        
    @property
    def toa_solarspectrum_1AU(self):
        """Solar spectrum at 1 AU"""
        if isinstance(self._solarspectrum, type(None)):
            self._solarspectrum = xr.open_dataset(self.path2solar_spectrum)
        return self._solarspectrum
    
    @property
    def toa_spectral_irradiance(self):
        """Top of the atmosphere spectral irradiances at channels in the dataset. 
        These are the actuall irradiances with sun earth distanced inbedded."""
        if "toa_spectral_irradiance" not in self.dataset:
            ds_solar = self.toa_solarspectrum_1AU# xr.open_dataset(self.path2solar_spectrum)
            rad_top = ds_solar.interp(wavelength=self.dataset.channel_wavelength, method='linear')# 
            self.tp_rad_topI = rad_top.copy()
            rad_top = rad_top / self.sun_position.sun_earth_distance.to_xarray()**2
            rad_top = rad_top.drop_vars('wavelength')
            # self.tp_rad_topII = rad_top.copy()
            self.dataset['toa_spectral_irradiance'] = rad_top.irradiance
        return self.dataset.toa_spectral_irradiance
            

        

    def unify_variables(self, dataset):
        """Seach for variable names containing global, diffuse and direct and 
        renames to global_horizontal, diffuse_horizontal, direct_normal"""
        #### variable cleaining
        # The exact variable name is sometimes
        ds = dataset.copy()
        for altvar in ['global_horizontal', 'diffuse_horizontal', 'direct_normal']:
            # altshort = altvar.split('_')[0] # This lead to problems when direct_horizontal and direct_normal are present
            match = [var for var in ds.variables if altvar in var]
            assert(len(match) < 2), f'There are multiple variables with {altvar} in it ({match}).'
            if len(match) == 0:
                # e.g. MFR has no direct
                continue
            ds = ds.rename({match[0]: altvar})
        return ds
        

    def register_file_in_database(self, database, overwrite = False):
        ds = self.dataset
        if 'creation_timestamp' in ds.attrs:
            date_created = ds.date_created
        else:
            date_created = 'NA'
            
        database.add_file(self.path2file, 
                             date=ds.datetime.values[0], 
                             date_created=date_created, 
                             site=ds.site, 
                             product_name = ds.product_name, 
                             product_version = ds.product_version, 
                             overwrite = overwrite)
    @property
    def sun_position(self):
        if isinstance(self._sun_position, type(None)):
            self._sun_position = self.site.get_sun_position(self.dataset.datetime)
        return self._sun_position
    
    def apply_calibration_langley(self, langley_calibration, calibrate_to = 'irradiance'):
        """Applies the langley calibration to the dataset. Make sure all other calibrations (spectral, head, cosine, etc.) 
        have been applied before. Note, the calibration to irradiance is still quite crued, as only the center wavelength 
        is considered instead of the full filter function! Make sure to keep the spectral irradiance at toa if you want to
        do further processing to transmission, OD, AOD ...
        
        Parameters
        ----------
        langley_calibration: 
            Instance of atmPy.radiation.retrievals.langley_calibration.Langley_Timeseries.
            This is the result of a langley calibration.
        calibrate_to: str
            Either 'irradiance' or 'transmission'. If 'irradiance' the output will be in W/m^2/nm.
            If 'transmission' the output will be unitless (0-1). See note above about top of the atmosphere irradiance.
        """

        assert(calibrate_to in ['irradiance', 'transmission'])
        # assert(isinstance(langley_calibration, atmlangcalib.Langley_Timeseries)), 
        assert(langley_calibration.__class__.__name__ == "Langley_Timeseries"), f'currently only the following instances are excepted: atmPy.radiation.retrievals.langley_calibration.Langley_Timeseries. I got {langley_calibration}'
        assert('calibrated_langley' not in self.dataset.attrs), 'it seems likt langley calibrations have already been applied.'
        lt = langley_calibration
        v0 = lt.V0_simple.V0 / self.sun_position.sun_earth_distance.to_xarray()**2
        toasi = self.toa_spectral_irradiance # this is already adjusted to sun earth distan
        v0 = v0.rename({'wavelength':'channel'})
        
        dsnew = self.dataset.copy()
        gh = self.dataset.global_horizontal / v0 # ratio of irradiance of what we had at the top of the atmosphere at normal
        if calibrate_to == 'irradiance':
            gh = gh * toasi # actual irradiance observed
        dsnew['global_horizontal'] = gh
        dsnew.global_horizontal.attrs['unit'] = 'W/m^2/nm'
        
        dh = self.dataset.diffuse_horizontal / v0 
        dh = dh * toasi
        dsnew['diffuse_horizontal'] = dh
        dsnew.diffuse_horizontal.attrs['units'] = 'W/m^2/nm'   
        
        dn = self.dataset.direct_normal / v0 
        dn = dn * toasi
        dsnew['direct_normal'] = dn
        dsnew.direct_normal.attrs['units'] = 'W/m^2/nm'
        

        dsnew.attrs['calibrated_langley'] = 'True'

        return self.__class__(dsnew)

    def apply_calibration_spectral(self, calibration):
        """
        This will assign the nominal channel wavelength and provide the exact
        channel central wavelengths.

        Parameters
        ----------
        calibration : xarray.Dataset
            Use atmPy.data_archives.NOAA_ESRL_GMD_GRAD.cal_facility.lab.read_mfrsr_cal 
            to open calibration file and use return as input here.

        Returns
        -------
        TYPE
            DESCRIPTION.
        TYPE
            DESCRIPTION.

        """
        assert(isinstance(calibration, xr.Dataset))
        assert('statistics' in calibration.variables), "I don't think the calibration file is a spectral calibration file for an MFR(SR). 'statistics' varible is missing"
        
        ds = self.dataset.copy()         
        ds['channel'] = calibration.channel
        ds['channel_wavelength'] = calibration.statistics.sel(stats = 'CENT', drop=True).astype(float)
        ds.attrs['calibrated_spectral'] = 'True'
        return self.__class__(ds) #returns the same class, allows for application to all subclasses
    
    def apply_calibration_responsivity(self, calibration, 
                                       varname_responsivity_spectral = 'responsivity_spectral',
                                       varname_dark_signal_spectral = 'dark_signal_spectral',
                                       ignore_has_been_applied_error = False):
        """
        This will calibrate for amplifier responsivity. This is sometimes
        applied multiple times, e.g. in MFR-type instruments where we have head
        and datalogger sensitivity. If executed multiple time the "has been 
        applied error" will trigger. Make sure to set ignore_has_been_applied_error
        to True

        Parameters
        ----------
        calibration : TYPE
            DESCRIPTION.
         : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.
        TYPE
            DESCRIPTION.

        """
        calibration = calibration.copy()
        calibration = calibration.rename({varname_responsivity_spectral:'responsivity_spectral',
                                          varname_dark_signal_spectral: 'dark_signal_spectral'})
        
        if not ignore_has_been_applied_error:
            if 'calibration_dark_signal' in self.dataset.attrs:
                assert(self.dataset.attrs['calibration_dark_signal'] != 'True'), 'Responds calibration already applied'        
            if 'calibration_responds' in self.dataset.attrs:
                assert(self.dataset.attrs['calibration_responds'] != 'True'), 'Responds calibration already applied'

        assert(isinstance(calibration, xr.Dataset))
        assert('dark_signal_spectral' in calibration.variables), "I don't think the calibration file is a spectral calibration file for an MFR(SR). 'dark_signal_spectral' varible is missing"
        
        ds = self.dataset.copy() 
        #### global horizontal
        
        #### dark signal
        da = ds.global_horizontal - calibration.dark_signal_spectral 
        #### responsivity
        da = da / calibration.responsivity_spectral
        
        da.attrs['unit'] = 'W * m^-2 * nm'
        da.attrs['calibration_responds'] = 'True'
        da.attrs['calibration_dark_signal'] = 'True'
        
        ds['global_horizontal'] = da
        
        #### diffuse horizontal
        if 'diffuse_horizontal' in ds.variables:
            da = ds.diffuse_horizontal - calibration.dark_signal_spectral
            da = da / calibration.responsivity_spectral
            
            da.attrs['unit'] = 'W * m^-2 * nm'
            da.attrs['calibration_dark_signal'] = 'True'
            da.attrs['calibration_responds'] = 'True'
            
            ds['diffuse_horizontal'] = da
            
        ds.attrs['calibration_dark_signal'] = 'True'
        ds.attrs['calibration_responds'] = 'True'
        return self.__class__(ds) #returns the same class, allows for application to all subclasses
    
    def apply_calibration_cosine(self, calibration):
        
        if 'clalibration_cosine' in self.dataset.attrs:
            assert(self.dataset.attrs['clalibration_cosine'] != 'True'), 'Responds calibration already applied'  
            
        ds = self.dataset.copy()
        
        #### for diffuse or global in case of an MFR
        cal_angle = 45
        ew = calibration.spectral_EW.interp(Angle = [cal_angle, -cal_angle]).sum(dim = 'Angle') / 2 
        ns = calibration.spectral_NS.interp(Angle = [cal_angle, -cal_angle]).sum(dim = 'Angle') / 2 
        cal = (ew + ns) / 2
        
        # The following should only happen when no global and direct is measured, ideally only for upwelling
        #### MFR
        if not 'diffuse_horizontal' in ds.variables: 
            ds['global_horizontal'] = ds.global_horizontal / cal
        
        # only if the direct component is resolved the following is relevant
        #### MFRSR
        else:
            assert('diffuse_horizontal' in ds.variables)
            
            #### - diffuse
            ds['diffuse_horizontal'] = ds.diffuse_horizontal / cal
            
            #### - direct
            sp = atmsol.SolarPosition(self.sun_position.azimuth, np.pi/2 - self.sun_position.elevation, unit = 'rad')
            
            # NS
            # interpolate the cosine respond with the particular angles resulting from the projetion
            # This results in :
            #     * calibration value as a function of time
            
            da = calibration.spectral_NS.interp(Angle = np.rad2deg(sp.projectionNS_angle) - 90)
            cos_cal_NS = da.rename({'Angle': 'datetime'}).assign_coords(datetime = ('datetime', self.sun_position.index))
            
            
            # The calibration value needs to be normalized with the relevant component of the solar radiation
            # With other words how much light is actually comming this way?
            cos_cal_NS_norm = cos_cal_NS * xr.DataArray(sp.projectionNS_norm)
            
            # Do the same for EW
            
            da = calibration.spectral_EW.interp(Angle = np.rad2deg(sp.projectionEW_angle) - 90)
            cos_cal_EW = da.rename({'Angle': 'datetime'}).assign_coords(datetime = ('datetime', self.sun_position.index))
            cos_cal_EW_norm = cos_cal_EW * xr.DataArray(sp.projectionEW_norm)
            
            # Sum NS and EW
            cos_cal_sum = cos_cal_EW_norm + cos_cal_NS_norm
            
            # Divide by the sum of **norms** (Not the calibration value! As we are dealing with vectors the sum is not automatically num
            sumofnorm = sp.projectionEW_norm + sp.projectionNS_norm
            cos_cal_sum_nom = cos_cal_sum / xr.DataArray(sumofnorm)
            
            # apply final cosine correction to the data
            
            
            ds['direct_horizontal'] = (ds.global_horizontal - ds.diffuse_horizontal)
            ds['direct_horizontal'] = ds.direct_horizontal / cos_cal_sum_nom
            ds['direct_normal'] = ds.direct_horizontal / xr.DataArray(np.sin(self.sun_position.elevation))
        

        
            #### - compose global based on cosine corrected direct and diffuse
            ds['global_horizontal'] = ds.direct_horizontal + ds.diffuse_horizontal
            
        
        ds.attrs['clalibration_cosine'] = 'True'
        out = self.__class__(ds) #returns the same class, allows for application to all subclasses
        if 'diffuse_horizontal' in ds.variables:
            out.tp_cos_cal_sum = cos_cal_sum
            out.tp_cos_cal_EW_norm = cos_cal_EW_norm
            out.tp_cos_cal_NS_norm = cos_cal_NS_norm
            out.tp_cos_cal_sum_nom = cos_cal_sum_nom
        return out
        

class GlobalHorizontalIrradiation(SolarIrradiation):
    def __init__(self, dataset):
        super().__init__(dataset)

class DiffuseHorizontalIrradiation(SolarIrradiation):
    def __init__(self, dataset):
        super().__init__(dataset)

class DirectNormalIrradiation(SolarIrradiation):
    def __init__(self, dataset, 
                 site = None, 
                 langley_fit_settings = None,
                 calibration_strategy = 'johns',
                 metdata = 'surfrad',
                #  raise_error_if_no_channel_wavelength = True, #this is for the actual wavelength, the calibrated ones not the nominal ones in the channel dimension
                 settings_ozone = None):
        
        """
        Parameters
        ----------
        dataset : xarray.Dataset
            The data.
        site : atmPy.general.measurement_site.Station, optional
            Site information should be in the dataset, if it is not or you want to overwite it (doulecheck if it actually does overwrite), set it here (forgot which format). The default is None.
        langley_fit_settings : dict, optional
            DESCRIPTION. The default is None.
        calibration_strategy : str, optional
            Depending on the strategie different files will be needed. The default is 'johns'.
                johns: uses the langley calibration files from John.
                sp02: uses the sp02 calibration files
                toa_radiation: This assumes calibrated spectral radiation in the dataset. It will then get transmission using the toa radiation
        metdata : str, optional
            Which metdata (or data format) to use?. The default is 'surfrad'.
                surfrad: (double check this) uses the surfrad met data (netcdf version)
        settings_ozone : various, optional
            A way to set the ozone. Currently most calibration strategies will pick there own way to get ozone (todo: this needs to change!). Eventually this should take a function that can read a file and aligns it to the dataset. Currently, for testing, only a single number is possible ... probly 300 DU.

        """

        super().__init__(dataset, site = site)
        self.raw_data = dataset #this is not exactly raw, it is cosine corrected voltag readings, thus, un-calibrated irradiances
        self.variable_name_alternatives_direct_normal = ['direct_normal', 'direct_normal_irradiation', 'direct_normal_irradiance']
        var_ops_direct_normal = [v for v in self.variable_name_alternatives_direct_normal if v in self.raw_data.data_vars]
        if len(var_ops_direct_normal)!= 1:
            if len(var_ops_direct_normal) == 0:
                raise ValueError('Could not find direct normal variable in dataset')
            if len(var_ops_direct_normal) > 1:
                raise ValueError('Multiple potential direct normal variables found in dataset: {}'.format(var_ops_direct_normal))
            else:
                raise RuntimeError('This should never happen (id2334442313)')
        self.variable_name_direct_normal = var_ops_direct_normal[0]

        # if raise_error_if_no_channel_wavelength:
        self.variable_name_alternatives_channel_wavelength = ['channel_wavelength', 'channel_center'] 
        variable_name_channel_wavelength = [v for v in self.variable_name_alternatives_channel_wavelength if v in self.raw_data.data_vars]
        if len(variable_name_channel_wavelength)!= 1:
            if len(var_channel_center) == 0:
                raise ValueError('Could not find channel center variable in dataset. This might mean that the data has not been spectrally calibrated yet. Use the apply_calibration_spectral method first.')
            if len(var_channel_center) > 1:
                raise ValueError('Multiple potential channel center variables found in dataset: {}'.format(var_channel_center))
            else:
                raise RuntimeError('This should never happen (id2334442313)')
        self.variable_name_channel_wavelength = variable_name_channel_wavelength[0]

        self.langley_fit_settings = langley_fit_settings
        self.settings_langley_airmass_limits = (2.5,4)
        self.settings_calibration = calibration_strategy 
        self.settings_metdata = metdata
        self.settings_ozone = settings_ozone
        self.mfrsr_history = '/home/grad/htelg/projects/AOD_redesign/MFRSR_History.xlsx'
        self.precipitable_water_varname = None
        self.settings_rayleigh = 'firstprinciple' #options: 'john', 'firstprinciple'. "john" might have an errror in it, not sure if the error happend appon translation or his code actually has that error. results in an over estimation of the rayleigh od.
        self.path2absorption_correction_ceoff_1625 = '1625nm_absorption_correction_coefficience.nc'
        self._am = None
        self._pm = None
        # self._transmission = None
        # self._od_rayleigh = None
        self._od_co2ch4h2o = None
        # self._tpw = None
        self._pwv_lut = None
        # self._aod = None
        # self._metdata = None
        # self._od_ozone = None
    
    @property
    def met_data(self):
        # if isinstance(self._metdata, type(None)):
        if 'pressure' not in self.raw_data:
            try:
                if isinstance(self.settings_metdata, types.FunctionType):
                    self._metdata = self.settings_metdata()
                else:
                    if self.settings_metdata == 'surfrad':
                        self._metdata = self._get_surfrad_met_data()
                    elif self.settings_metdata == 'baseline':
                        self._metdata = self._get_baseline_met_data()
                    else:
                        assert(False), 'moeeeep!'
            except FileNotFoundError as e:
                # print('Could not find met data. Make sure to set this value and or that the "settings_metdata" attribute is set correctly! The following error was raised:')
                raise FileNotFoundError(f'Could not find met data. Make sure the "settings_metdata" attribute is set correctly! The following error was raised: {e}') from e
        return self.raw_data[['pressure','temperature']]
    
    @met_data.setter
    def met_data(self, value):
        if isinstance(value, xr.Dataset):
            dsmet = value
        elif pl.Path(value).is_file():
            dsmet = xr.open_dataset(value)
        else:
            ValueError(f'Expected xr.Dataset, or valid Path (as str or pathlib.Path). Got: {value}')
        dsmet_clean = xr.Dataset()
        if 'datetime' not in dsmet.coords:
            altcoords = ['time',]
            found = [coo for coo in dsmet.coords if coo in altcoords]    
            assert(len(found) == 1), f'Problem with time coordinate. Expeced one of: {altcoords + ['datetime']}. Available coordinates: {list(dsmet.coords)}'
            dsmet = dsmet.rename_dims({found[0]:'datetime'})
            dsmet = dsmet.rename_vars({found[0]:'datetime'})

        if fill_value := dsmet.attrs.get('fill_value', False):
            dsmet = dsmet.where(dsmet != fill_value)

        if 'pressure' not in dsmet:
            altcoords = []
            found = [coo for coo in dsmet if coo in altcoords]    
            assert(len(found) == 1), f'Problem with pressure variable. Expeced one of: {altcoords + ['pressure']}. Available coordinates: {list(dsmet.coords)}'
            dsmet = dsmet.rename_vars({found[0]:'pressure'})
        dsmet_clean['pressure'] = dsmet['pressure']

        if 'temperature' not in dsmet:
            altcoords = ['air_temperature', 'temp']
            found = [coo for coo in dsmet if coo in altcoords]    
            assert(len(found) == 1), f'Problem with temperature variable. Expeced one of: {altcoords + ['pressure']}. Available coordinates: {list(dsmet.coords)}'
            dsmet = dsmet.rename_vars({found[0]:'temperature'})
        dsmet_clean['temperature'] = dsmet['temperature']
      
        dsmet_clean = dsmet_clean.interp({'datetime':self.raw_data.datetime}, method = 'linear')
        self.tp_dsmet_clean = dsmet_clean
        # dst = dsmet[['pressure','temperature']].interp({'datetime':self.raw_data.datetime})
        # dst = dst.where(dst > -90)

        self.raw_data[['pressure','temperature']] = dsmet_clean[['pressure','temperature']]
        return
    
    def _get_baseline_met_data(self):
        site = self.raw_data.attrs['site']
        dtf = pd.to_datetime(self.raw_data.datetime.values[0])
        p2metfld_base = '/nfs/iftp/aftp/g-rad/baseline/'
        p2metfld_base = pl.Path(p2metfld_base)
        p2f = p2metfld_base.joinpath(f'{site}/{dtf.year}/{site}{str(dtf.year)[2:]}{dtf.day_of_year:03d}.dat')
        ds = atmbl.read(p2f)
        pt_interp = ds[['pressure','temp']].interp({'datetime':self.raw_data.datetime})
        return pt_interp
    
    def _get_surfrad_met_data(self):
        """
        Loads the surfrad data (from my netcdf version), interpolates and adds 
        to the dataset. Location is currently hardcoded, this will likey cause problems sooner or later

        Returns
        -------
        None.

        """
        start_dt = pd.to_datetime(self.raw_data[self.variable_name_direct_normal].datetime.min().values)
        end_dt = pd.to_datetime(self.raw_data[self.variable_name_direct_normal].datetime.max().values)
        
        # open relevant files
        site = self.site.abb
        fns = []
        for dt in [start_dt, end_dt]:
            fns.append(f'/nfs/grad/surfrad/products_level1/radiation_netcdf/{site}/srf_rad_full_{site}_{dt.year:04d}{dt.month:02d}{dt.day:02d}.nc')
        
        ds = xr.open_mfdataset(np.unique(fns)) # unique in case start and end is the same
        
        # interpolate to mfrsr datetime index
        pt_interp = ds[['pressure','temp']].interp({'datetime':self.raw_data.datetime})
        # self._metdata = pt_interp
        # add to dataset
        # for var in pt_interp:
        #     self.raw_data[var] = pt_interp[var]
        return pt_interp
            
            
    def _apply_calibration_sp02(self):
        """
        Loads the calibration files and applies by interpolation.

        Returns
        -------
        None.

        """
        calibrations = atmlangcalib.load_calibration_history()
        cal = calibrations[int(self.raw_data.attrs['serial_no'])]
        self.calibration_inst = cal
        # use the mean and only the actual channels, other channels are artefacts
        cal = cal.results['mean'].loc[:,self.raw_data.channel_center.values].sort_index()
    
        #### interpolate and resample calibration (V0s)
        dt = self.raw_data.datetime.to_pandas()
        calib_interp = pd.concat([cal,dt]).drop([0], axis = 1).sort_index().interpolate().reindex(dt.index)
        #### correct VOs for earth sun distance see functions above
        calib_interp_secorr = calib_interp.divide(self.sun_position.sun_earth_distance**2, axis = 0)
        self.calibration_data = calib_interp_secorr
        
        #### match channels for operation
        channels = self.raw_data.channel_center.to_pandas()
        raw_data = self.raw_data.direct_normal_irradiation.to_pandas().rename(columns = channels)
        raw_data.columns.name = 'wl'
        
        #### get transmission
        transmission = raw_data/calib_interp_secorr
        transmission = xr.DataArray(transmission)
        transmission = transmission.rename({'wl':'channel'})
        self.calibration = 'sp02'
        self._transmission =transmission
        
        
    def _apply_calibration_atm_gam(self,p2fld='/export/htelg/data/grad/surfrad/mfrsr/langleys_concat/tbl/',
                                   th=0.02,
                                    order_stderr=2,
                                    lam_overtime=2.5e4,
                                    ns_overtime=2,
                                    lam_season=1e4,
                                    ns_season=6,
                                    lands = None):
        if isinstance(lands, type(None)):
            print('langleys.load.', end = '', flush = True)
            lands = atmlangcalib.open_langley_dailys(#end = '20211001',
                                                  p2fld=p2fld,)

            print('gampredict.', end = '')
            lands.predict_until = pd.Timestamp.now()
            lands._get_v0_gam(th=th,
                      order_stderr=order_stderr,
                      lam_overtime=lam_overtime,
                      ns_overtime=ns_overtime,
                      lam_season=lam_season,
                      ns_season=ns_season,
                    )
            # Days after this date will be considered as temporary. These days will be recalculated every day until this date is no longer before the date under question
            # the parameters where selected based on a rough guess, no validation or evaluation was conducted
            th_predict = 0.003
            req_no_good_langleys = 15
            wl = 500
            istderr = lands.dataset.sel(fit_results = 'intercept_stderr', wavelength = wl).langley_fitres
            date_predict = pd.to_datetime(lands.dataset.sel(fit_results = 'intercept', wavelength = wl).langley_fitres.where(istderr < th_predict).dropna('datetime').datetime[-req_no_good_langleys].values)
            lands.date_predict = date_predict
            print('done')
        else:
            pass
        
        self.lands = lands
        v0 = lands.v0prediction.dataset
        v0 = v0.rename({'wavelength':'channel'})
        
        #### FIXME: the below should no longer be required
        # v0 = v0.assign_coords({'channel': [ 415,  500, 1625,  670,  870,  940]})
        
        v0_interp = v0.interp(datetime = self.raw_data.direct_normal_irradiation.datetime)
        sedistcorr = self.sun_position.sun_earth_distance**2
        v0_interp_secorr = v0_interp.V0 / sedistcorr.to_xarray()
        self.tp_v0_interp_secorr = v0_interp_secorr
        transmission = self.raw_data.direct_normal_irradiation / v0_interp_secorr
        self.calibration = 'atm_gam'
        self._transmission =transmission

        return
        
    def _apply_calibration_singleV0(self):
        # 
        pass
        

    def _apply_calibration_johns(self):
        def datetime2angulardate(dt):
            if dt.is_leap_year:
                noofday = 366
            else:
                noofday = 365
            fract_year = dt.day_of_year / noofday
            yyyy_frac = dt.year+fract_year
            angular_date = yyyy_frac * 2 * np.pi
            return angular_date
        
        #### open file that contains the calibration parameters (the fit parameters from john)
        site = self.site.abb
        p2f=f'/home/grad/surfrad/aod/{site}_mfrhead'
        if not pl.Path(p2f).is_file():
            raise FileNotFoundError(f'Calibration parameter file {p2f} not found. This might be due to the wrong value in settings_calibration. Current value is {self.settings_calibration}.')
        langley_params = atmlangcalib.read_langley_params(p2f = p2f)
        
        #### get V0 for that date
        # get closesed calibration parameter set based on first timestamp. 
        # This was the wrong approach that I assumed 
        # falsely. The right way is to use the coefficients from before that date
        
        #old:
        # dt = pd.to_datetime(self.raw_data.datetime.values[0])
        # assert((langley_params.datetime.to_pandas()-dt).abs().min() < pd.to_timedelta(1.5*365, 'days')), f'The closest fit the Langleys is more than 1 year out ({(langley_params.datetime.to_pandas()-dt).abs().min()}). This probably means that John needs update fit parameter file.'
        # idxmin = (langley_params.datetime.to_pandas()-dt).abs().argmin()
        # langley_params_sel = langley_params.isel(datetime = idxmin)
        
        #new:using the coefficients from the last datetime before the date of this .ccc file.
        dt = pd.to_datetime(self.raw_data.datetime.values[0])
        
        dtmin = (dt - langley_params.datetime.to_pandas()) / pd.to_timedelta(1, 'day')
        dtmin[dtmin < 0] = np.nan
        idxmin = dtmin.argmin()
        langley_params_sel = langley_params.isel(datetime = idxmin)
        
        
        # turn date into that johns angular_date (see definition of get_izero in aod_analysis.f 
        angular_date = datetime2angulardate(dt)
        
        # apply the parameters to get the V0 and V0_err
        # V0 =
        cons_term = langley_params_sel.V0.sel(V0_params = 'const')
        lin_term = langley_params_sel.V0.sel(V0_params = 'lin') * angular_date 
        sin_term = langley_params_sel.V0.sel(V0_params = 'sin') * np.sin(angular_date)
        cos_term = langley_params_sel.V0.sel(V0_params = 'cos') * np.cos(angular_date)
        V0 = cons_term + lin_term + sin_term + cos_term
        
        V0df = pd.DataFrame(V0.to_pandas())
        
        #### correct for earth sun distance
        sedistcorr = pd.DataFrame(self.sun_position.sun_earth_distance**2)
        sedistcorr.columns = [0]
        
        # self.tp_V0df = V0df
        # self.tp_sedistcorr = sedistcorr
        calib_interp_secorr = V0df.dot(1/sedistcorr.transpose()).transpose()
        calib_interp_secorr.rename({936:940}, axis = 1, inplace=True)
        calib_interp_secorr.columns.name = 'channel'
        # self.tp_calib_interp_secorr = calib_interp_secorr
        self.tp_calib_interp_secorr = calib_interp_secorr
        raw_data = self.raw_data.direct_normal_irradiation.to_pandas()
        transmission = raw_data/calib_interp_secorr
        transmission = xr.DataArray(transmission)     
        self.calibration = 'johns'
        self._transmission =transmission

    @property
    def od_rayleigh(self):
        if 'od_rayleigh' not in self.raw_data:
            if self.settings_rayleigh == 'john':
                odr = xr.concat([atmraylab.rayleigh_od_johnsmethod(self.met_data.pressure, chan.values) for chan in self.raw_data[self.variable_name_channel_wavelength]], 'channel')
            elif self.settings_rayleigh == 'firstprinciple':
                odr = xr.concat([atmraylab.rayleigh_od_first_principles(chan.values, pressure = self.met_data.pressure) for chan in self.raw_data[self.variable_name_channel_wavelength]], 'channel')
            self._od_rayleigh = odr
            self.tp_odr = odr
            self.raw_data['od_rayleigh'] = odr
        return self.raw_data['od_rayleigh']


    @property
    def ozone_absorption_spectrum(self):
        if 'ozone_absorption_spectrum' not in self.raw_data:
           #try the default folder
           self.ozone_absorption_spectrum = '/home/grad/surfrad/aod/ozone.coefs'
            # else:
            #     assert(False), 'set ozone_absorption_spectrum!'
        return self.raw_data['ozone_absorption_spectrum']
    

    @ozone_absorption_spectrum.setter
    def ozone_absorption_spectrum(self, value):
        if pl.Path(value).is_file():
            ozone_abs_coef = pd.read_csv(value, index_col= 0, sep = ' ', names=['wavelength', 'coeff'])
            ozone_abs_coef = ozone_abs_coef.to_xarray()
            ozone_abs_coef = ozone_abs_coef.coeff
            
        elif isinstance(value, xr.Dataset):
            ozone_abs_coef = value
        else:
            assert(False), 'no idea what to do'
        self.tp_ozone_absorption_spectrum = ozone_abs_coef
        self.raw_data['ozone_absorption_spectrum'] = ozone_abs_coef
        return
  
    @property
    def ozone_absorption_by_channel(self):
        if 'ozon_absoption_by_channel' not in self.raw_data:
            ozon2bychannel = self.ozone_absorption_spectrum.interp(wavelength = self.raw_data[self.variable_name_channel_wavelength])
            self.tp_ozon2bychannel = ozon2bychannel.copy()
            ozon2bychannel = ozon2bychannel.drop('wavelength')

            # 368 and 1050 are outside the ozon spectrum. Values should be very small
            ozon2bychannel = ozon2bychannel.fillna(0)

            self.raw_data['ozon_absoption_by_channel'] = ozon2bychannel
        return self.raw_data['ozon_absoption_by_channel']

    @property
    def ozone_data(self) -> xr.DataArray:
        """
        Returns
        -------
        xr.DataArray
        Returns the ozone data as an xarray.DataArray. 
        
        Setter
        ------
        Can be set with: xr.Dataset, xr.DataArray, int, float, str, pathlib.Path
            If str is given it can be:
                - a stragy, e.g. 'johns' or 'sp02'
                - a path to a file (e.g. netcdf or csv)
            if int or float is given it is assumed to be a constant value (in DU) for all times.      
        """
        if 'ozone_data' not in self.raw_data:
            assert(False), 'Set ozone_data'
        return self.raw_data['ozone_data']
    
    @ozone_data.setter
    def ozone_data(self, value: typing.Union[xr.Dataset, xr.DataArray, int, float, str, pl.Path]) -> None:
        if isinstance(value, (xr.Dataset, xr.DataArray)):
            ozone_data = value
            print('whatup')
            self.tp_ozon_data1 = ozone_data.copy()
        elif isinstance(value, str):
            if pl.Path(value).is_file():
                ozone_data = xr.open_dataset(value)
            elif value == 'johns':
                p2f = f'/home/grad/surfrad/aod/ozone/{self.site.abb}_ozone.dat'
                ozone = atmsrf.read_ozon(p2f)
                dt = pd.to_datetime(pd.to_datetime(self.raw_data.datetime.values[0]).date())
                ozone_data = float(ozone.interp(datetime = dt).ozone)
            elif value == 'sp02':
                assert(self.site.abb == 'brw'), 'only works for barrow so far'

                #### read ozon and interpolate
                dt = pd.to_datetime(self.raw_data.datetime.values[0])
                p2fld= pl.Path('/nfs/stu3data2/Model_data/merra_2/barrow/merra_M2I1nxasm_5.12.4/')
                
                p2f_list = []
                # load day before and after for interolation purposes
                for i in list(range(-1,2)):
                    dtt = dt + pd.to_timedelta(i, 'd')
                
                    pattern = f'MERRA2_*.inst1_2d_asm_Nx.{dtt.year}{dtt.month:02d}{dtt.day:02d}.nc'
                
                    res = list(p2fld.glob(pattern))
                
                    if len(res) == 0:
                        assert(False), 'not matching ozon file'
                
                    assert(len(res) ==1), 'There is more than 1 matching ozon file ... inconceivable!'
                
                    p2f_list.append(res[0])
                
                ozone = xr.open_mfdataset(p2f_list)
                ozone = ozone.rename({'time': 'datetime'})
                
                ozone_data = ozone.interp(datetime = self.raw_data.datetime).TO3

        elif isinstance(value, (int, float)):
            ozone_data = xr.DataArray(value, dims=('datetime',), coords={'datetime': self.raw_data.datetime})

        else:
            assert(False), 'no idea what to do with this~!!?!?'

        if isinstance(ozone_data, xr.Dataset):
            varname_options = ['ozone_data', 'ozone']
            varname_found = [v for v in varname_options if v in ozone_data]
            assert(len(varname_found) == 1), f'Could not find ozone variable in dataset. Expected one of: {varname_options}. Available variables: {list(value.data_vars)}'
            ozone_data = ozone_data[varname_found[0]]

        if 'datetime' not in ozone_data.coords:
            if self.verbose:
                print('Renaming time coordinate to datetime for ozone data.')
            time_coord_options = ['time',]
            found = [coo for coo in ozone_data.coords if coo in time_coord_options]
            assert(len(found) == 1), f'Problem with time coordinate. Expeced one of: {time_coord_options + ['datetime']}. Available coordinates: {list(ozone_data.coords)}'
            ozone_data = ozone_data.rename({found[0]:'datetime'})

        # interp to match times
        ozone_data = ozone_data.interp(datetime = self.raw_data.datetime)

        if ozone_data.sum() == 0:
            raise ValueError('Ozone data is all none or zero! Try setting it with a constant value.')
        self.raw_data['ozone_data'] = ozone_data
        self.raw_data.ozone_data.attrs = dict(long_name = "Total column ozone",
                                             units     = "Dobson Units (DU)",)
        return

    @property
    def od_ozone(self):
        if 'od_ozone' not in self.raw_data:
        # if isinstance(self._od_ozone, type(None)):
            #### read spectral function of ozon absorption coeff 
            # p2f = '/home/grad/surfrad/aod/ozone.coefs'
            # ozone_abs_coef = pd.read_csv(p2f, index_col= 0, sep = ' ', names=['wavelength', 'coeff'])
            # ozone_abs_coef = ozone_abs_coef.to_xarray()
            # ozon2bychannel = ozone_abs_coef.interp(wavelength = self.raw_data[self.variable_name_channel_wavelength])
            # self.tp_ozon2bychannel = ozon2bychannel.copy()
            # ozon2bychannel = ozon2bychannel.drop('wavelength')

            # 368 and 1050 are outside the ozon spectrum. Values should be very small
            # ozon2bychannel.coeff[ozon2bychannel.coeff.isnull()] = 0
            
            # if self.settings_calibration in ['johns', 'atm_gam']:
                #### read ozon concentration file for particular site and extract data for this day
                # p2f = f'/home/grad/surfrad/aod/ozone/{self.site.abb}_ozone.dat'
                # ozone = atmsrf.read_ozon(p2f)
                # dt = pd.to_datetime(pd.to_datetime(self.raw_data.datetime.values[0]).date())
                # total_ozone = float(ozone.interp(datetime = dt).ozone)
                
            
            # elif self.settings_calibration == 'sp02':
                # assert(self.site.abb == 'brw'), 'only works for barrow so far'

                # #### read ozon and interpolate
                # dt = pd.to_datetime(self.raw_data.datetime.values[0])
                # p2fld= pl.Path('/nfs/stu3data2/Model_data/merra_2/barrow/merra_M2I1nxasm_5.12.4/')
                
                # p2f_list = []
                # # load day before and after for interolation purposes
                # for i in list(range(-1,2)):
                #     dtt = dt + pd.to_timedelta(i, 'd')
                
                #     pattern = f'MERRA2_*.inst1_2d_asm_Nx.{dtt.year}{dtt.month:02d}{dtt.day:02d}.nc'
                
                #     res = list(p2fld.glob(pattern))
                
                #     if len(res) == 0:
                #         assert(False), 'not matching ozon file'
                
                #     assert(len(res) ==1), 'There is more than 1 matching ozon file ... inconceivable!'
                
                #     p2f_list.append(res[0])
                
                # ozone = xr.open_mfdataset(p2f_list)
                # ozone = ozone.rename({'time': 'datetime'})
                
                # total_ozone = ozone.interp(datetime = self.raw_data.datetime).TO3
            
                
            # else:
            #     if isinstance(self.settings_ozone, (int, float)):
            #         total_ozone = self.settings_ozone
            #     else:
            #         assert(False), f'not a known calibration strategy: {self.settings_calibration}. And unknown ozone setting: {self.settings_ozone}. Make sure to set those.'
                
            #### scale abs. coef. to total ozone to get OD_ozone
            od_ozone = self.ozone_absorption_by_channel * (self.ozone_data/1000) #any idea why we devide by 1000?
            od_ozone = od_ozone.fillna(0)
            # self._od_ozone = od_ozone
            self.raw_data['od_ozone'] = od_ozone
        return self.raw_data['od_ozone']

        
    def _retreive_pwv_from_940_channel(self, lookuptable: xr.DataArray) -> xr.DataArray:
        """Retrieve precipitable water from 940nm channel using provided lookup table.
        """
        ds_pwdlut = lookuptable
        ammin = ds_pwdlut.air_mass.min()
        ammax = ds_pwdlut.air_mass.max()
        def row2pwv(row):
            if  ammin<= row.airmass <= ammax:
                # print('buba')
                aodt = self.aod.sel(channel = 940, datetime = row.name)
                # the lut is in optical depth along the path, so not normalize to air mass
                aodt *= row.airmass
                lut_at_am = ds_pwdlut.interp(air_mass = row.airmass, method='linear')
                lut_at_am_inv = xr.DataArray(lut_at_am.pwv.values, coords = {'od': lut_at_am.optical_depth.values})
                pwv = float(lut_at_am_inv.interp(od = aodt).values)
            else:
                pwv = np.nan
            return pwv
        df = self.sun_position.apply(row2pwv, axis = 1)
        da = df.to_xarray().rename('precipitable_water')
        da.attrs = dict(long_name = "Precipitable water",
                        units     = "cm",
                        comment   = "Equivalent liquid water depth; not CF-standard, convertible to 0.1 kg m-2 per mm.",
                        )
        return da

    @property
    def precipitable_water(self):
        """Set/Returns precipitable water in cm. Make sure to set it first using the setter!
        Returns
        -------
        xr.DataArray
        The precipitable water in cm.

        Setter
        ------
        value : xr.Dataset, str, int, float
            Either set the precipitable water directly as a constant, or dataarray (e.i. through xr.Dataset), or provide a lookup table that can be used to retreive the precipitable water from the 940nm channel. If a path is provided as str or pathlib.Path, the file will be opened as an xr.Dataset. 
            If xr.Dataset:
                - it should either contain a variable with precipitable water (e.g. from a sounding).
                - or a lookup table that can be used to retreive the precipitable water from the 940nm channel.
            If float or int: Value is used as constant precipitable water for all timestamps, unit is assumed to be cm.
            If str or pathlib.Path: Path to a netcdf file that contains either the precipitable water variable or a lookup table that can be used to retreive the precipitable water from the 940nm channel.
        -----------
        """

        
        if 'precipitable_water' not in self.raw_data:
            if isinstance(self._pwv_lut, xr.Dataset):
                if self.verbose: 
                    print('Retrieving precipitable water from 940nm channel using lookup table.')
                da =  self._retreive_pwv_from_940_channel(self._pwv_lut)
                self.raw_data['precipitable_water'] = da
            else:
                raise AttributeError('precipitable water is not set. Use the setter to set it first!')
            
        return self.raw_data['precipitable_water']

    # @property
    # def precipitable_water(self):
    #     if isinstance(self._tpw, type(None)):
    #         sitedict = {'Bondville': 'BND', 
    #                     'Fort Peck': 'FPK',
    #                     'Goodwin Creek': 'GWN', 
    #                     'Table Mountain': 'TBL',
    #                     'Desert Rock': 'DRA',
    #                     'Penn State': 'PSU', 
    #                     'ARM SGP': 'sgp', 
    #                     'Sioux Falls': 'SXF',
    #                     'Canaan Valley': 'CVA',
    #                    }
            
    #         #### get path to relevant soundings files
    #         files = pd.DataFrame()
    #         datetime = pd.to_datetime(self.raw_data.datetime.values[0])
    #         for dt in [datetime, datetime + pd.to_timedelta(1, 'day')]:
    #             # pass
            
    #             p2fld = pl.Path(f'/nfs/grad/surfrad/sounding/{dt.year}/')
    #             searchstr = f'{dt.year}{dt.month:02d}{dt.day:02d}*.int'
    #             df = pd.DataFrame(p2fld.glob(searchstr), columns=['p2f'])
    #             files = pd.concat([files, df])
            
    #         files.sort_values('p2f', inplace = True)
            
            
    #         tpwts = []
    #         # tpsoundis = []
    #         for p2f in files.p2f:
    #             # pass
            
    #             soundi = atmsrf.read_sounding(p2f)
    #             # tpsoundis.append({"sounding" : soundi, "fn":p2f})
    #             tpwts.append(soundi.precipitable_water.expand_dims({'datetime': [soundi.data.attrs['datetime'],]}))
                
    #         # self.tp_soundi = tpsoundis
    #         tpw = xr.concat(tpwts, 'datetime')
    #         # self.tp_tpw_all_1 = tpw.copy()
    #         tpw = tpw.assign_coords(site = [sitedict[str(s.values)].lower() for s in tpw.site])
    #         # self.tp_tpw_all_2 = tpw.copy()
    #         tpw = tpw.sel(site = self.raw_data.site)
    #         tpw = tpw.drop(['site'])
            
    #         tpw = tpw.interp({'datetime': self.raw_data.datetime})
    #         self._tpw =tpw
    #     return self._tpw


    @precipitable_water.setter
    def precipitable_water(self, value):
        if isinstance(value, xr.Dataset):
            ds = value
        elif isinstance(value, str): 
            if pl.Path(value).is_file():
                ds = xr.open_dataset(value)
            else:
                raise FileNotFoundError(f'File {value} not found')
        elif isinstance(value, (int, float)): 
            ds = xr.Dataset()
            ds['precipitable_water'] =xr.DataArray(value, dims=('datetime',), coords={'datetime': self.raw_data.datetime})
            # return self.raw_data['precipitable_water']
        else:
            ValueError(f'Expected xr.Dataset, or valid Path (as str or pathlib.Path). Got: {value}')

        # ds_clean = xr.Dataset()
        if 'product' in ds.attrs:
            if ds.attrs["product"] == "mfrsr_pwv_lut":
                # this is a lookup table to retreive pwv from 940nm channel, we will process this when the getter is called.
                self._pwv_lut = ds
                return
            else:
                pass
        
        if 'datetime' not in ds.coords:
            altcoords = ['time',]
            found = [coo for coo in ds.coords if coo in altcoords]    
            assert(len(found) == 1), f'Problem with time coordinate. Expeced one of: {altcoords + ['datetime']}. Available coordinates: {list(dsmet.coords)}'
            ds = ds.rename_dims({found[0]:'datetime'})
            ds = ds.rename_vars({found[0]:'datetime'})

        if fill_value := ds.attrs.get('fill_value', False):
            ds = ds.where(ds != fill_value)

        if not isinstance(self.precipitable_water_varname, type(None)):
            ds_clean = ds[self.precipitable_water_varname]
        else:
            varlist = [var for var in ds if 'precipitable' in var]
            if len(varlist) > 0:
                assert(len(varlist) == 1), f'more then one variable found that indicate precipitable water ({varlist}). Try setting self.precipitable_water_varname with the appropriete variable name.'
                ds_clean = ds[varlist[0]]
        
        ds_clean = ds_clean.interp({'datetime':self.raw_data.datetime}, method = 'linear')
        # self.tp_dsmet_clean = ds_clean
        if 'units' in ds_clean.attrs:
            if ds_clean.attrs['units'] == 'mm':
                ds_clean = ds_clean / 10 # convert to cm
            elif ds_clean.attrs['units'] in ['cm', 'centimeter', 'centimeters']:
                pass
            else:
                assert(False), f'precipitable water unit {ds_clean.attrs["units"]} not recognized. Set the precipitable_water_varname attribute to the correct variable name and make sure it is in cm.'
        else:
            print('precipitable water variable has no unit attribute. Assuming it is in cm.')
            
        ds_clean.attrs['units'] = 'cm'
        ds_clean.attrs['long_name'] = 'Precipitable water'
        self.raw_data['precipitable_water'] = ds_clean
        self.raw_data.precipitable_water.attrs = dict(long_name = "Precipitable water",
                                                    units     = "cm",
                                                    comment   = "Equivalent liquid water depth; not CF-standard, convertible to 0.1 kg m-2 per mm.")
        return
        

    @property
    def od_co2_ch4_h2o(self):
        if isinstance(self._od_co2ch4h2o, type(None)):
            # get the 1625 filter info based on the MFRSR instrument serial no
            if 1625 in self.raw_data.channel:
                #### TODO: this file should be stored somewhere more meaning full
                if not isinstance(self.mfrsr_history, type(None)):
                    fn = '/home/grad/htelg/projects/AOD_redesign/MFRSR_History.xlsx'
                    mfrsr_info = pd.read_excel(self.mfrsr_history, sheet_name='Overview')
                    inst_info = mfrsr_info[mfrsr_info.Instrument == self.raw_data.serial_no]
                    assert(len(inst_info) == 1), f'Could not find unique MFRSR instrument info for serial no {self.raw_data.serial_no} in MFRSR history file {self.mfrsr_history}. Make sure the correct serial number is stored in raw_data.attrs["serial_no"]'
                    self.tp_inst_info = inst_info.copy()
                    fab = inst_info.Filter_1625nm.iloc[0]
                    self.tp_fab = fab
                    filter_no = int(''.join([i for i in fab if i.isnumeric()]))
                    filter_batch = ''.join([i for i in fab if not i.isnumeric()]).lower()
                    self.tp_filter_no = filter_no
                    self.tp_filter_batch = filter_batch
                    # open the lookup dabel for the Optical depth correction
                    # correction_info = xr.open_dataset(self.path2absorption_correction_ceoff_1625)
                    with xr.open_dataset(self.path2absorption_correction_ceoff_1625) as correction_info:
                        correction_info.load()
                
                    ds = xr.Dataset()
                    params_dict = {}
                    for molecule in ['co2', 'ch4', 'h2o_5cm']:
                        params = correction_info.sel(filter_no = filter_no, batch = filter_batch, molecule = molecule)
                        params_dict[molecule] = params
                        # apply the airmass dependence
                        da = params.OD.interp({'airmass': self.sun_position.airmass})
                        da = da.assign_coords(airmass = self.sun_position.index.values)
                        da = da.rename({'airmass': 'datetime'})
                        da = da.drop(['filter_no', 'molecule', 'batch'])
                        ds[molecule] = da
                    
                    # self.tp_params = params_dict
                    # self.tp_params_interp = ds.copy()
                    
                    ds = ds.expand_dims({'channel': [1625,]})
                    # normalize to the ambiant pressure ... less air -> less absorption, 
                    # only for ch4 and c02, water is scaled by the precipitable water
                    # TODO, this can probably done better? 
                    ds[['ch4', 'co2']] = ds[['ch4', 'co2']] * (self.met_data.pressure/1013.25)
                
                    # self.tp_ds = ds.copy()
                    # self.tp_tpw = self.tpw.copy()
                    tpw = self.precipitable_water
                    ds['h2o_5cm'] = ds.h2o_5cm / 5 * tpw
                
                    #### add 0 for all other channels
                    dstlist = [ds]
                    for cha in self.raw_data.channel:
                        if int(cha) == 1625:
                            continue    
                        
                        dst = copy.deepcopy(ds)
                        dst = dst.assign_coords({'channel': [cha]})
                    
                        for var in dst:
                            dst[var][:] = 0
                    
                        dstlist.append(dst)
                    ds = xr.concat(dstlist, 'channel')
                else:
                    assert(False), 'Set the mfrsr_history attribute to the path of the MFRSR history file'
                
            else:
                # generate a ds with zeros
                ds =  xr.Dataset()
                ds['co2'] = xr.DataArray(np.zeros(tuple(self.raw_data.dims[d] for d in ['datetime', 'channel'])), 
                                         coords={'datetime':self.raw_data.datetime, 'channel':self.raw_data.channel})
                ds['ch4'] = ds.co2.copy()
                ds['h2o_5cm'] = ds.co2.copy()
                
            self._od_co2ch4h2o = ds
        return self._od_co2ch4h2o

    @property
    def od_total(self):
        od = - np.log(self.transmission)/self.sun_position.airmass.to_xarray()
        return od
    
    @property
    def aod(self):
        if 'aod' not in self.raw_data:
            odt = self.od_total
            odr = self.od_rayleigh
            self.raw_data['aod'] = odt - odr #this generates a temporary aod variable, which is needed in case the pwv retrieva:l uses the aod at 940nm channel.
            odch4 = self.od_co2_ch4_h2o.ch4
            odco2 = self.od_co2_ch4_h2o.co2
            odh2o = self.od_co2_ch4_h2o.h2o_5cm
            odozone = self.od_ozone
            aod = odt - odr - odch4 - odco2 - odh2o - odozone
            # self._aod = aod
            self.raw_data['aod'] = aod
            
        return self.raw_data.aod
        

    # This is now split up in the SolarIrradiatin instance and here somewhere
    def _transmission_from_toa_radiation(self):
        """
        This assumes that the dataset already contains calibrated spectral radiation data
        in the direct_normal_irradiance variable. It will then get transmission using the toa radiation

        Returns
        -------
        None.

        """
        #### get solar spectrum at top of atmosphere
        if 'toa_spectral_irradiance' in self.raw_data:
            trans = self.raw_data.direct_normal / self.raw_data.toa_spectral_irradiance
            self.raw_data['transmission'] = trans
        else:
            ds_solar = self.toa_solarspectrum_1AU# xr.open_dataset(self.path2solar_spectrum)
            rad_top = ds_solar.interp(wavelength=self.dataset.channel_wavelength, method='linear')
            direct_es_corr = self.dataset.direct_normal * self.sun_position.sun_earth_distance.to_xarray()**2 #correct for earth sun distance, this requires mulitplication not division!!!!
            self.tp_direct_es_corr = direct_es_corr
            trans = direct_es_corr/rad_top #transmission
            trans = trans.drop_vars('wavelength')
            trans = trans.rename_vars({'irradiance':'transmission'})
            self.raw_data['transmission'] = trans.transmission
        return

    @property
    def transmission(self):
        # if isinstance(self._transmission, type(None)):
        if 'transmission' not in self.raw_data:
            if self.raw_data.attrs.get('calibrated_langley', False) and self.raw_data.global_horizontal.attrs.get('unit', 'bla') == 'W/m^2/nm':
                self._transmission_from_toa_radiation()
            elif self.settings_calibration == 'johns':
                self._apply_calibration_johns()
            elif self.settings_calibration == 'atm_gam':
                self._apply_calibration_atm_gam()
            elif self.settings_calibration == 'sp02':
                self._apply_calibration_sp02()
            elif self.settings_calibration == 'toa_radiation':
                self._transmission_from_toa_radiation()
            else:
                assert(False), f'Unknown calibration strategy - {self.settings_calibration}'
            if 0:
                #### Deprecated!!! below is the old sp02 retrieval, remove when sp02 retrieval is adapted
                #### load calibrations
                #### TODO: remove once SP02 is working again
                calibrations = atmlangcalib.load_calibration_history()
                cal = calibrations[int(self.raw_data.serial_no.values)]
                # use the mean and only the actual channels, other channels are artefacts
                cal = cal.results['mean'].loc[:,self.raw_data.channle_wavelengths.values].sort_index()
                
                #### interpolate and resample calibration (V0s)
                dt = self.raw_data.datetime.to_pandas()
                calib_interp = pd.concat([cal,dt]).drop([0], axis = 1).sort_index().interpolate().reindex(dt.index)
                
                #### correct VOs for earth sun distance see functions above
                calib_interp_secorr = calib_interp.divide(self.sun_position.sun_earth_distance**2, axis = 0)
                
                #### match channels for operation
                channels = self.raw_data.channle_wavelengths.to_pandas()
                raw_data = self.raw_data.raw_data.to_pandas().rename(columns = channels)
                raw_data.columns.name = 'wl'
                
                #### get transmission
                self._transmission = raw_data/calib_interp_secorr
        return self.raw_data.transmission
    
    
    # @property
    # def sun_position(self):
    #     if isinstance(self._sun_position, type(None)):
    #         self._sun_position = self.site.get_sun_position(self.raw_data.datetime)
    #     return self._sun_position
    
    @property
    def langley_am(self):
        if isinstance(self._am, type(None)):
            self._get_langley_from_raw() 
        return self._am
    
    @property
    def langley_pm(self):
        if isinstance(self._pm, type(None)):
            self._get_langley_from_raw() 
        return self._pm
    
    # def tp_get_rdl(self):
    #     raw_df = self.raw_data.raw_data.to_pandas()
        
    #     # changing to local time
    #     raw_df_loc = raw_df.copy()
    #     index_local = raw_df.index + pd.to_timedelta(self.site.time_zone[1], 'h')
    #     raw_df_loc.index = index_local
    #     self.raw_df_loc = raw_df_loc
        
    
    def _get_langley_from_raw(self, verbose = False):
        """Converts the current direct normal data into Langleys. Including:
        - correction for Sun-Earth distance
        - separaton into AM and PM
        Make sure the entire AM and PM are in the dataset, merge multiple files if necessary

        Returns
        -------
        None.

        Yields
        ------
        data of the langley_am and langley_pm properties, both of which are Langley instances.

        """

        raw_df = self.raw_data[self.variable_name_direct_normal].to_pandas()
        
        if verbose:
            print('change to local time')
        #### changing to local time
        raw_df_loc = raw_df.copy()
        index_local = raw_df.index + pd.to_timedelta(self.site.time_zone['diff2UTC_of_standard_time'], 'h')
        raw_df_loc.index = index_local
        # self.tp_rdl = raw_df_loc.copy()
        
        if verbose:
            print('get the one day')
        ##### getting the one day
        sunpos = self.sun_position.copy()
        start = raw_df_loc.index[0]
        if sunpos.iloc[0].airmass > 0:
            start = pd.to_datetime(f'{start.year}{start.month:02d}{start.day:02d}') + pd.to_timedelta(1,'d')
        end = start + pd.to_timedelta(1, 'd')
        raw_df_loc = raw_df_loc.truncate(start, end)

        if verbose:
            print('localize and cut day ')
        #### localize and cut day for sunposition
        sunpos.index = index_local
        sunpos = sunpos.truncate(start, end)
        if verbose:
            print('remove night')
        #### remove the night
        # return sunpos
        # sunpos[sunpos.airmass < 0] = np.nan
        sunpos = sunpos[sunpos.airmass > 0]
        if verbose:
            print('get min airmass')
        #### get the minimum airmass befor I start cutting it out
        noon = sunpos.airmass.idxmin()
        
        if verbose:
            print('normalize to sun earth dist.')
        #### normalize to the sun_earth_distance
        raw_df_loc = raw_df_loc.multiply(sunpos.sun_earth_distance**2, axis=0)
        
        if verbose:
            print('remove negatives')
        # langleys are the natural logarith of the voltage over the AMF ... -> log
        # to avoid warnings and strange values do some cleaning before log
        raw_df_loc[raw_df_loc <= 0] = np.nan
#         self.tp_raw_df = raw_df.copy()
        raw_df_loc = np.log(raw_df_loc)    
    
        if verbose:
            print('keep whats in right airmasses')
        # keep only what is considered relevant airmasses
        amf_min = 2.2 
        amf_max = 4.7
        # sunpos[sunpos.airmass < amf_min] = np.nan
        # sunpos[sunpos.airmass > amf_max] = np.nan
        sunpos = sunpos[sunpos.airmass > amf_min]
        sunpos = sunpos[sunpos.airmass < amf_max]
        if verbose:
            print('split into am/pm')
        sunpos_am = sunpos.copy()
        sunpos_pm = sunpos.copy()

        # sunpos_am[sunpos.index > noon] = np.nan
        # sunpos_pm[sunpos.index < noon] = np.nan
        sunpos_am = sunpos_am[sunpos.index < noon]
        sunpos_pm = sunpos_pm[sunpos.index > noon]

        # langley_am = raw_df_loc.copy()
        # langley_pm = raw_df_loc.copy()
        langley_am = raw_df_loc.loc[sunpos_am.index].copy()
        langley_pm = raw_df_loc.loc[sunpos_pm.index].copy()

        # self.tp_sp_am = sunpos_am
        # self.tp_sp_pm = sunpos_pm
        # self.tp_df_am = langley_am[~sunpos_am.airmass.isna()].copy()
        # self.tp_df_pm = langley_am[~sunpos_pm.airmass.isna()].copy()
        if verbose:
            print('set the values')
        # return langley_am, sunpos_am
        langley_am.index = sunpos_am.airmass
        langley_am = langley_am[~langley_am.index.isna()]
        langley_am.sort_index(ascending=False, inplace=True)
        
        langley_pm.index = sunpos_pm.airmass
        langley_pm = langley_pm[~langley_pm.index.isna()]
        langley_pm.sort_index(ascending=False, inplace=True)

        self._am = atmlangcalib.Langley(self,langley_am, langley_fit_settings = self.langley_fit_settings, airmass_limits = self.settings_langley_airmass_limits, when = 'am')
        self._pm = atmlangcalib.Langley(self,langley_pm, langley_fit_settings = self.langley_fit_settings, airmass_limits = self.settings_langley_airmass_limits, when = 'pm')
        if verbose:
            print('done')
        return True

class CombinedGlobalDiffuseDirect(SolarIrradiation):
    def __init__(self, dataset):
        """
        A class that combines the three irradiation components into one dataset and provides some specific functions
        Parameters
        ----------
        dataset : xr.Dataset
            A dataset that contains the three components as variables named 'global_horizontal', 'diffuse_horizontal', 'direct_normal'
            and the coordinates 'datetime' and 'channel'
        -----------
        2024-06-10, HTelg   
        """

        super().__init__(dataset)
        self._global_horizontal_irradiation = None
        self._diffuse_horizontal_irradiation = None
        self._direct_normal_irradiation = None
        # self.global_horizontal_irradiation = GlobalHorizontalIrradiation(dataset)
        # self.diffuse_horizontal_irradiation = DiffuseHorizontalIrradiation(dataset)
        # self.direct_normal_irradiation = DirectNormalIrradiation(dataset)
        
    # def apply_calibration(self, 
    #                       path2cal_spectral,
    #                       path2cal_logger = None,
    #                       path2cal_head = None,
    #                       pat2cal_cos = None,
    #                       ):
    #     #### spectral calibration
    #     cal = atmcal.read_mfrsr_cal(path2cal_spectral)
    #     si = self.apply_calibration_spectral(cal)
        
    #     #### logger amplifyer calibration
    #     cal = atmcal.read_factory_cal(path2cal_logger)
    #     si = self.apply_calibration_responsivity(cal)
        
    #     #### head calibration
        
    #     # remaining calibrations for lamp calibration
    #     mfrsr_fact = mfrsr.apply_calibration_responsivity(mfrsr_cal_fac, 
    #                                                            varname_responsivity_spectral = 'logger_gain',
    #                                                            varname_dark_signal_spectral = 'logger_darksignal',
    #                                                            # ignore_has_been_applied_error = True,
    #                                                           ) #logger calibration from factory cal ==> counts to voltage
    #     mfrsr_fact = mfrsr_fact.apply_calibration_responsivity(mfrsr_cal_fac, 
    #                                                            varname_responsivity_spectral = 'head_gain',
    #                                                            varname_dark_signal_spectral = 'head_darksignal',
    #                                                            ignore_has_been_applied_error = True,
    #                                                             ) #head calibration with lamp
    
    @property
    def global_horizontal_irradiation(self):
        if isinstance(self._global_horizontal_irradiation, type(None)):
            self._global_horizontal_irradiation = GlobalHorizontalIrradiation(self.dataset)
        return self._global_horizontal_irradiation
    
    @property
    def diffuse_horizontal_irradiation(self):
        if isinstance(self._diffuse_horizontal_irradiation, type(None)):
            self._diffuse_horizontal_irradiation = DiffuseHorizontalIrradiation(self.dataset)
        return self._diffuse_horizontal_irradiation

    @property
    def direct_normal_irradiation(self):
        if isinstance(self._direct_normal_irradiation, type(None)):
            self._direct_normal_irradiation = DirectNormalIrradiation(self.dataset)
        return self._direct_normal_irradiation

    def plot_overview(self, channel = 500, ax = None, 
                      show_alltime = False,
                      show_sunelevation = False):
        
        if isinstance(ax, type(None)):
            f, a= _plt.subplots()    
        else:
            a = ax
            f = a.get_figure()
            
        dssel = self.dataset.sel(channel = channel)
        
        if show_alltime:
            dssel.alltime.plot(ax = a, label = 'alltime')
        
        dssel.global_horizontal.plot(ax = a, label = 'global_horizontal')
        dssel.diffuse_horizontal.plot(ax = a, label = 'diffuse_horizontal')
        dssel.direct_normal.plot(ax = a, label = 'direct')
        
        if show_sunelevation:
            at = a.twinx()
            np.rad2deg(self.sun_position.elevation).plot(ax = at, color = 'black', ls = '--')
            # at.set_ylim(top = 0.9, bottom = 0)
        
        # a.set_xlim(left = pd.to_datetime('20220103 14:00:00'))
        a.grid()
        a.legend()
        return f,a
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
