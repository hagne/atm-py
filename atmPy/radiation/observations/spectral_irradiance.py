#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 16:30:35 2022

@author: hagen
"""
import numpy as np
import pandas as pd
import atmPy.general.measurement_site as atmms
import atmPy.radiation.observations.langley_calibration as atmlangcalib

class GlobalHorizontalIrradiation(object):
    def __init__(self, dataset):
        self.dataset = dataset

class DiffuseHorizontalIrradiation(object):
    def __init__(self, dataset):
        self.dataset = dataset

class DirectNormalIrradiation(object):
    def __init__(self, dataset, site = None, langley_fit_settings = None):
        self.raw_data = dataset #this is not exactly raw, it is cosine corrected voltag readings, thus, un-calibrated irradiances
        if isinstance(site, type(None)):
            assert('site' in dataset.attrs.keys()), 'If site is None, then the dataset has have lat,lon,site, site_name, attributes'
            self.site = atmms.Station(lat= dataset.attrs['site_latitude'], 
                                      lon = dataset.attrs['site_longitude'], 
                                      alt = dataset.attrs['site_elevation'], 
                                      name = dataset.attrs['site_name'], 
                                      abbreviation = dataset.attrs['site'],)
        else:
            self.site = site
        self.langley_fit_settings = langley_fit_settings
        self._sun_position = None
        self._am = None
        self._pm = None
        self._transmission = None
        # self._langleys_am = None
        # self._langleys_pm = None
        # self._langley_fitres_am = None
        # self._langley_fitres_pm = None
    

    #### TODO: New/changed, make it work!
    # this is more like apply calibration -> Make this a function that sets the self._transmission property
    def appply_calibration(calib):
        assert(False), 'work to be done here'
    
    @property
    def transmission(self):
        if isinstance(self._transmission, type(None)):
            assert(False), 'apply a langley calibration first'
            #### load calibrations
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
        return self._transmission
    
    @property
    def sun_position(self):
        if isinstance(self._sun_position, type(None)):
            self._sun_position = self.site.get_sun_position(self.raw_data.datetime)
        return self._sun_position
    
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
        
    
    def _get_langley_from_raw(self):
        raw_df = self.raw_data.direct_normal_irradiation.to_pandas()
        
        #### changing to local time
        raw_df_loc = raw_df.copy()
        index_local = raw_df.index + pd.to_timedelta(self.site.time_zone['diff2UTC_of_standard_time'], 'h')
        raw_df_loc.index = index_local
        # self.tp_rdl = raw_df_loc.copy()
        
        ##### getting the one day
        sunpos = self.sun_position.copy()
        start = raw_df_loc.index[0]
        if sunpos.iloc[0].airmass > 0:
            start = pd.to_datetime(f'{start.year}{start.month:02d}{start.day:02d}') + pd.to_timedelta(1,'d')
        end = start + pd.to_timedelta(1, 'd')
        raw_df_loc = raw_df_loc.truncate(start, end)

        #### localize and cut day for sunposition
        sunpos.index = index_local
        sunpos = sunpos.truncate(start, end)

        #### remove the night
        sunpos[sunpos.airmass < 0] = np.nan

        #### get the minimum airmass befor I start cutting it out
        noon = sunpos.airmass.idxmin()

        #### normalize to the sun_earth_distance
        raw_df_loc = raw_df_loc.multiply(sunpos.sun_earth_distance**2, axis=0)
    
        # langleys are the natural logarith of the voltage over the AMF ... -> log
        # to avoid warnings and strange values do some cleaning before log
        raw_df_loc[raw_df_loc <= 0] = np.nan
#         self.tp_raw_df = raw_df.copy()
        raw_df_loc = np.log(raw_df_loc)    
    
        # keep only what is considered relevant airmasses
        amf_min = 2.2 
        amf_max = 4.7
        sunpos[sunpos.airmass < amf_min] = np.nan
        sunpos[sunpos.airmass > amf_max] = np.nan

        sunpos_am = sunpos.copy()
        sunpos_pm = sunpos.copy()

        sunpos_am[sunpos.index > noon] = np.nan
        sunpos_pm[sunpos.index < noon] = np.nan
        

        langley_am = raw_df_loc.copy()
        langley_pm = raw_df_loc.copy()

        self.tp_sp_am = sunpos_am
        self.tp_sp_pm = sunpos_pm
        self.tp_df_am = langley_am[~sunpos_am.airmass.isna()].copy()
        self.tp_df_pm = langley_am[~sunpos_pm.airmass.isna()].copy()

        langley_am.index = sunpos_am.airmass
        langley_pm.index = sunpos_pm.airmass

        self._am = atmlangcalib.Langley(self,langley_am[~langley_am.index.isna()], langley_fit_settings = self.langley_fit_settings)
        self._pm = atmlangcalib.Langley(self,langley_pm[~langley_pm.index.isna()], langley_fit_settings = self.langley_fit_settings)
        return True

class CombinedGlobalDiffuseDirect(object):
    def __init__(self, dataset):
        self.dataset = dataset
        self.global_horizontal_irradiation = GlobalHorizontalIrradiation(dataset)
        self.diffuse_horizontal_irradiation = DiffuseHorizontalIrradiation(dataset)
        self.direct_normal_irradiation = DirectNormalIrradiation(dataset)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
