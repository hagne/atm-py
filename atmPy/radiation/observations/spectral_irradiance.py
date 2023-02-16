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
import atmPy.radiation.observations.langley_calibration as atmlangcalib
import atmPy.radiation.rayleigh.lab as atmraylab
import pathlib as pl
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.surfrad.surfrad as atmsrf
import copy

class GlobalHorizontalIrradiation(object):
    def __init__(self, dataset):
        self.dataset = dataset

class DiffuseHorizontalIrradiation(object):
    def __init__(self, dataset):
        self.dataset = dataset

class DirectNormalIrradiation(object):
    def __init__(self, dataset, 
                 site = None, 
                 langley_fit_settings = None,
                 calibration_strategy = 'johns',
                 metdata = 'surfrad'):
        
        self.raw_data = dataset #this is not exactly raw, it is cosine corrected voltag readings, thus, un-calibrated irradiances
        if isinstance(site, type(None)):
            assert('site' in dataset.attrs.keys()), 'If site is None, then the dataset has to have lat,lon,site, site_name, attributes'
            self.site = atmms.Station(lat= dataset.attrs['site_latitude'], 
                                      lon = dataset.attrs['site_longitude'], 
                                      alt = dataset.attrs['site_elevation'], 
                                      name = dataset.attrs['site_name'], 
                                      abbreviation = dataset.attrs['site'],)
        else:
            self.site = site
        self.langley_fit_settings = langley_fit_settings
        self.settings_calibration = calibration_strategy #
        self.settings_metdata = metdata
        self._sun_position = None
        self._am = None
        self._pm = None
        self._transmission = None
        # self._langleys_am = None
        # self._langleys_pm = None
        # self._langley_fitres_am = None
        # self._langley_fitres_pm = None
        self._od_rayleigh = None
        self._od_co2ch4h2o = None
        self._tpw = None
        self._aod = None
        self._metdata = None
        self._od_ozone = None
        self.path2absorption_correction_ceoff_1625 = '1625nm_absorption_correction_coefficience.nc'
    
    @property
    def met_data(self):
        if isinstance(self._metdata, type(None)):
            if self.settings_metdata == 'surfrad':
                self._metdata = self._get_surface_met_data()
            else:
                assert(False), 'moeeeep!'
        return self._metdata
    
    def _get_surface_met_data(self):
        """
        Loads the surfrad data (from my netcdf version), interpolates and adds 
        to the dataset. Location is currently hardcoded, this will likey cause problems sooner or later

        Returns
        -------
        None.

        """
        start_dt = pd.to_datetime(self.raw_data.direct_normal_irradiation.datetime.min().values)
        end_dt = pd.to_datetime(self.raw_data.direct_normal_irradiation.datetime.max().values)
        
        # open relevant files
        site = self.site.abb
        fns = []
        for dt in [start_dt, end_dt]:
            fns.append(f'/mnt/telg/data/grad/surfrad/radiation/{site}/srf_rad_full_{site}_{dt.year:04d}{dt.month:02d}{dt.day:02d}.nc')
        
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
        
        #### match channels for operation
        channels = self.raw_data.channel_center.to_pandas()
        raw_data = self.raw_data.direct_normal_irradiation.to_pandas().rename(columns = channels)
        raw_data.columns.name = 'wl'
        
        #### get transmission
        transmission = raw_data/calib_interp_secorr
        transmission = xr.DataArray(transmission)
        self.calibration = 'sp02'
        self._transmission =transmission
            
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
        
        # self.tp_dt = dt
        # self.tp_langley_params = langley_params
        
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
        
        raw_data = self.raw_data.direct_normal_irradiation.to_pandas()
        transmission = raw_data/calib_interp_secorr
        transmission = xr.DataArray(transmission)     
        self.calibration = 'johns'
        self._transmission =transmission

    @property
    def od_rayleigh(self):
        if isinstance(self._od_rayleigh, type(None)):
            odr = xr.concat([atmraylab.rayleigh_od_johnsmethod(self.met_data.pressure, chan) for chan in self.raw_data.channel_center], 'channel')
            self._od_rayleigh = odr
        return self._od_rayleigh

    @property
    def od_ozone(self):
        if isinstance(self._od_ozone, type(None)):
            #### read ozon concentration file for particular site and extract data for this day
            p2f = f'/home/grad/surfrad/aod/ozone/{self.site.abb}_ozone.dat'
            ozone = atmsrf.read_ozon(p2f)
            dt = pd.to_datetime(pd.to_datetime(self.raw_data.datetime.values[0]).date())
            total_ozone = float(ozone.interp(datetime = dt).ozone)
            
            #### read spectral function of ozon absorption coeff and interpolate to exact filter wavelength
            p2f = '/home/grad/surfrad/aod/ozone.coefs'
            ozone_abs_coef = pd.read_csv(p2f, index_col= 0, sep = ' ', names=['wavelength', 'coeff'])
            ozone_abs_coef = ozone_abs_coef.to_xarray()
            
            ozon2bychannel = ozone_abs_coef.interp(wavelength = self.raw_data.channel_center)
            ozon2bychannel = ozon2bychannel.drop('wavelength')
            
            #### scale abs. coef. to total ozone to get OD_ozone
            od_ozone = ozon2bychannel * (total_ozone/1000)
            od_ozone = od_ozone.fillna(0)
            self._od_ozone = od_ozone.coeff
        return self._od_ozone
    
    @property
    def precipitable_water(self):
        if isinstance(self._tpw, type(None)):
            sitedict = {'Bondville': 'BND', 
                        'Fort Peck': 'FPK',
                        'Goodwin Creek': 'GWN', 
                        'Table Mountain': 'TBL',
                        'Desert Rock': 'DRA',
                        'Penn State': 'PSU', 
                        'ARM SGP': 'sgp', 
                        'Sioux Falls': 'SXF',
                        'Canaan Valley': 'CVA',
                       }
            
            #### get path to relevant soundings files
            files = pd.DataFrame()
            datetime = pd.to_datetime(self.raw_data.datetime.values[0])
            for dt in [datetime, datetime + pd.to_timedelta(1, 'day')]:
                # pass
            
                p2fld = pl.Path(f'/nfs/grad/surfrad/sounding/{dt.year}/')
                searchstr = f'{dt.year}{dt.month:02d}{dt.day:02d}*.int'
                df = pd.DataFrame(p2fld.glob(searchstr), columns=['p2f'])
                files = pd.concat([files, df])
            
            files.sort_values('p2f', inplace = True)
            
            
            tpwts = []
            # tpsoundis = []
            for p2f in files.p2f:
                # pass
            
                soundi = atmsrf.read_sounding(p2f)
                # tpsoundis.append({"sounding" : soundi, "fn":p2f})
                tpwts.append(soundi.precipitable_water.expand_dims({'datetime': [soundi.data.attrs['datetime'],]}))
                
            # self.tp_soundi = tpsoundis
            tpw = xr.concat(tpwts, 'datetime')
            # self.tp_tpw_all_1 = tpw.copy()
            tpw = tpw.assign_coords(site = [sitedict[str(s.values)].lower() for s in tpw.site])
            # self.tp_tpw_all_2 = tpw.copy()
            tpw = tpw.sel(site = self.raw_data.site)
            tpw = tpw.drop(['site'])
            
            tpw = tpw.interp({'datetime': self.raw_data.datetime})
            self._tpw =tpw
        return self._tpw
        

    @property
    def od_co2_ch4_h2o(self):
        if isinstance(self._od_co2ch4h2o, type(None)):
            # get the 1625 filter info based on the MFRSR instrument serial no
            fn = '/mnt/telg/projects/AOD_redesign/MFRSR_History.xlsx' 
            mfrsr_info = pd.read_excel(fn, sheet_name='Overview')
            inst_info = mfrsr_info[mfrsr_info.Instrument == self.raw_data.serial_no]
            
            fab = inst_info.Filter_1625nm.iloc[0]
            filter_no = int(''.join([i for i in fab if i.isnumeric()]))
            filter_batch = ''.join([i for i in fab if not i.isnumeric()]).lower()
            
            # open the lookup dabel for the Optical depth correction
            correction_info = xr.open_dataset(self.path2absorption_correction_ceoff_1625)
            
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
            self._od_co2ch4h2o = ds
        return self._od_co2ch4h2o

    @property
    def od_total(self):
        od = - np.log(self.transmission)/self.sun_position.airmass.to_xarray()
        return od
    
    @property
    def aod(self):
        if isinstance(self._aod, type(None)):
            odt = self.od_total
            odr = self.od_rayleigh
            odch4 = self.od_co2_ch4_h2o.ch4
            odco2 = self.od_co2_ch4_h2o.co2
            odh2o = self.od_co2_ch4_h2o.h2o_5cm
            odozone = self.od_ozone
            aod = odt - odr - odch4 - odco2 - odh2o - odozone
            self._aod = aod
            
        return self._aod
        
    @property
    def transmission(self):
        if isinstance(self._transmission, type(None)):
            if self.settings_calibration == 'johns':
                self._apply_calibration_johns()
            elif self.settings_calibration == 'sp02':
                self._apply_calibration_sp02()
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

        # self.tp_sp_am = sunpos_am
        # self.tp_sp_pm = sunpos_pm
        # self.tp_df_am = langley_am[~sunpos_am.airmass.isna()].copy()
        # self.tp_df_pm = langley_am[~sunpos_pm.airmass.isna()].copy()

        langley_am.index = sunpos_am.airmass
        langley_am = langley_am[~langley_am.index.isna()]
        langley_am.sort_index(ascending=False, inplace=True)
        
        langley_pm.index = sunpos_pm.airmass
        langley_pm = langley_pm[~langley_pm.index.isna()]
        langley_pm.sort_index(ascending=False, inplace=True)

        self._am = atmlangcalib.Langley(self,langley_am, langley_fit_settings = self.langley_fit_settings)
        self._pm = atmlangcalib.Langley(self,langley_pm, langley_fit_settings = self.langley_fit_settings)
        return True

class CombinedGlobalDiffuseDirect(object):
    def __init__(self, dataset):
        self.dataset = dataset
        self.global_horizontal_irradiation = GlobalHorizontalIrradiation(dataset)
        self.diffuse_horizontal_irradiation = DiffuseHorizontalIrradiation(dataset)
        self.direct_normal_irradiation = DirectNormalIrradiation(dataset)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
