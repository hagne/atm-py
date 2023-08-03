#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 10:45:55 2023

@author: htelg
"""
import pathlib as pl
import xarray as xr
import pandas as pd
import requests
import configparser
import pyorbital
import pyorbital.orbital
import numpy as np

satellite_info_list = [dict(name = 'NOAA 20', id = 43013, sensors = ['VIIRS']),]
sensor_info_list = [dict(name = 'VIIRS', viewing_angle = 56)]

def get_satellite_info(satellite):
    satt = [sat for sat in satellite_info_list if sat['name'] == satellite]
    assert(len(satt) == 1), f'found {len(satt)} satellites with that name! Choose from following list or update the satellite_info list.\n{[s["name"] for s in satellite_info_list]}'
    return satt[0]

def get_sensor_info(sensor, satellite = None):
    sens = [s for s in sensor_info_list if s['name'] == sensor]
    assert(len(sens) == 1), f'found {len(sens)} satellites with that name! Choose from following list or update the sensor_info list.\n{[s["name"] for s in sensor_info_list]}'
    if not isinstance(satellite, type(None)):
        sat = get_satellite_info(satellite)
        assert(sensor in sat['sensors']), f'The {sensor} sensor is not on the {satellite} satellite. Correct or update stored info.'
    return sens[0]

class Orbit(object):
    def __init__(self, satellite = 'NOAA 20', year = 2019, download_when_missing = False, verbose = True):
        """
        https://www.space-track.org/documentation#howto-api_python

        Parameters
        ----------
        satellite : TYPE, optional
            DESCRIPTION. The default is 'NOAA_20'.
        year: int
            the year to use, currently the orbit information is stored in yearly files
            FIXIT, the above will probably have changed
        download_when_missing : TYPE, optional
            DESCRIPTION. The default is False.
        verbose : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        None.

        """
        self.satellite = satellite
        self.verbose = verbose
        self.year = 2019
        self.start = f'{self.year}-01-01'
        self.end = f'{self.year + 1}-01-01'
        self.download_when_missing = download_when_missing
        self.fname = pl.Path(f"orbit_tle-data_{satellite.replace(' ', '-')}_{self.year}.nc")
        


        ### get satellite info
        # satt = [sat for sat in self._satellite_info_list if sat['name'] == satellite]
        # assert(len(satt) == 1), f'found {len(satt)} satellites with that name! Choose from following list or update the satellite_info list.\n{[s["name"] for s in self._satellite_info_list]}'
        self.satellite_info = get_satellite_info(satellite)
        
        self._tle_data = None
        

    @property
    def tle_data(self):
        if self.fname.is_file():
            self._tle_data = xr.open_dataarray(self.fname)
        else:
            if self.download_when_missing:
                if self.verbose:
                    print('downloading tle data')
                # self._tle_data = download_tle_data()
            else:
                assert(False), 'tle_data file not found and download_when_missing is set to False. Change the latter to download the tle data.'
        return self._tle_data
            
        
    def download_tle_data(self, save = True):

        
        ### some settings that should not change
        uriBase                = "https://www.space-track.org"
        requestLogin           = "/ajaxauth/login"
        tle_baseurl = "https://www.space-track.org/basicspacedata/query/class/gp_history/NORAD_CAT_ID/{satellite_id}/orderby/TLE_LINE1%20ASC/EPOCH/{start}--{end}/format/tle"
        
        

        # satellite_info
        

        
        ### format and check time range
        dtstart = pd.to_datetime(self.start)
        dtend = pd.to_datetime(self.end)
        ddtd = (dtend - dtstart)/pd.to_timedelta(1, 'd')
        assert(ddtd >= 1), f'differnt has to be at least 1 day, is: {ddtd:0.2f} days' 
        
        # satellite_id = 43013
        satellite_id = self.satellite_info['id']
        
        ### Use configparser package to pull in the ini file (pip install configparser)
        config = configparser.ConfigParser()
        config.read("./SLTrack.ini")
        configUsr = config.get("configuration","username")
        configPwd = config.get("configuration","password")
        # configOut = config.get("configuration","output")
        siteCred = {'identity': configUsr, 'password': configPwd}
        
        # format start and end time for url, just in case the times got more complicated
        dtstartstr = f'{dtstart.year:04d}-{dtstart.month:02d}-{dtstart.day:02d}'
        dtendstr = f'{dtend.year:04d}-{dtend.month:02d}-{dtend.day:02d}'
        
        ### get the tle data
        with requests.Session() as session:
            #login
            resp = session.post(uriBase + requestLogin, data = siteCred)
            if resp.status_code != 200:
                assert(False) ,"POST fail on login"
            url_to_download = tle_baseurl.format(satellite_id = satellite_id, start = dtstartstr, end = dtendstr)    
            # url_to_download = "https://www.space-track.org/basicspacedata/query/class/gp_history/NORAD_CAT_ID/43013/orderby/TLE_LINE1%20ASC/EPOCH/2023-08-01--2023-08-02/format/tle"
            # url_to_download = "https://www.space-track.org/basicspacedata/query/class/gp_history/NORAD_CAT_ID/43013/orderby/TLE_LINE1%20ASC/EPOCH/2023-08-01--2023-08-02/format/tle"
            response = session.get(url_to_download)
        
        ### format the for later usage
        text = response.text.replace('\r\n', '\n') # window to linux
        text = text.strip() #prevents last line in df to be empty
        
        lines = text.split('\n')
        
        # lines = [l.strip() for l in lines]
        # lines = [l.replace('+', ' ') for l in lines]
        
        df = pd.DataFrame(lines, columns = ['tle'])
        
        df
        
        df['line'] = df.apply(lambda row: int(row.tle.split()[0]), axis = 1)
        def get_timestamp(row):
            if row.line == 1:
                dtt = row.tle.split()[3]
                timestamp = pd.to_datetime(dtt.split('.')[0], format = '%y%j') + pd.to_timedelta(float(f"0.{dtt.split('.')[1]}"), 'd')
            else:
                timestamp = pd.NaT
            return timestamp
        
        dts = df.apply(get_timestamp, axis = 1)
        df.index = dts.ffill()
        df.index.name = 'datetime'
        df.index = pd.MultiIndex.from_arrays([df.index, df.line])
        # df = df.drop_duplicates() # for some reasons there where duplicates
        df = df[~df.index.duplicated(keep='first')]
        tle_data = df.tle.to_xarray()
        if save:
            tle_data.to_netcdf(self.fname)
        return tle_data
    
    

class OverPasses(object):
    def __init__(self, satellite = 'NOAA 20', sensor = 'VIIRS', start = '20220802', end = None,
                 site = dict(lon = -105.2705, lat = 40.0150, alt = 1.5,), #boulder 
                 verbose = True):

        self.satellite = satellite        
        self.satellite_info = get_satellite_info(satellite)
        self.sensor = sensor
        self.sensor_info = get_sensor_info(sensor, satellite = satellite)
        self.start = pd.to_datetime(start)
        
        
        self.site = type('adhoc_site', (), site)
        
        if isinstance(end, type(None)):
            self.end = self.start + pd.to_timedelta(1, 'd')
        else:
            self.end = pd.to_datetime(end)
            
        assert(self.start.year == self.end.year), 'start and end in different years, currently not supported, fix it!!'
        
        self._orbit = None
        self._overpasses = None
            
    @property 
    def orbit(self):
        if isinstance(self._orbit, type(None)):
            self._orbit = Orbit(satellite = self.satellite)
        return self._orbit
    
    @property
    def overpasses(self):
        if isinstance(self._overpasses, type(None)):
            # out = {}
            tle_data = self.orbit.tle_data
            tstart = self.start
            tend = self.end
            
            ddt = pd.to_datetime(tle_data.datetime) - tstart
            ddt_closest = abs(ddt).min()/ pd.to_timedelta(1,'d')
            assert(ddt_closest <= 1), f'There is no tle information within 1 day to the date of interes. Closes tle file is {ddt_closest:0.2f} days'
            
            dt_idx = abs(ddt).argmin()
            
            tle_line1 = str(tle_data.isel(datetime = dt_idx).sel(line = 1).values)
            tle_line2 = str(tle_data.isel(datetime = dt_idx).sel(line = 2).values)
            
            sat = pyorbital.orbital.Orbital("NOAA-20", line1=tle_line1, line2=tle_line2)
            
            # sat.get_equatorial_crossing_time(tstart, tend, node='ascending', local_time=False, rtol=1e-09)
            
            utc_time = tstart
            length = (tend - tstart) / pd.to_timedelta(1, 'hour')
            length = int(np.ceil(length))
            
            overpasses = sat.get_next_passes(utc_time, length, self.site.lon, self.site.lat, self.site.alt, tol=0.001, horizon=0)
            #### select highest point (ignore rise and setting time) 
            overpasses = [op[2] for op in overpasses]
            overpasses = pd.DataFrame(overpasses, columns=['overpass_time_utc'])            
            overpasses['overpass_time_local'] = overpasses.apply(lambda row: sat.utc2local(row.overpass_time_utc), axis = 1)
            
            #### get observer viewing angles
            obs = np.array([sat.get_observer_look(utc_time, self.site.lon, self.site.lat, self.site.alt) for utc_time in overpasses.overpass_time_utc])
            overpasses = pd.concat([overpasses,pd.DataFrame(obs, columns=['obs_azimus_angle', 'obs_elevation_angle'])], axis=1)
            # out['obs'] = obs
            
            #### select overpasses if observer is in viewing angle. Note, very simple estimate, curvature of earth is ignored!!!
            if not isinstance(self.sensor, type(None)):
            #     overpasses = [[utc_time,(90 - obsel)[e]]   for e,utc_time in enumerate(overpasses) if (90 - obsel[e]) < self.sensor_info['viewing_angle']]
                overpasses['observed_by_sensor'] = overpasses.apply(lambda row: (90 - row.obs_elevation_angle) < self.sensor_info['viewing_angle'], axis = 1)

            # out['overpasses'] = overpasses
            # out['sat'] = sat
            
            overpasses.index.name = 'overpass_idx'
            overpasses = overpasses.to_xarray()
            self._overpasses = overpasses
        return self._overpasses 
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                