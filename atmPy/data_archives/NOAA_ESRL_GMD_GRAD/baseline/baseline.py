# -*- coding: utf-8 -*-
import pandas as pd
import xarray as xr
import atmPy.general.measurement_site as _measurement_site
import pathlib as pl

_locations = [{'name': 'Barrow',
              'state' :'AK',
              'abbreviation': ['BRW',],
              'lon': -156.6114,
              'lat': 71.3230,
              'alt' : 11,
              'timezone': 9},
              {'name': 'Mauna Loa',
              'state' :'HI',
              'abbreviation': ['MLO',],
              'lon': -155.5763,
              'lat': 19.5362,
              'alt' : 3397.,
              'timezone': 10},
              {'name': 'American Samoa',
              'state' :'Samoa',
              'abbreviation': ['SMO',],
              'lon': -170.5644,
              'lat': -14.2474,
              'alt' : 42.,
              'timezone': 11},
              {'name': 'South Pole',
              'state' :'',
              'abbreviation': ['SPO', ],
              'lon': 59,
              'lat': -90.00,
              'alt' : 2840,
              'timezone': -12}]

network = _measurement_site.Network(_locations)
network.name = 'GML Observatory Operations (OBOP)'

    

def read(path2file):
    #### get the units... taken from readme in  /nfs/iftp/aftp/g-rad/baseline/brw
    text = """station_name    character       station name, e. g., Goodwin Creek
latitude                real    latitude in decimal degrees (e. g., 40.80)
longitude               real    longitude in decimal degrees (e. g., 105.12)
elevation               integer elevation above sea level in meters
year                    integer year, i.e., 1995
jday                    integer Julian day (1 through 365 [or 366])
month                   integer number of the month (1-12)
day                     integer day of the month(1-31)
hour                    integer hour of the day (0-23)
min                     integer minute of the hour (0-59)
dt                      real    decimal time (hour.decimalminutes, e.g., 23.5 = 2330)
zen                     real    solar zenith angle (degrees)
dw_solar                real    downwelling global solar (Watts m^-2)
uw_solar                real    upwelling global solar (Watts m^-2)
direct_n                real    direct-normal solar (Watts m^-2)
diffuse                 real    downwelling diffuse solar (Watts m^-2)
dw_ir                   real    downwelling thermal infrared (Watts m^-2)
dw_casetemp             real    downwelling IR case temp. (K)
dw_dometemp             real    downwelling IR dome temp. (K)
uw_ir                   real    upwelling thermal infrared (Watts m^-2)
uw_casetemp             real    upwelling IR case temp. (K)
uw_dometemp             real    upwelling IR dome temp. (K)
uvb                     real    global UVB (milliWatts m^-2)
par                     real    photosynthetically active radiation (Watts m^-2)
netsolar                real    net solar (dw_solar - uw_solar) (Watts m^-2)
netir                   real    net infrared (dw_ir - uw_ir) (Watts m^-2)
totalnet                real    net radiation (netsolar+netir) (Watts m^-2)
temp                    real    10-meter air temperature (?C)
rh                      real    relative humidity (%)
windspd                 real    wind speed (ms^-1)
winddir                 real    wind direction (degrees, clockwise from north)
pressure                real    station pressure (mb)"""

    lines = [l.split() for l in text.split('\n')]
    lines = [[l[0],' '.join(l[2:])] for l in lines]
    units = pd.DataFrame(lines[4:], columns = ['col_name', 'unit'])
    units.index = units.col_name
    units = units.drop('col_name', axis = 1)
    # units = units.iloc[7:]
    
    
    #### get column labels, from /nfs/iftp/aftp/g-rad/baseline/brw
    text = 'year,jday,month,day,hour(i),min(i),dt(i),zen(i),dw_solar(i),qc_dwsolar(i),uw_solar(i),qc_uwsolar(i),direct_n(i),qc_direct_n(i),diffuse(i),qc_diffuse(i),dw_ir(i),qc_dwir(i),dw_casetemp(i),qc_dwcasetemp(i),dw_dometemp(i),qc_dwdometemp(i),uw_ir(i),qc_uwir(i),uw_casetemp(i),qc_uwcasetemp(i),uw_dometemp(i),qc_uwdometemp(i),uvb(i),qc_uvb(i),par(i),qc_par(i),netsolar(i),qc_netsolar(i),netir(i),qc_netir(i),totalnet(i),qc_totalnet(i),temp(i),qc_temp(i),rh(i),qc_rh(i),windspd(i),qc_windspd(i),winddir(i),qc_winddir(i),pressure(i),qc_pressure(i)'
    column_labels = [col.strip().replace('(i)', '') for col in text.split(',')]
    column_labels = ['minute' if col == 'min' else col for col in column_labels]
    
    #### open file
    p2f = pl.Path(path2file)
    site = p2f.name[:3]
    data = pd.read_csv(p2f, delim_whitespace=True, skiprows=2,
                names = column_labels,
               )
    
    data.index = data.apply(lambda row: pd.to_datetime(f'{row.year:0.0f}-{int(row.month):02d}-{int(row.day):02d} {int(row.hour):02d}:{int(row.minute):02d}:00'), axis = 1)
    data = data.drop(['year', 'jday', 'month', 'day', 'hour', 'minute','dt'], axis = 1)
    data.index.name = 'datetime'
    # return data
    ds = xr.Dataset(data)
    
    for var in ds.variables:
        if var == 'datetime':
            continue
        # break
        elif 'qc' in var:
            continue
    
        ut = units.loc[var,'unit']
        ds[var].attrs['unit'] = ut
    
    site_inst = network.stations.find_site(site)
    ds.attrs['site_latitude'] = site_inst.lat
    ds.attrs['site_longitude'] = site_inst.lon
    ds.attrs['site_elevation'] =  site_inst.alt
    ds.attrs['site_name'] = site_inst.name
    ds.attrs['site'] =site
    ds.attrs['instrument_type'] = 'sp02'    
    return ds
    
    