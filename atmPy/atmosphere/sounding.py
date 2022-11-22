from atmPy.general import timeseries as _timeseries
from atmPy.general import flightpath as _flightpath
import xarray as xr
import pandas as pd
import metpy
import metpy.calc

class BalloonSounding(object):
    def __init__(self, data, column_lat='Lat', column_lon='Lon', column_altitude='Altitude'):
        if isinstance(data, xr.core.dataarray.Dataset):
            self.data = data
        else:
            # some old code ...not sure how valuable
            self.timeseries = _timeseries.TimeSeries(data)
            self.vertical_profile = self.timeseries.convert2verticalprofile()
            self.flight_path = _flightpath.FlightPath(self.timeseries, column_lat=column_lat, column_lon=column_lon, column_altitude=column_altitude)
            
        self._tpw = None
            
            
    @property
    def precipitable_water(self):
        if isinstance(self._tpw, type(None)):
            
            dser = pd.Series(index = self.data.site, dtype=float)
    
            for site in self.data.site:
                ds_sel = self.data.sel(site = site)
                # ds_sel = ds_sel.dropna('index') # just to double check, this made no difference
                pressure = ds_sel.pressure * metpy.units.units.hPa
                dewpoint = ds_sel.dewpoint * metpy.units.units.degC
                tpw = metpy.calc.precipitable_water(pressure, dewpoint)
                tpw = tpw.to('cm')
                dser.loc[str(site.values)] = tpw.magnitude
            
            datpw = dser.to_xarray()
            datpw = datpw.rename({'index': 'site'})
            datpw.attrs['unit'] = 'cm'
            self._tpw = datpw
        return self._tpw