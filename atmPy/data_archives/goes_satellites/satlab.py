# -*- coding: utf-8 -*-
import xarray as xr
# import pathlib as pl
import numpy as np
# import cartopy.crs as ccrs
# import metpy 
# from scipy import interpolate
# from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
from pyproj import Proj
import matplotlib.pyplot as plt
import plt_tools
from matplotlib import colors


class GeosSatteliteProducts(object):
    def __init__(self,file):
        if type(file) == xr.core.dataset.Dataset:
            ds = file
        else:
            ds = xr.open_dataset(file)
        
        self.ds = ds
        
#         self._varname4test = 'CMI_C02'
        self._lonlat = None
        
    @property
    def lonlat(self):
        if isinstance(self._lonlat, type(None)):
            # Satellite height
            sat_h = self.ds['goes_imager_projection'].perspective_point_height

            # Satellite longitude
            sat_lon = self.ds['goes_imager_projection'].longitude_of_projection_origin

            # Satellite sweep
            sat_sweep = self.ds['goes_imager_projection'].sweep_angle_axis

            # The projection x and y coordinates equals the scanning angle (in radians) multiplied by the satellite height
            # See details here: https://proj4.org/operations/projections/geos.html?highlight=geostationary
            x = self.ds['x'][:] * sat_h
            y = self.ds['y'][:] * sat_h

            # Create a pyproj geostationary map object to be able to convert to what ever projecton is required
            p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)

            # Perform cartographic transformation. That is, convert image projection coordinates (x and y)
            # to latitude and longitude values.
            XX, YY = np.meshgrid(x, y)
            lons, lats = p(XX, YY, inverse=True)
            
            # Assign the pixels showing space as a single point in the Gulf of Alaska
#             where = np.isnan(self.ds[self._varname4test].values)
            where = np.isinf(lons)
            lats[where] = 57
            lons[where] = -152

            self._lonlat = (lons, lats) #dict(lons = lons, 
#                                     lats = lats)
        return self._lonlat
            # Assign the pixels showing space as a single point in the Gulf of Alaska
    #             where = np.isnan(channels_rgb['red'])
    #             lats[where] = 57
    #             lons[where] = -152
    
class OR_ABI_L2_AODC_M6(GeosSatteliteProducts):
    def plot(self, bmap = None, colorbar = True, max_aod = 1, norm = 'linear'):
        cb = None
        lons, lats = self.lonlat
        if isinstance(bmap, type(None)):
            bmap = Basemap(resolution='c', projection='aea', area_thresh=5000, 
                             width=5e6, height=3e6, 
    #                          lat_1=38.5, lat_2=38.5, 
                             lat_0=38.5, lon_0=-97.5)

            bmap.drawcoastlines(linewidth = 0.7)
            bmap.drawcountries(linewidth = 0.7)
            bmap.drawstates(linewidth = 0.3)
        pc = bmap.pcolormesh(lons, lats, self.ds.AOD.values, linewidth=0, latlon=True, zorder = 10)
        pc.set_clim(vmax = max_aod)
        pc.set_cmap(plt.cm.inferno)
        pc.set_alpha(0.3)
        
        if norm == 'log':    
            pc.set_norm(colors.LogNorm(0.01,1))
       
        if colorbar:
            # f = plt.gcf()
            # f.colorbar(pc)
            
            a= plt.gca()
            cb, a = plt_tools.colorbar.colorbar_axis_split_off(pc,a)
            cb.set_label('Aerosol optical depth')
            cb.set_alpha(1)
            cb.draw_all()
        return bmap, pc, cb

class OR_ABI_L2_MCMIPC(GeosSatteliteProducts):
    def __init__(self, *args):
        super().__init__(*args)
#         self._varname4test = 'CMI_C02'

    
    def plot_true_color(self, 
                        gamma = 1.8,#2.2, 
                        contrast = 130, #105
                        projection = None,
                        bmap = None
                       ):
        channels_rgb = dict(red = self.ds['CMI_C02'].data.copy(),
                            green = self.ds['CMI_C03'].data.copy(),
                            blue = self.ds['CMI_C01'].data.copy())
        
        channels_rgb['green_true'] = 0.45 * channels_rgb['red'] + 0.1 * channels_rgb['green'] + 0.45 * channels_rgb['blue']

        
        
        for chan in channels_rgb:
            col = channels_rgb[chan]
            # Apply range limits for each channel. RGB values must be between 0 and 1
            new_col = col / col[~np.isnan(col)].max()
            
            # apply gamma
            if not isinstance(gamma, type(None)):
                new_col = new_col**(1/gamma)
            
            # contrast
            #www.dfstudios.co.uk/articles/programming/image-programming-algorithms/image-processing-algorithms-part-5-contrast-adjustment/
            if not isinstance(contrast, type(None)):
                cfact = (259*(contrast + 255))/(255.*259-contrast)
                new_col = cfact*(new_col-.5)+.5
            
            channels_rgb[chan] = new_col
            
        rgb_image = np.dstack([channels_rgb['red'],
                             channels_rgb['green_true'],
                             channels_rgb['blue']])
        rgb_image = np.clip(rgb_image,0,1)
            
        a = plt.subplot()
        if isinstance(projection, type(None)) and isinstance(bmap, type(None)):
            
            a.imshow(rgb_image)
#             a.set_title('GOES-16 RGB True Color', fontweight='semibold', loc='left', fontsize=12);
#             a.set_title('%s' % scan_start.strftime('%d %B %Y %H:%M UTC '), loc='right');
            a.axis('off')
    
        else:          
            lons,lats = self.lonlat
            
            # Make a new map object Lambert Conformal projection
            if not isinstance(bmap,Basemap):
                bmap = Basemap(resolution='i', projection='aea', area_thresh=5000, 
                             width=3000*3000, height=2500*3000, 
    #                          lat_1=38.5, lat_2=38.5, 
                             lat_0=38.5, lon_0=-97.5)

                bmap.drawcoastlines()
                bmap.drawcountries()
                bmap.drawstates()

            # Create a color tuple for pcolormesh

            # Don't use the last column of the RGB array or else the image will be scrambled!
            # This is the strange nature of pcolormesh.
            rgb_image = rgb_image[:,:-1,:]

            # Flatten the array, becuase that's what pcolormesh wants.
            colortuple = rgb_image.reshape((rgb_image.shape[0] * rgb_image.shape[1]), 3)

            # Adding an alpha channel will plot faster, according to Stack Overflow. Not sure why.
            colortuple = np.insert(colortuple, 3, 1.0, axis=1)

            # We need an array the shape of the data, so use R. The color of each pixel will be set by color=colorTuple.
            pc = bmap.pcolormesh(lons, lats, channels_rgb['red'], color=colortuple, linewidth=0, latlon=True, zorder = 0)
            pc.set_array(None) # Without this line the RGB colorTuple is ignored and only R is plotted.

#             plt.title('GOES-16 True Color', loc='left', fontweight='semibold', fontsize=15)
#             plt.title('%s' % scan_start.strftime('%d %B %Y %H:%M UTC'), loc='right');
