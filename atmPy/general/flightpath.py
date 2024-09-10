from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import warnings
# try:
#     from mpl_toolkits.basemap import Basemap
# except KeyError:
#     warnings.warn("An error accured while trying to import mpl_toolkits.basemap.Basemap. Plotting of maps will not work!")
# try:
#     # from geopy.distance import vincenty
#     import geopy
# except ModuleNotFoundError:
#     warnings.warn('geopy not installed. You might encounter some functionality limitations.')
# import matplotlib.pylab as plt
# import pandas as pd
# import atmPy.general.timeseries as timeseries
# from functools import wraps as _wraps
from atmPy.opt_imports import mpl_toolkits_basemap as Basemap
from atmPy.opt_imports import geopy







# def get_path(hdf, datetime):
#     lon = hdf.select('Longitude')[:][:, 0]
#     lat = hdf.select('Latitude')[:][:, 0]
#     loc_ts = Path(pd.DataFrame({'Lon': lon, 'Lat': lat}, index=datetime))
#     return loc_ts
#
# def get_datetime(hdf):
#     time = hdf.select('Profile_UTC_Time')
#     time.dimensions()
#     td = time[:].transpose()[0]
#     day, date = np.modf(td)
#     datetime = pd.to_datetime(date.astype(int), format = '%y%m%d') + pd.to_timedelta(day, unit='d')
#     return datetime

# @_wraps(timeseries.TimeSeries())
# class Path(timeseries.TimeSeries):
#     @_wraps(timeseries.TimeSeries.__init__)
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
class FlightPath(object):
    def __init__(self, parent_ts, column_lat ='Lat', column_lon ='Lon', column_altitude ='Altitude'):
        self._parent_ts = parent_ts
        self._column_lat = column_lat
        self._column_lon = column_lon
        self._column_alt = column_altitude

        self.data = parent_ts.data.loc[:, [column_altitude, column_lat, column_lon]]

    def get_distance2location(self, location):
        """Ads a collumn to path.data.
        Args
        ----
        location: tuple
            (lat,lon)"""

        def get_dist(row, location):
            """
            Args
            ----
            location: tuple
                (lat,lon)
                """
            # dist = vincenty(location, (row['Lat'], row['Lon']))
            dist = geopy.distance.geodesic(location, (row['Lat'], row['Lon'])) # not sure if the syntax is right here 
            #     print(row)
            return dist.km

        dist = self.data.apply(lambda x: get_dist(x, location), axis=1)
        self.data['distance'] = dist
        return


    def get_closest2location(self, location):
        """
        Args
        ----
        location: tuple
            (lat,lon)"""

        def get_dist(row, location):
            """
            Args
            ----
            location: tuple
                (lat,lon)
                """
            dist = geopy.distance.geodesic(location, (row['Lat'], row['Lon']))
            #     print(row)
            return dist.km

        self.get_distance2location(self, location)
        data = self.data.copy()
        closest = data.sort_values('distance').iloc[[0]]
        return closest

    def plot_on_map(self, bmap_kwargs=None, estimate_extend=True, costlines=True, continents=True, background=None,
                    points_of_interest=None, three_d=False, verbose=False):
        """

        Args:
            bmap_kwargs:
            estimate_extend:
            costlines:
            background: str [None]
                options are 'bluemarble', 'shadedrelief', 'etopo'
            points_of_interest:
            three_d:
            verbose:

        Returns:

        """

        # data = self.data.copy()
        # data = data.loc[:, ['Lon', 'Lat']]
        data = self.data.dropna()

        if not bmap_kwargs:
            bmap_kwargs = {'projection': 'aeqd',
                           'resolution': 'c'}

        if estimate_extend:
            if bmap_kwargs['projection'] == 'aeqd':
                lon_center = (data[self._column_lat].values.max() + data[self._column_lon].values.min()) / 2.
                lat_center = (data.Lat.values.max() + data.Lat.values.min()) / 2.

                points = np.array([data.Lat.values, data.Lon.values]).transpose()
                distances_from_center_lat = np.zeros(points.shape[0])
                distances_from_center_lon = np.zeros(points.shape[0])
                for e, p in enumerate(points):
                    distances_from_center_lat[e] = geopy.distance.geodesic(p, (lat_center, p[1])).m
                    distances_from_center_lon[e] = geopy.distance.geodesic(p, (p[0], lon_center)).m

                lat_radius = distances_from_center_lat.max()
                lon_radius = distances_from_center_lon.max()
                scale = 1
                border = scale * 2 * np.array([lat_radius, lon_radius]).max()

                height = border + lat_radius
                width = border + lon_radius
                bmap_kwargs['lat_0'] = lat_center
                bmap_kwargs['lon_0'] = lon_center
                bmap_kwargs['width'] = width
                bmap_kwargs['height'] = height

            elif bmap_kwargs['projection'] == 'mill':
                bmap_kwargs['llcrnrlat'] = -90
                bmap_kwargs['urcrnrlat'] = 90
                bmap_kwargs['llcrnrlon'] = -180
                bmap_kwargs['urcrnrlon'] = 180

            if verbose:
                print('bmap_kwargs: {}'.format(bmap_kwargs))

        bmap = Basemap(**bmap_kwargs)
        if not three_d:
            wcal = np.array([161., 190., 255.]) / 255.

            if background == 'bluemarbel':
                bmap.bluemarble()
            elif background == 'etopo':
                bmap.etopo()
            elif background == 'shadedrelief':
                bmap.shadedrelief()
            elif not background:
                # Fill the globe with a blue color
                bmap.drawmapboundary(fill_color=wcal)
            else:
                raise ValueError('"{}" is not a valid background!'.format(background))

            if continents:
                grau = 0.9
                bmap.fillcontinents(color=[grau, grau, grau], lake_color=wcal)

            if costlines:
                bmap.drawcoastlines()

            x, y = bmap(data.Lon.values, data.Lat.values)
            path = bmap.plot(x, y,
                             color='m')

        # return bmap

        else:
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.add_collection3d(bmap.drawcoastlines())
            x, y = bmap(self.data.Lon.values, self.data.Lat.values)
            # ax.plot(x, y,self.data.Altitude.values,
            #           color='m')
            N = len(x)
            for i in range(N - 1):
                color = plt.cm.jet(i / N)
                ax.plot(x[i:i + 2], y[i:i + 2], self.data.Altitude.values[i:i + 2],
                        color=color)

        if points_of_interest:
            for point in points_of_interest:

                lat, lon = point.pop('loc_lat_lon')
                x, y = bmap(lon, lat)
                ax = plt.gca()

                try:
                    annotation = point.pop('annotation')
                    annotation_kwargs = point.pop('annotations_kwargs')
                except KeyError:
                    annotation = None

                g, = bmap.plot(x, y)
                g.set(**point)

                if annotation:
                    anno = ax.annotate(annotation, (x, y), **annotation_kwargs)
        return bmap


    # def __init__(self, parent_ts):)
    #     self._closest_location = None
    #     self._closest = None
    #     self._radius = None
    #     self._in_radius = None
    #
    # plot_on_map = plot_on_map
    #
    # def get_closest2location(self, location):
    #     if location == self._closest_location:
    #         return self._closest
    #     else:
    #         self._closest_location = location
    #         self._closest = get_closest2location(self, location)
    #         return self._closest
    #
    # def get_all_in_radius(self, location, radius):
    #     self.get_closest2location(location)
    #     if (location != self._closest_location) or (radius != self._radius):
    #         self._closest_location = location
    #         self._radius = radius
    #         self._in_radius = self.data.distance <= radius
    #     return self._in_radius
