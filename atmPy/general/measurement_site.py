import warnings
try:
    from mpl_toolkits.basemap import Basemap as _Basemap
except:
    warnings.warn('There seams to be an issue with importing mpl_toolkits.basemap. Make sure it is installed and working (try: "from mpl_toolkits.basemap import Basemap as _Basemap"). For now plotting on map will not be possible')
import matplotlib.pylab as _plt
import os as _os
import numpy as _np
import pandas as _pd
from atmPy.radiation import solar as _solar
# import datetime as _datetime
import timezonefinder as _tzf
import pytz as _pytz
# The following is to ensure that one can use as large of an image as one desires
try:
    from PIL import Image as _image
    _image.MAX_IMAGE_PIXELS = None
except ModuleNotFoundError:
    warnings.warn('PIL not installed. This is needed to plot images of arbitrary resolution, e.g. when plotting satellite images.')
from matplotlib import path as _path
try:
    import geopy as _geopy
    # import geopy.distance as _gd #otherwise distance will not be available
except ModuleNotFoundError:
    warnings.warn('geopy not installed. You might encounter some functionality limitations.')
import plt_tools as _plt_tools

default_colors = _plt.rcParams['axes.prop_cycle'].by_key()['color']

class NetworkStations(object):
    def __init__(self):
        self._stations_list = []
        
    def __iter__(self):
        for st in self._stations_list:
            yield st
            
    def find_site(self, site):
        """
        Searches for "site" in site abbriviations (programming requrired if seach in other keys is desired) and returns the site instance

        Parameters
        ----------
        site : TYPE
            DESCRIPTION.

        Returns
        -------
        res : TYPE
            DESCRIPTION.

        """
        res = [si for si in self._stations_list if site.upper() in si.abb]
        if len(res) == 0:
            res = [si for si in self._stations_list if site.upper() in si.abb]
            
        assert(len(res) != 0), 'not found'
        return res[0]
    
    @property
    def list(self):
        return self._stations_list

class SubNetworks(object):
    def __init__(self):
        self._network_list = []

    def plot(self, colors = None, **kwargs):
        if isinstance(colors, type(None)):
            site_label_colors = default_colors
        else:
            site_label_colors = colors
        # if isinstance(kwargs['station_symbol_kwargs'], type(None)):
        if ('station_symbol_kwargs' not in kwargs) or (isinstance(kwargs['station_symbol_kwargs'], type(None))):
            kwargs['station_symbol_kwargs'] = {}
        # else:
        #     site_label_colors = colors
        # if 'site_label_marker_size' not in kwargs:
        #     kwargs['site_label_marker_size'] = 8

        # labels = []
        for e,network in enumerate(self._network_list):
            kwargs['station_symbol_kwargs']['color'] = site_label_colors[e]
            a, bmap = network.plot(**kwargs)
            kwargs['bmap'] = bmap
            g = a.get_lines()[-1]
            g.set_label(network.name)
            # label = _plt.Line2D([0], [0], linestyle='', marker='o', label=network.name, markerfacecolor=site_label_colors[e], markersize= )
            # labels.append(label)

        # if 'site_label_font_size' in kwargs:
        #     fontsize = kwargs['site_label_font_size']
        # else:
        #     fontsize = None
        # a.legend(handles=labels, fontsize = fontsize)
        return a, bmap

class Network(object):
    def __init__(self, network_stations, generate_subnetworks = True, name = None, abbr = None, by_keys = ['name']):
        """Generate a network of stations
        Parameters
        ----------
        network_stations: list of dicts
            make sure keys of dicts follow the rules outlined in the Station class which is called for
            each station in network_stations
        by_key: str
            stations within the network are givin by name by the self.stations attribute. In addition stations can be
            given by a given key which must be a key in network_stations"""
        self.name = name
        self.abbr = abbr
        self.stations = NetworkStations()
        self._networkstation_by_key_list = [{'networkstation_inst': self.stations, 'key': 'name'}]
        if len(by_keys) > 1:
            for key in by_keys:
                nwsinst_name = 'station_by_{}'.format(key)
                setattr(self, nwsinst_name, NetworkStations())
                nwsinst = getattr(self, nwsinst_name)
                self._networkstation_by_key_list.append({'networkstation_inst': nwsinst, 'key':key, 'name': nwsinst_name})

        if isinstance(network_stations, list):
            if isinstance(network_stations[0], dict):
                network_stations = network_stations.copy()
                for station in network_stations:
                    station = station.copy()
                    if isinstance(station['abbreviation'], list):
                        station['abbreviation'] = station['abbreviation'][0]
                    # else:
                    #     abb = station['abbreviation']
                    station['name'] = station['name']#.replace(' ', '_').replace('.', '')
                    # site = Station(lat=station['lat'],
                    #                lon=station['lon'],
                    #                alt=station['alt'],
                    #                name=name,
                    #                abbreviation=abb,
                    #                info=None)
                    site = Station(**station)
                    self.add_station(site)
                    # if len(by_keys) > 1:
                    #     for key in by_keys:
                    #         self.add_station(site, key = key)

            else:
                for station in network_stations:
                    self.add_station(station)

        if generate_subnetworks:
            self._operation_period2sub_network_active()

    @property
    def extend(self):
        extend = {}
        st_lat_max = max(self.stations._stations_list, key=lambda x: x.lat)
        st_lat_min = min(self.stations._stations_list, key=lambda x: x.lat)
        st_lon_max = max(self.stations._stations_list, key=lambda x: x.lon)
        st_lon_min = min(self.stations._stations_list, key=lambda x: x.lon)
        center = (st_lat_max.lat + st_lat_min.lat) / 2, (st_lon_max.lon + st_lon_min.lon) / 2
        extend['center'] = center

        height = _geopy.distance.geodesic((st_lat_min.lat, center[1]), (st_lat_max.lat, center[1]))
        width = _geopy.distance.geodesic((center[0], st_lon_min.lon), (center[0], st_lon_max.lon))

        extend['width_m'] = width.m
        extend['height_m'] = height.m
        return extend

    def add_station(self, site_instance, key = 'name'):
        site_instance.parent_network = self
        attrname = getattr(site_instance, key)
        attrname = attrname.replace(' ', '_').replace('.', '')
        setattr(self.stations, attrname, site_instance)
        self.stations._stations_list.append(site_instance)
        # not sure what the folling is to achieve
        # for network in self._networkstation_by_key_list:
        #     inst = network['networkstation_inst']
        #     key = network['key']
        #     # name = network['name']
        #     setattr(inst, getattr(site_instance,key), site_instance)

    def add_subnetwork(self, network_instance):
        if not hasattr(self, 'subnetworks'):
            self.subnetworks = SubNetworks()
        setattr(self.subnetworks, network_instance.name, network_instance)
        self.subnetworks._network_list.append(network_instance)

    def plot(self, zoom = 1.1, network_label_format = '{abbr}', network_label_kwargs = False,  **kwargs):
        """Plots all station in network on map

        Parameters
        ----------
        network_label_kwargs: dict or False
            font size, bbox, postion ... of the label. If False then no label is printed
            defaults = dict(xytext=(0, 0),
                            size=18,
                            ha="center",
                            va='center',
                            textcoords='offset points',
                            bbox=dict(boxstyle="round", fc=[1, 1, 1, 0.5], ec='black'),
                            )
        kwargs: dictionary with arguments that are passed to station.plot

            """
        stl = self.stations._stations_list
        # a, bmap = stl[0].plot(plot_only_if_on_map = True, **kwargs)
        # kwargs['bmap'] = bmap

        if 'center' not in kwargs.keys():
            kwargs['center'] = 'auto'
        if kwargs['center'] == 'auto':
            # st = self.stations._stations_list[0]
            # lat_max, lat_min = max([st.lat for st in self.stations._stations_list]), min(
            #     [st.lat for st in self.stations._stations_list])
            # lon_max, lon_min = max([st.lon for st in self.stations._stations_list]), min(
            #     [st.lon for st in self.stations._stations_list])
            # kwargs['center'] = ((lat_max+lat_min)/2, (lon_max+lon_min)/2)
            kwargs['center'] = self.extend['center']

        if 'width' not in kwargs.keys():
            kwargs['width'] = 'auto'
        if kwargs['width'] == 'auto':
            kwargs['width'] = self.extend['width_m'] * zoom

        if 'height' not in kwargs.keys():
            kwargs['height'] = 'auto'
        if kwargs['height'] == 'auto':
            kwargs['height'] = self.extend['height_m'] * zoom

        if 'station_symbol_kwargs' not in kwargs.keys():
            kwargs['station_symbol_kwargs'] = {}

        for e, station in enumerate(stl):
            if 'color' not in kwargs['station_symbol_kwargs']:
                kwargs['station_symbol_kwargs']['color'] = default_colors[1]
            a, bmap = station.plot(plot_only_if_on_map = True, **kwargs)
            kwargs['bmap'] = bmap

        if network_label_kwargs != False:
            annodefaults = dict(xytext=(0, 0),
                                size=12,
                                ha="center",
                                va='center',
                                textcoords='offset points',
                                bbox=dict(boxstyle="round", fc=[1, 1, 1, 0.5], ec='black'),
                                )
            if isinstance(network_label_kwargs, type(None)):
                network_label_kwargs = {}

            for ak in annodefaults:
                if ak not in network_label_kwargs:
                    network_label_kwargs[ak] = annodefaults[ak]

            label = network_label_format.format(abbr=self.abbr, name=self.name)
            center_lat = sum([station.lat for station in stl]) / len(stl)
            center_lon = sum([station.lon for station in stl]) / len(stl)
            xpt, ypt = bmap(center_lon, center_lat)
            a.annotate(label, xy=(xpt, ypt), **network_label_kwargs
                       #                 xycoords='data',
                       # xytext=(10 ,-10),
                       # size = station_annotation_kwargs['size'],
                       # ha="left",
                       # va = 'top',
                       # textcoords='offset points',
                       # bbox=dict(boxstyle="round", fc=[1 ,1 ,1 ,0.5], ec='black'),
                       )



        return a,bmap

    def _operation_period2sub_network_active(self):
        has_operation_period = _np.array([hasattr(sta, 'operation_period') for sta in self.stations._stations_list])
        # test if any/all have operation_period
        if not _np.any(has_operation_period):
            return False
        assert (_np.all(has_operation_period))

        # skip if all or none have 'present' in period
        has_present = _np.array(['present' in sta.operation_period for sta in self.stations._stations_list])
        if ((_np.all(has_present)) | (has_present.sum() == 0)):
            return False

        lt = [sta for sta in self.stations._stations_list if 'present' in sta.operation_period]
        active = Network(lt, generate_subnetworks=False)
        active.name = 'active'
        self.add_subnetwork(active)

        lt = [sta for sta in self.stations._stations_list if 'present' not in sta.operation_period]
        inactive = Network(lt, generate_subnetworks=False)
        inactive.name = 'inactive'
        self.add_subnetwork(inactive)


class Station(object):
    def __init__(self, lat = None, lon = None, alt = None, name = None, abbreviation = None, active = None,
                 operation_period = None, info = None, state = '', country = '', parent_network = None, **kwargs):
        """
        Generates a Station instance
        Parameters
        ----------
        lat
        lon
        alt
        name
        abbreviation
        active
        operation_period
        info
        kwargs
        """
        self.lat = lat
        self.lon = lon
        self.alt = alt
        self.name = name
        self.abb = abbreviation
        self.state = state
        self.country = country

        self.parent_network = parent_network

        if info:
            self.info = info
        if active:
            self.active = active
        if operation_period:
            self.operation_period = operation_period
        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.name = self.name.strip().replace(',','_').replace('(','_').replace(')','_').strip('_')
        
        self._time_zone = None
        
    def __repr__(self):
        txt = (f'name: {self.name} ({self.abb})\n'
               f'lat/lon/alt: {self.lat}/{self.lon}/{self.alt}')
        # print(txt)
        return txt
    
    @property
    def time_zone(self, date ='now'):
        # get timezone
        tz = _tzf.TimezoneFinder()
        tz_str = tz.timezone_at(lng = self.lon, lat = self.lat)
        
        # now = _datetime.datetime.now(_pytz.timezone(tz_str))
        # tz_hr = now.utcoffset() / _datetime.timedelta(hours = 1)
        if date == 'now':
            now = _pd.Timestamp.now(_pytz.timezone(tz_str))
        else:
            now = _pd.to_datetime(date).tz_localize(_pytz.timezone(tz_str))
        out = {}
        out['tz_str'] = tz_str
        out['tz_name'] = now.tzname()
        out['diff2UTC'] = now.utcoffset() / _pd.to_timedelta(1, unit = 'h')
        out['diff2UTC_of_standard_time'] = (now.utcoffset() / _pd.to_timedelta(1,'h')) - (now.dst()/_pd.to_timedelta(1,'h'))
        return out
        
    
    def get_sun_position(self, datetime):
        out = _solar.get_sun_position(self.lat, self.lon, datetime, elevation = self.alt)
        return out
        

    def plot(self,
             projection='lcc',
             center = 'auto',
             width=400000 * 7,
             height=500000 * 3,
             station_label_format ='{abbr}',
             station_label_kwargs = None,
             resolution='c',
             background='blue_marble',
             station_symbol_kwargs = None,
             # site_label_marker_size = 8,
             # site_label_font_size = 18,
             # site_label_color='auto',
             bmap = None,
             plot_only_if_on_map = False,
             ax = None,
             verbose = False):
        """

        Parameters
        ----------
        projection
        center: 'auto' or (lat, lon)
        width
        height
        station_label_format: format str ('abbr', 'name', 'state')
            This takes a fromat string with the given optional arguments. E.g. '{name}, {state}'.
        station_label_kwargs: dict or bool
            This defines the position, size ... of the label. If False no label will
            be shown. See doc of plt.annotate() for details.
            defaults = dict(xytext = (10, -10),
                            size = 18,
                            ha = "left",
                            va = 'top',
                            textcoords = 'offset points',
                            bbox = dict(boxstyle="round", fc=[1, 1, 1, 0.5], ec='black'),
                             )
        resolution: str ('c','i','h'....
        background: str
            blue_marble: use the blue_marble provided by basemap
            "path_to_file_name": This will use the warpimage function to use the image in the filename ... works with the blue marble stuff (https://visibleearth.nasa.gov/view_cat.php?categoryID=1484)
        plot_only_if_on_map: bool
            as the name says
        ax

        Returns
        -------

        """
        if isinstance(station_symbol_kwargs, type(None)):
            station_symbol_kwargs = {}
        if 'marker' not in station_symbol_kwargs:
            station_symbol_kwargs['marker'] = 'o'
        if 'markersize' not in station_symbol_kwargs:
            station_symbol_kwargs['markersize'] = 4
        if 'color' not in station_symbol_kwargs:
            station_symbol_kwargs['color'] = default_colors[1]
        if 'zorder' not in station_symbol_kwargs:
            station_symbol_kwargs['zorder'] = 100
            

        if bmap:
            # a = bmap.ax
            a = _plt.gca()
        else:
            if ax:
                a = ax
            else:
                f ,a = _plt.subplots()

            #         self.lat = 36.605700
            #         self.lon = -97.487846
            # self.lat = 78.5
            # self.lon = 18
            # width = 34000
            # height = 22000

            #         width = 400000 * 7
            #         height = 500000 * 3

            if center == 'auto':
                lat = self.lat
                lon = self.lon
            else:
                lat, lon = center

            bmap = _Basemap  (  # projection='laea',
                projection=projection,
                lat_0=lat,
                lon_0=lon,
                width=width,
                height=height,
                resolution=resolution,
                ax = a
            )
            # bmap.drawcoastlines()
            if background == 'blue_marble':
                bmap.bluemarble()
            elif isinstance(background, type(None)):
                pass
            else:
                assert(_os.path.isfile(background))
                bmap.warpimage(image=background)

            bmap.drawcountries(linewidth=2)
            bmap.drawstates()

        # out = bmap.fillcontinents()
        # wcal = np.array([161., 190., 255.]) / 255.
        # boundary = bmap.drawmapboundary(fill_color=wcal)

        lon, lat = self.lon, self.lat # Location Ny-Alesund 78°13′N 15°33′E
        # convert to map projection coords.
        # Note that lon,lat can be scalars, lists or numpy arrays.

        if plot_only_if_on_map:
            map_bound_path = _path.Path(_np.array([bmap.boundarylons, bmap.boundarylats]).transpose())
            if not map_bound_path.contains_point((lon, lat)):
                if verbose:
                    txt = 'Station {} was not plot, since it is not on the map.'.format(self.name)
                    print(txt)
                return a,bmap

        xpt ,ypt = bmap(lon ,lat)
        p, = bmap.plot(xpt ,ypt ,linestyle = '',**station_symbol_kwargs)
        # if site_label_color == 'auto':
        #     col = colors[1]
        # else:
        #     col = site_label_color

        # p.set_color(station_symbol_kwargs['color'])
        # p.set_markersize()

        # if station_label == 'abbr':
        #     label = self.abb
        # elif station_label == 'name':
        #     label = self.name
        # elif station_label == 'label':
        #     label = self.label
        # else:
        #     raise ValueError('{} is not an option for station_label'.format(station_label))


    # Station Label
        if station_label_kwargs != False:
            annodefaults = dict(xytext = (10, -10),
                                size = 10,
                                ha = "left",
                                va = 'top',
                                textcoords = 'offset points',
                                bbox = dict(boxstyle="round", fc=[1, 1, 1, 0.5], ec='black'),
                                zorder = 100
                                 )
            if isinstance(station_label_kwargs, type(None)):
                station_label_kwargs = {}

            for ak in annodefaults:
                if ak not in station_label_kwargs:
                    station_label_kwargs[ak] = annodefaults[ak]

            label = station_label_format.format(abbr=self.abb, name=self.name, state=self.state, country=self.country)


            a.annotate(label, xy=(xpt, ypt), **station_label_kwargs
                           #                 xycoords='data',
                           # xytext=(10 ,-10),
                           # size = station_annotation_kwargs['size'],
                           # ha="left",
                           # va = 'top',
                           # textcoords='offset points',
                           # bbox=dict(boxstyle="round", fc=[1 ,1 ,1 ,0.5], ec='black'),
                           )

    # return
        self._bmap = bmap
        self._xpt, self._ypt = xpt, ypt
        return a,bmap

    def plot_square_around_site(self, edgelength=10, print_distance = True, color = 'white',
                                print_label = None):
        """
        Parameters
        ----------
        edgelength: in km
        print_label: list
            prints the distance to the edge of the box.
            default: ['left', 'bottom']
        """
        col1 = color
        if isinstance(print_label, type(None)):
            print_label = ['left', 'bottom']
        lat = self.lat
        lon = self.lon
        #     bmap = ...
        point_c = (lat, lon)

        dist = _geopy.distance.distance()  # thats just to get an instance of distance, can change the actual distance later in the destination function

        # calculate distance from center to corner of square
        #     el = 10 # edge length of square
        dist2corner = _np.sqrt(2) * edgelength / 2

        # compass bearings to the corners
        bearings = [(i * 90) - 45 for i in (1, 2, 3, 4)]

        # calculate the coordinates of the corners
        corners = []
        for bear in bearings:
            point = dist.destination(point_c, bear, distance=dist2corner)
            corners.append((point.longitude, point.latitude))

        # convert coordinates to x, y
        corners_pt = [self._bmap(*cco) for cco in corners]

        #######
        # repeat for distance text
        # compass bearings to the corners
        label_bearings = [(i * 90)  for i in (0, 1, 2, 3)]

        # calculate the coordinates of the corners
        label_pos = []
        for bear in label_bearings:
            point = dist.destination(point_c, bear, distance=edgelength/2)
            label_pos.append((point.longitude, point.latitude))

        # convert coordinates to x, y
        label_pos_pt = [self._bmap(*cco) for cco in label_pos]

        # plot it
        square = _plt.Polygon(corners_pt)
        square.set_linestyle('--')
        square.set_edgecolor(col1)
        square.set_facecolor([1, 1, 1, 0])
        a = _plt.gca()
        a.add_artist(square)

        #######
        ## plot the distanced on the edges of the square
        labels = []
        label_kwargs = dict(color = col1,
                            fontsize = 'medium',
                            textcoords = 'offset points',
                            )
        text_dist = 5
        txt = f'{edgelength} km'
        if 'top' in print_label:
            lab = label_pos_pt[0]
            anno = a.annotate(txt, (lab[0], lab[1]),
                        xytext=(0, text_dist),va = 'bottom', ha = 'center', **label_kwargs)
            labels.append(anno)

        if 'bottom' in print_label:
            lab = label_pos_pt[2]
            anno = a.annotate(txt, (lab[0], lab[1]),
                              xytext=(0, - text_dist), va='top', ha='center', **label_kwargs)
            labels.append(anno)

        if 'left' in print_label:
            lab = label_pos_pt[3]
            anno = a.annotate(txt, (lab[0], lab[1]),
                              xytext=(- text_dist, 0), va='center', ha='right',
                              rotation=90,
                              **label_kwargs)
            labels.append(anno)

        if 'right' in print_label:
            lab = label_pos_pt[1]
            anno = a.annotate(txt, (lab[0], lab[1]),
                              xytext=(text_dist, 0), va='center', ha='left',
                              rotation=90,
                              **label_kwargs)
            labels.append(anno)
            # a.text(lab[0], lab[1],
            #        # 'bla',
            #        f'{edgelength} km',
            #        color = col1,
            #        va = 'bottom', ha = 'center',
            #        xytext=(0, 35))

        # for lab in label_pos_pt:
        #     a.text(lab[0], lab[1], 'bla')
        # a.plot(x,y, zorder = 100, marker = 'o')
        return square

    def plot_distance2coordinate(self, point):
        """ Plot a line between the site and an arbitrary point

        Parameters
        ----------
        point: (lat, lon) or Station instance

        Returns
        -------
        graph, text, bbox
        """
        if isinstance(point, Station):
            pass
        elif len(point) == 2:
            point = Station(lat = point[0], lon=[1])

        g, = self._bmap.plot([self._xpt, point._xpt], [self._ypt, point._ypt], color='white', zorder=1,
                             linestyle='--',
                             linewidth=2)
        col = g.get_color()
        dist = round(_geopy.distance.geodesic((self.lat, self.lon), (point.lat, point.lon)).km)
        txt = f'{dist} km'
        txt, fbp = _plt_tools.text.add_text_along_graph(g, txt, (self._xpt + point._xpt) / 2,
                                                       txtkwargs=dict(va='bottom',
                                                                      ha='center',
                                                                      color=col,
                                                                      fontsize='x-large'),
                                                       bbox=False
                                                       )
        return g, txt, fbp