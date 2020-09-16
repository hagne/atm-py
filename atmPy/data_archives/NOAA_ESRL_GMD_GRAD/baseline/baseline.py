# -*- coding: utf-8 -*-
# import numpy as _np
# import pandas as _pd
# import os as _os
# import atmPy.general.timeseries as _timeseries
# import atmPy.aerosols.physics.column_optical_properties as _column_optical_properties
import atmPy.general.measurement_site as _measurement_site
# import pathlib
# import warnings as _warnings

# from atmPy.general import measurement_site as _measurement_site
# from atmPy.radiation import solar as _solar

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
              'lat': 14.2474,
              'alt' : 42.,
              'timezone': 11},
              {'name': '',
              'state' :'',
              'abbreviation': ['', ''],
              'lon': 59,
              'lat': 90.00,
              'alt' : 2840,
              'timezone': -12}]

network = _measurement_site.Network(_locations)
network.name = 'GML Observatory Operations (OBOP)'
