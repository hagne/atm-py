from atmPy.general import measurement_site as _measurement_site




_stations = [{'name': 'Alert',
              'state': 'Canada',
              'lat': 82.45,
              'lon': -62.52,
              'abbreviation': 'ALT',
              'alt': 210.0,
              'operation_period': '2004-present'},
             {'name': 'Anmyeon-do',
              'state': 'S. Korea',
              'lat': 36.54,
              'lon': 126.33,
              'abbreviation': 'AMY',
              'alt': 45.0,
              'operation_period': '2009-2015,2017-present'},
             {'name': 'Boone',
              'state': 'N. Carolina',
              'lat': 36.2,
              'lon': -81.7,
              'abbreviation': 'APP',
              'alt': 1100.0,
              'operation_period': '2009-present'},
             {'name': 'El Arenocillo',
              'state': 'Spain',
              'lat': 37.1,
              'lon': -6.73,
              'abbreviation': 'ARN',
              'alt': 41.0,
              'operation_period': '2009-present'},
             {'name': 'BEO-Moussala',
              'state': 'Hungary',
              'lat': 42.18,
              'lon': 23.59,
              'abbreviation': 'BEO',
              'alt': 2925.0,
              'operation_period': '2010-present'},
             {'name': 'Bondville',
              'state': 'Illinois',
              'lat': 40.05,
              'lon': -88.37,
              'abbreviation': 'BND',
              'alt': 230.0,
              'operation_period': '1996-present'},
             {'name': 'Barrow',
              'state': 'Alaska',
              'lat': 71.32,
              'lon': -156.6,
              'abbreviation': 'BRW',
              'alt': 11.0,
              'operation_period': '1976-present'},
             # {'name': 'Cape Grim',
             #  'state': 'Australia',
             #  'lat': -40.68,
             #  'lon': 144.69,
             #  'abbreviation': 'CGO',
             #  'alt': 94.0,
             #  'operation_period': '2015-present'},
             {'name': 'Cape San Juan',
              'state': 'Puerto Rico',
              'lat': 18.48,
              'lon': -66.13,
              'abbreviation': 'CPR',
              'alt': 230.0,
              'operation_period': '2004-present'},
             {'name': 'Cape Point',
              'state': 'South Africa',
              'lat': -34.35,
              'lon': 18.49,
              'abbreviation': 'CPT',
              'alt': 230.0,
              'operation_period': '2004-2005'},
             {'name': 'Egbert',
              'state': 'Canada',
              'lat': 44.23,
              'lon': -79.783,
              'abbreviation': 'EGB',
              'alt': 253.0,
              'operation_period': '2009-present'},
             {'name': 'BERMS',
              'state': 'Canada',
              'lat': 54.35,
              'lon': -104.99,
              'abbreviation': 'ETL',
              'alt': 540.0,
              'operation_period': '2008-present'},
             {'name': 'AMF Black Forest',
              'state': 'Germany',
              'lat': 48.54,
              'lon': 8.4,
              'abbreviation': 'FKB',
              'alt': 511.0,
              'operation_period': '2007'},
             {'name': 'AMF Graciosa Isl.',
              'state': 'Portugal',
              'lat': 39.09,
              'lon': -28.03,
              'abbreviation': 'GRW',
              'alt': 15.24,
              'operation_period': '2009-2010'},
             {'name': 'Gosan',
              'state': 'S. Korea',
              'lat': 33.28,
              'lon': 126.17,
              'abbreviation': 'GSN',
              'alt': 72.0,
              'operation_period': '2009-present'},
             {'name': 'AMF Shouxian',
              'state': 'China',
              'lat': 32.56,
              'lon': 116.78,
              'abbreviation': 'HFE',
              'alt': 22.7,
              'operation_period': '2008'},
             {'name': "K'puszta",
              'state': 'Hungary',
              'lat': 46.96,
              'lon': 19.583,
              'abbreviation': 'KPS',
              'alt': 125.0,
              'operation_period': '1994-1996,2008-present'},
             {'name': 'Lulin',
              'state': 'Taiwan',
              'lat': 23.47,
              'lon': 120.873,
              'abbreviation': 'LLN',
              'alt': 2862.0,
              'operation_period': '2008-present'},
             {'name': 'AMF Manacapuro',
              'state': 'Brazil',
              'lat': -3.21,
              'lon': -60.0,
              'abbreviation': 'MAN',
              'alt': 50.0,
              'operation_period': '2014-2015'},
             {'name': 'Mount Bachelor',
              'state': 'Oregon',
              'lat': 43.979,
              'lon': -121.687,
              'abbreviation': 'MBO',
              'alt': 2743.0,
              'operation_period': '2014-2015, 2018-present'},
             {'name': 'Mauna Loa',
              'state': 'Hawaii',
              'lat': 19.54,
              'lon': -155.58,
              'abbreviation': 'MLO',
              'alt': 3397.0,
              'operation_period': '1974-present'},
             {'name': 'Montsec',
              'state': 'Spain',
              'lat': 42.05,
              'lon': 0.73,
              'abbreviation': 'MSA',
              'alt': 3397.0,
              'operation_period': '2000-2021'},
             {'name': 'Montseny',
              'state': 'Spain',
              'lat': 41.78,
              'lon': 2.36,
              'abbreviation': 'MSY',
              'alt': 3397.0,
              'operation_period': '2009-2021'},
             {'name': 'AMF Niamey',
              'state': 'Niger',
              'lat': 13.48,
              'lon': 2.18,
              'abbreviation': 'NIM',
              'alt': 205.0,
              'operation_period': '2006'},
             {'name': 'AMF Nainital',
              'state': 'India',
              'lat': 29.36,
              'lon': 79.46,
              'abbreviation': 'PGH',
              'alt': 1951.0,
              'operation_period': '2011-2012'},
             {'name': 'AMF Cape Cod',
              'state': 'Massachusetts',
              'lat': 42.07,
              'lon': -70.2,
              'abbreviation': 'PVC',
              'alt': 1.0,
              'operation_period': '2012-2013'},
             {'name': 'AMF Pt Reyes',
              'state': 'California',
              'lat': 38.09,
              'lon': -122.96,
              'abbreviation': 'PYE',
              'alt': 5.0,
              'operation_period': '2005'},
             {'name': 'Resolute',
              'state': 'Canada',
              'lat': 74.71,
              'lon': -94.98,
              'abbreviation': 'RSL',
              'alt': 64.0,
              'operation_period': '2013-2017'},
             {'name': 'Seoul',
              'state': 'South Korea',
              'lat': 37.57,
              'lon': 126.966,
              'abbreviation': 'SEL',
              'alt': 85.8,
              'operation_period': '2015-2017'},
             {'name': 'Southern Great Plains',
              'state': 'OK',
              'lat': 36.61,
              'lon': -97.49,
              'abbreviation': 'SGP',
              'alt': 315.0,
              'operation_period': '1997-2017'},
             {'name': 'Samoa',
              'state': 'American Samoa',
              'lat': -14.23,
              'lon': -170.56,
              'abbreviation': 'SMO',
              'alt': 77.0,
              'operation_period': '1977-2017'},
             {'name': 'Hyytiala',
              'state': 'Finland',
              'lat': 61.847,
              'lon': 24.295,
              'abbreviation': 'SMR',
              'alt': 181.0,
              'operation_period': '2018-present'},
             {'name': 'Sierra Nevada Station',
              'state': 'Spain',
              'lat': 37.096,
              'lon': -3.387,
              'abbreviation': 'SNS',
              'alt': 2501.0,
              'operation_period': '2016-present'},
             {'name': 'Storm Peak',
              'state': 'Colorado',
              'lat': 40.45,
              'lon': -106.73,
              'abbreviation': 'SPL',
              'alt': 3220.0,
              'operation_period': '2011-present'},
             {'name': 'South Pole',
              'state': 'Antarctica',
              'lat': -89.98,
              'lon': -24.8,
              'abbreviation': 'SPO',
              'alt': 2410.0,
              'operation_period': '1974-present'},
             {'name': 'Summit',
              'state': 'Greenland',
              'lat': 72.58,
              'lon': -38.48,
              'abbreviation': 'SUM',
              'alt': 32.38,
              'operation_period': '2011-present'},
             {'name': 'Trinidad Head',
              'state': 'California',
              'lat': 41.05,
              'lon': -124.15,
              'abbreviation': 'THD',
              'alt': 107.0,
              'operation_period': '2002-2017'},
             {'name': 'Tiksi',
              'state': 'Russia',
              'lat': 71.6,
              'lon': 128.9,
              'abbreviation': 'TIK',
              'alt': 7.0,
              'operation_period': '2009'},
             {'name': 'Whistler',
              'state': 'Canada',
              'lat': 50.01,
              'lon': -122.95,
              'abbreviation': 'WHI',
              'alt': 2182.0,
              'operation_period': '2008-present'},
             {'name': 'Mt Waliguan',
              'state': 'China',
              'lat': 36.28,
              'lon': 100.9,
              'abbreviation': 'WLG',
              'alt': 3810.0,
              'operation_period': '2005-present'},
             {'name': 'Sable Island',
              'state': 'Nova Scotia',
              'lat': 43.93,
              'lon': -60.01,
              'abbreviation': 'WSA',
              'alt': 5.0,
              'operation_period': '1992-2000'},
             {'name': 'Ny-Alesund',
              'state': 'Norway',
              'lat': 78.907,
              'lon': 11.889,
              'abbreviation': 'ZEP',
              'alt': 475.0,
              'operation_period': '2018-present'},
             {'name': 'Zugspitze',
              'state': 'Germany',
              'lat': 47.42,
              'lon': 10.98,
              'abbreviation': 'ZSF',
              'alt': 2671.0,
              'operation_period': '2018-present'},
             {'name': 'Helmos Mountain',
              'state': 'Greece',
              'lat': 37.984265,
              'lon': 22.196262,
              'abbreviation': 'HAC',
              'alt': 2314.0,
              'operation_period': '2018-present'},
             {'name': 'Table Mountain',
              'state': 'Colorado',
              'lat': 40.1250,
              'lon': -105.2369995117,
              'abbreviation': 'BOS',
              'alt': 1689,
              'operation_period': '2019-present'},
             {'name': 'Ragged Point',
              'state': '',
              'lat': 13.1650,
              'lon': -59.432,
              'abbreviation': 'RPB',
              'alt': 15,
              'operation_period': '2023-present'},
             {'name': 'University of Granada',
              'state': 'Spain',
              'lat': 37.1640,
              'lon': -3.605,
              'abbreviation': 'UGR',
              'alt': 680.00,
              'operation_period': '2013-present'},
             {'name': 'Varrio',
              'state': 'Finland',
              'lat': 67.7550,
              'lon': 29.6090,
              'abbreviation': 'VAR',
              'alt': 400.00,
              'operation_period': '2018-present'},
             {'name': 'Wright Patterson Air Force Base',
              'state': 'Ohio',
              'lat': 39.82306,
              'lon': -84.04944,
              'abbreviation': 'WPB',
              'alt': 250.8,
              'operation_period': '2022-present'},
             {'name': 'Pond Inlet',
              'state': 'Canada',
              'lat': 72.699,
              'lon': -77.959,
              'abbreviation': 'PON',
              'alt': 55,
              'operation_period': '2018-2020'},
             # {'name': '',
             #  'state': '',
             #  'lat': ,
             #  'lon': ,
             #  'abbreviation': '',
             #  'alt': ,
             #  'operation_period': ''},
            ]

network = _measurement_site.Network(_stations)
network.name = 'aero'