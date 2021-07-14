# -*- coding: utf-8 -*-

import urllib as _ul
import pandas as _pd
import atmPy.general.timeseries as _ts

def get_oceanic_nino_index():
    # get the data from the internet (by parsing)
    url = 'https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php'

    html = _ul.request.urlopen(url).read()
    tabels = _pd.read_html(html)

    i=8
    noi = tabels[i]

    # format the table
    noi_ts = _pd.DataFrame()

    for idx, row in noi.iterrows():
        if row.iloc[0] == 'Year':
            continue
    #     break

        year = row[0]
        values = _pd.DataFrame(row[1:])
        values.index = values.apply(lambda x:  _pd.to_datetime(f'{year}-{x.name:02d}-15'), axis = 1)
        values.index.name = 'datetime'
        values.columns = ['noi']

        noi_ts = noi_ts.append(values.astype(float), sort = True)
        
    noi_ts = _ts.TimeSeries(noi_ts)
    return noi_ts