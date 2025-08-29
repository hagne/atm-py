# -*- coding: utf-8 -*-
import urllib as _ul
import pandas as _pd
from atmPy.general import measurement_site as _ms

def get_all_sites(url = "https://www.esrl.noaa.gov/gmd/dv/site/", return_type = 'network'):
    '''
    

    Parameters
    ----------
    url : TYPE, optional
        DESCRIPTION. The default is "https://www.esrl.noaa.gov/gmd/dv/site/".
    return_type : TYPE, optional
        Choose from ['DataFrame', 'network']. The default is 'network'.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return_type_options = ['DataFrame', 'network']
    if return_type not in return_type_options:
        raise ValueError(f'Unknown return_type: {return_type}. Choose from {return_type_options}')
    
    
    html = _ul.request.urlopen(url).read()
    
    sites = _pd.read_html(html)[0]
    
    #remove discontinued
    sites = sites[~sites.apply(lambda row: '*' in row.Code, axis = 1)].copy()
    
    sites['name'] = sites.apply(lambda row: row.Name.split(',')[0], axis=1)
    
    def gstate(row):
        try: 
            out = row.Name.split(',')[1] 
        except IndexError: 
            out = None
        return out
    sites['state'] = sites.apply(gstate, axis=1)
    
    sites.rename({'Longitude': 'lon',
                  'Latitude': 'lat', 
                  'Elevation (meters)': 'alt',
                  'Code': 'abbreviation',
                  'Country': 'country'}, axis=1, inplace = True)
    
    if return_type == 'DataFrame':
        return sites
    
    elif return_type == 'network':
        site_dict_list = []
        for idx, si in sites.iterrows():
            site_dict_list.append(si.to_dict())
        
        gml_sites = _ms.Network(site_dict_list)
        return gml_sites
    else:
        assert(False), 'not possible (35454222)'