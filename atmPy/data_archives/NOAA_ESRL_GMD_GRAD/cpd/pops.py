#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:21:13 2021

@author: hagen
"""

import pandas as pd
import atmPy.aerosols.size_distribution.sizedistribution as sd
import atmPy.aerosols.size_distribution.diameter_binning as db


header_dict = [  #{'cpd3': 'DateTimeUTC',   'verbose': 'Date String (YYYY-MM-DD hh:mm:ss) UTC'              ,'mylabel': 'datetime'},
                 {'cpd3': 'F1_N21',        'verbose': 'Instrument flags'                                   ,'mylabel': 'flag_instrument'},
                 {'cpd3': 'N_N21',         'verbose': 'Total concentration (cm⁻³)'                         ,'mylabel': 'total_concentration'},
                 {'cpd3': 'P_N21',         'verbose': 'Board pressure sensor (hPa)'                        ,'mylabel': 'ptu_pressure'},
                 {'cpd3': 'I_N21',         'verbose': 'Baseline value'                                     ,'mylabel': 'baseline'},
                 {'cpd3': 'Ig_N21',        'verbose': 'Baseline standard deviation'                        ,'mylabel': 'baseline_std'},
                 {'cpd3': 'Q_N21',         'verbose': 'Sample flow (lpm)'                                  ,'mylabel': 'flow'},
                 {'cpd3': 'T1_N21',        'verbose': 'Temperature of pressure sensor (°C)'                ,'mylabel': 'ptu_temp'},
                 {'cpd3': 'T2_N21',        'verbose': 'Laser temperature (°C)'                             ,'mylabel': 'laser_temp'},
                 {'cpd3': 'T3_N21',        'verbose': 'Internal temperature (°C)'                          ,'mylabel': 'temp_extra'},
                 {'cpd3': 'ZLASERMON_N21', 'verbose': 'Laser monitor'                                      ,'mylabel': 'laser_monitor'},
                 {'cpd3': 'ZMEANWIDTH_N21','verbose': 'Mean width of the detection peak in sampling cycles','mylabel': 'peak_width'},
                 {'cpd3': 'ZPUMPFB_N21',   'verbose': 'Pump feedback'                                      ,'mylabel': 'pump_feedback'},
              ]



def read_file(fn):
    out = {}
    df = pd.read_csv(fn)
    df.index = pd.to_datetime(df.DateTimeUTC)
    df.drop('DateTimeUTC', axis = 1, inplace = True)
    # df.shape

    dist = df.loc[:,[i for i in df.columns if i[:2]=='Nb']].copy().astype(float)
    dist.columns = df.loc[:,[i for i in df.columns if i[:2]=='Ns']].iloc[0].astype(float) * 1000

#     dist.index = pd.to_datetime(df.DateTimeUTC)



    dist = sd.SizeDist_TS(dist,db.bincenters2binsANDnames(dist.columns.values)[0] , 'numberConcentration')
    dist = dist.convert2dNdlogDp()
    out['size_distribution'] = dist
    if 1:
        rest = df.drop([i for i in df.columns if i[:2]=='Nb'], axis = 1)
        rest = rest.drop([i for i in rest.columns if i[:2]=='Ns'], axis = 1)
        rest = rest.rename({h['cpd3']: h['mylabel'] for h in header_dict}, axis = 1)
        out['rest'] = rest
    return out




    


