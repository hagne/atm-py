#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:32:08 2021

@author: hagen
"""
import pandas as pd
import atmPy.aerosols.size_distribution.sizedistribution as sd
import atmPy.aerosols.size_distribution.diameter_binning as db


def read_file(fn):
    out = {}
    df = pd.read_csv(fn)
    df.index = pd.to_datetime(df.DateTimeUTC)
    df.drop('DateTimeUTC', axis = 1, inplace = True)
    # df.shape

    dist = df.loc[:,[i for i in df.columns if i[:2]=='Nn']].copy().astype(float)
    dist.columns = df.loc[:,[i for i in df.columns if i[:2]=='Ns']].iloc[0].astype(float) * 1000

    #     dist.index = pd.to_datetime(df.DateTimeUTC)



    dist = sd.SizeDist_TS(dist,db.bincenters2binsANDnames(dist.columns.values)[0] , 'dNdlogDp')
    dist = dist.convert2dNdlogDp()
    out['size_distribution'] = dist

    rest = df.drop([i for i in df.columns if i[:2]=='Nn'], axis = 1)
    rest = rest.drop([i for i in df.columns if i[:2]=='Ns'], axis = 1)
#     rest = rest.rename({h['cpd3']: h['mylabel'] for h in header_dict}, axis = 1)
    out['rest'] = rest
    return out