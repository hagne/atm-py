#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:24:16 2021

@author: hagen
"""
import pandas as pd
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.cpd.cpd3_lab as cpd3_lab


def read_file(fn):
#     out = {}
    df = pd.read_csv(fn)
    cpd3_lab.fixtime(df)
    return df