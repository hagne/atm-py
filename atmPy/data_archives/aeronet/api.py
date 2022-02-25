#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 09:58:35 2022

@author: hagen
"""

import atmPy.data_archives.aeronet.file_io.aod as fileioaod
import atmPy.data_archives.aeronet.file_io.inversion as fileioinv


def read_file(path, verbose = False):
    
    header = ''
    with open(path) as rein:
        for i in range(3):
            header += rein.readline()
    
    if 'Direct Sun' in header:
        out = fileioaod.read_file(path)
        if verbose:
            print('read file as direct sun retrieval')
    elif 'Almucantar' in header:
        out = fileioinv.read_file(path)
        if verbose:
            print('read file as almucantar retrieval')
    return out