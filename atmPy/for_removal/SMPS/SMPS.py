# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 16:39:10 2015

@author: htelg
"""

import sys
sys.path.append('/Users/htelg/projecte/POPS/prog/')
from POPS_lib import sizedistribution
import pandas as pd
import numpy as np



    
def bincenters2BinStuff(bincenters):
    if type(bincenters) != np.ndarray:
        raise TypeError('Parameter "bincenters" has to be numpy.ndarray. Given object is %s'%type(bincenters))
    noEnds = (bincenters[:-1]+bincenters[1:])/2.
    firstEdge = bincenters[0] - (noEnds[0]-bincenters[0])
    lastEdge = bincenters[-1] + (noEnds[-1]-bincenters[-1])
    binedges = np.append(firstEdge,noEnds)
    binedges = np.append(binedges,lastEdge)
    minusses = binedges[:-1].copy().astype(str)
    minusses[:] = ' - '
    a = binedges[:-1].astype(str) 
    b = binedges[1:].astype(str)
    newColnames = np.core.defchararray.add(a,minusses)
    newColnames = np.core.defchararray.add(newColnames,b)
    return binedges,newColnames

def load_PMEL_SMPS(fname):
    """returns an aerosolsizedistributio instance from the POPS_lib libary"""

    tab = pd.read_csv(fname,sep = '\t',header = 18)
    tab.index = pd.to_datetime(tab.Date+' '+tab['Start Time'])
    reducedTab = tab.iloc[:,8:-26]
    col = reducedTab.columns
    bincenters = col.values.astype(float)
    binedges,newCol = bincenters2BinStuff(bincenters)
    reducedTab.columns = newCol
    dist = sizedistribution.aerosolSizeDistribution(reducedTab,binedges, 'dNdlogDp', bincenters = bincenters)
    return dist
    
#def load_PMEL_SMPS(test_data_folder):
#    """returns an aerosolsizedistributio instance from the POPS_lib libary"""
#
#    tab = pd.read_csv(test_data_folder,sep = '\t',header = 18)
#    tab.index = pd.to_datetime(tab.Date+' '+tab['Start Time'])
#    reducedTab = tab.iloc[:,8:-26]
#    col = reducedTab.columns
#    bincenters = col.values.astype(float)
#    noEnds = (bincenters[:-1]+bincenters[1:])/2.
#    firstEdge = bincenters[0] - (noEnds[0]-bincenters[0])
#    lastEdge = bincenters[-1] + (noEnds[-1]-bincenters[-1])
#    binedges = np.append(firstEdge,noEnds)
#    binedges = np.append(binedges,lastEdge)
#    minusses = binedges[:-1].copy().astype(str)
#    minusses[:] = ' - '
#    a = binedges[:-1].astype(str) 
#    b = binedges[1:].astype(str)
#    newCol = np.core.defchararray.add(a,minusses)
#    newCol = np.core.defchararray.add(newCol,b)
#    reducedTab.columns = newCol
#    dist = sizedistribution.aerosolSizeDistribution(reducedTab,binedges, 'dNdlogDp', bincenters = bincenters)
#    return dist