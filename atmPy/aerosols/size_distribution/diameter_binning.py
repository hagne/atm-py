"""
This module contains function that help with dealing with binning of the diameters
in a sizedistribution.
"""

import numpy as np
import pandas as pd
import numba

def bincenters2binsANDnames(bincenters, round=None):
    """This creates bin edges from bincenters

    TODO
    ----
    One might want to distinguisch between logaritmically spaced and linearly spaced
    size distributions.

    Arguments
    ---------
    bincenters: array.
    round: int
        if rounding is desired in the returned column names

    Returns
    -------
    array: binedges
    array: array of strings that can be used as collumn names in sizedistributions
    """
    if type(bincenters) != np.ndarray:
        raise TypeError('Parameter "bincenters" has to be numpy.ndarray. Given object is %s'%type(bincenters))
    noEnds = (bincenters[:-1]+bincenters[1:])/2.
    firstEdge = bincenters[0] - abs(noEnds[0]-bincenters[0])
    lastEdge = bincenters[-1] + abs(noEnds[-1]-bincenters[-1])
    binedges = np.append(firstEdge,noEnds)
    binedges = np.append(binedges,lastEdge)
    minusses = binedges[:-1].copy().astype(str)
    minusses[:] = ' - '
    a = binedges[:-1]
    b = binedges[1:]
    if isinstance(round, int):
        a = a.round(round)
        b = b.round(round)
    a = a.astype(str)
    b = b.astype(str)
    newColnames = np.core.defchararray.add(a,minusses)
    newColnames = np.core.defchararray.add(newColnames,b)
    return binedges,newColnames


def re_bin(dist, number_of_bins = 50, spaced = 'log', bins = None):
    
    @numba.jit
    def match_bins(index, columns, df_match):
        for idx_row, edo in enumerate(index):
            # old bin is fully inside or outside the particular new bin
    #         pass
    #     return
            for idx_col,edn in enumerate(columns):
    #             pass
                if (edo[0] >= edn[0]) & (edo[1] <= edn[1]):
                    df_match[idx_row, idx_col] = 1
                elif (edo[0] < edn[0]) & (edo[1] <= edn[0]):
                    df_match[idx_row, idx_col] = 0
                elif (edo[0] >= edn[1]) & (edo[1] > edn[1]):
                    df_match[idx_row, idx_col] = 0
            # new bin is partially in new bin... get fraction
                else:
                    bwo = edo[1] - edo[0]
                    if (edn[1] > edo[0]) & (edn[1] < edo[1]):
                        bwr = edn[1] - edo[0]
                        df_match[idx_row, idx_col] = bwr/bwo
                    elif (edn[0] > edo[0]) & (edn[0] < edo[1]):
                    #     print('bla')
                    #     bwo = edo[1] - edo[0]
                        bwr = edo[1] - edn[0]
                        df_match[idx_row, idx_col] = bwr/bwo
                    else:
                        assert(False), 'This should not be possible'
    
    @numba.jit
    def match2data(data, df_match,new_data):
        for idx_row,data_row in enumerate(data):
            for idx_col, match_col in enumerate(df_match.transpose()):
    #             pass
                new_data[idx_row, idx_col] = (data_row * match_col).sum()
    
    dist = dist.convert2numberconcentration()

    # old bins
    bo = dist.bins.astype(np.float32) # float32 is required to avoid those cases of xy.0000000001 which mess up the cases below 1.0 == 1.0000000001 -> False
    # new bins
    if isinstance(bins, type(None)):
        if spaced == 'log':
            bn = np.logspace(np.log10(bo[0]),np.log10(bo[-1]), number_of_bins).astype(np.float32)
        else:
            assert(False), 'only spaced = log available right now'
    else:
        bn = bins
    
    # generate the matching table ... match old bin to new bin

    match_index = np.column_stack((bo[:-1], bo[1:]))
    match_columns = np.column_stack((bn[:-1], bn[1:]))
    match_data = np.zeros((match_index.shape[0], match_columns.shape[0]))
    match_bins(match_index, match_columns, match_data)
    
    # project match table onto old data to get the new data
    data = dist.data.copy()
    new_data = np.zeros((data.shape[0], match_data.shape[1]))
    match2data(data.values, match_data, new_data)
    new_data = pd.DataFrame(new_data, index = data.index)

    dist_new = type(dist)(new_data, bn, 'numberConcentration')
    dist_new.housekeeping = dist.housekeeping
    return dist_new
