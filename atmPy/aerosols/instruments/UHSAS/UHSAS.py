# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 11:43:10 2014

@author: htelg
"""

import datetime
import warnings
from io import StringIO as io

import numpy as np
import pandas as pd
import pylab as plt
from scipy.interpolate import UnivariateSpline

from atmPy.general import timeseries
from atmPy.aerosols.size_distribution import sizedistribution


def read_csv(fname, norm2time = True, norm2flow = True):
    uhsas_file_types = ['.xls']
    first = True
    if type(fname).__name__ == 'list':
        for file in fname:
            for i in uhsas_file_types:
                if i in file:
                    right_file_format = True
                else:
                    right_file_format = False
            if right_file_format:
                sdt, hkt= _read_csv(file, norm2time = norm2time, norm2flow = norm2flow)

                if first:
                    sd = sdt.copy()
                    hk = hkt.copy()
                    first = False
                else:
                    if not np.array_equal(sd.bincenters, sdt.bincenters):
                        txt = 'the bincenters changed between files! No good!'
                        raise ValueError(txt)

                    sd.data = pd.concat((sd.data,sdt.data))
                    hk.data = pd.concat((hk.data,hkt.data))

        if first:
            txt = """Either the prvided list of names is empty, the files are empty, or none of the file names end on
the required ending (*.xls)"""
            raise ValueError(txt)
    else:
        sd, hk= _read_csv(fname, norm2time = norm2time, norm2flow = norm2flow)
    return sd, hk


def _read_csv(fname, norm2time = True, norm2flow = True):
    uhsas = _readFromFakeXLS(fname)
#     return uhsas
    sd,hk = _separate_sizedist_and_housekeep(uhsas, norm2time = norm2time, norm2flow = norm2flow)
    hk = timeseries.TimeSeries(hk)
#     return size_distribution,hk
    bins = _get_bins(sd)
#     return bins
    dist = sizedistribution.SizeDist_TS(sd, bins, "numberConcentration")
    return dist, hk

def _readFromFakeXLS(fname):
    """reads and shapes a XLS file produced by the uhsas instruments"""
    fr = pd.read_csv(fname, sep='\t')
    newcolname = [fr.columns[e] + ' ' + str(fr.values[0][e]) for e, i in enumerate(fr.columns)]
    fr.columns = newcolname
    fr = fr.drop(fr.index[0])
    bla = pd.Series(fr['Date -'].values + ' ' + fr['Time -'].values)
#     return bla
    try:
        fr.index = bla.map(lambda x: datetime.datetime.strptime(x, '%m/%d/%Y %H:%M:%S.%f'))
    except ValueError:
        fr.index = bla.map(lambda x: datetime.datetime.strptime(x, '%m/%d/%Y %I:%M:%S.%f %p'))
    fr = fr.drop(['Date -', 'Time -'], axis=1)
    return fr

def _separate_sizedist_and_housekeep(uhsas, norm2time = True, norm2flow = True):
    """Beside separating size distribution and housekeeping this
    function also converts the data to a numberconcentration (#/cc)

    Parameters
    ----------
    uhsas: pandas.DataFrame"""

#     size_distribution = uhsas.copy()
#     hk = uhsas.copy()
# #     return size_distribution,hk

    first = False
    for e,col in enumerate(uhsas.columns):
        cola = col.split(' ')[0]

        try:
            float(cola)
            float(col.split(' ')[1])
        except ValueError:
            continue
        else:
            last = e
            if not first:
                first = e



    # k = size_distribution.keys()
    # where = np.argwhere(k == 'Valve 0=bypass') + 1

    hk = uhsas.iloc[:,:first]
    sd = uhsas.iloc[:,first:last+1]
    # khk = k[: first]
    # size_distribution = size_distribution.drop(khk, axis=1)
    # hsd = k[where:]
    # hk = hk.drop(hsd, axis=1)
#     return size_distribution,hk
    hk['Sample sccm'] = hk['Sample sccm'].astype(float)

    hk['Accum. Secs'] = hk['Accum. Secs'].astype(float)

    # normalize to time and flow
    if norm2time:
        sd = sd.mul(1 / hk['Accum. Secs'], axis = 0 )
    if norm2flow:
        sd = sd.mul(60./hk['Sample sccm'], axis = 0 )
    return sd,hk

def _get_bins(frame, log=False):
    """
    get the bins from the column labels of the size distribution DataFrame.
    """
    frame = frame.copy()

    bins = np.zeros(frame.keys().shape[0]+1)
    for e, i in enumerate(frame.keys()):
        bin_s, bin_e = i.split(' ')
        bin_s = float(bin_s)
        bin_e = float(bin_e)
        bins[e] = bin_s
        bins[e+1] = bin_e
    return bins #binCenters



def _string2dataframe(data):
    sb = io(data)
    dataFrame = pd.read_csv(sb,
#                             sep=' ',
                            names=('d', 'bin_no')
                           ).sort('d')
    return dataFrame

def read_calibration_fromString(data):
    '''
    unit of diameter must be nm
    e.g.:
data = """120., 19.5
130., 22.5
140., 25
150., 27.6
173., 33.
200., 38.
233., 43.4
270., 47.5
315., 53.
365., 58.
420., 62.5
490., 67.
570., 71.
660., 75.
770., 78.
890., 79.
1040., 84."""
    '''

    dataFrame = _string2dataframe(data)
#     return dataFrame
    calibrationInstance = calibration(dataFrame)
    return calibrationInstance

class calibration:
    def __init__(self,dataTabel):
        self.data = dataTabel
        self.calibrationFunction = self.get_calibrationFunctionSpline()

    def save_csv(self,fname):
#         save_Calibration(self,test_data_folder)
        self.data.to_csv(fname, index = False)
        return

    def get_calibrationFunctionSpline(self, fitOrder=1):
        """
        Performes a spline fit/smoothening (scipy.interpolate.UnivariateSpline) of d over amp (yes this way not the other way around).

        Returns (generates): creates a function self.spline which can later be used to calculate d from amp

        Optional Parameters:
        \t s: int - oder of the spline function
        \t noOfPts: int - length of generated graph
        \t plot: boolean - if result is supposed to be plotted
        """

        # The following two step method is necessary to get a smooth curve.
        #When I only do the second step on the cal_curve I get some wired whiggles
        ##### First Step
        if (self.data.bin_no.values[1:]-self.data.bin_no.values[:-1]).min() < 0:
            warnings.warn('The data represent a non injective function! This will not work. plot the calibration to see what I meen')

        sf = UnivariateSpline(self.data.d.values, self.data.bin_no.values, s=fitOrder)
        d = np.logspace(np.log10(self.data.d.values.min()), np.log10(self.data.d.values.max()), 500)
        bin_no = sf(d)

        # second step
        cal_function = UnivariateSpline(bin_no, d, s=fitOrder)
        return cal_function

    def plot_calibration(self):
        """Plots the calibration function and data
        Arguments
        ------------
            cal: calibration instance

        Returns
        ------------
            figure
            axes
            calibration data graph
            calibration function graph
        """
        cal_function = self.calibrationFunction
        bin_no = np.logspace(np.log10(self.data.bin_no.min()), np.log10(self.data.bin_no.max()), 500)
        d = cal_function(bin_no)

        f, a = plt.subplots()

        cal_data, = a.plot(self.data.d,  self.data.bin_no, 'o', label='data',)
        cal_func, = a.plot(d, bin_no, label='function')

        a.loglog()

        a.set_xlim(0.9*self.data.d.min(), 1.1*self.data.d.max())
        a.set_xlabel('Diameter (nm)')

        a.set_ylim(0.9*self.data.bin_no.min(), 1.1*self.data.bin_no.max())
        a.set_ylabel('bin number')

        a.set_title('Calibration curve')
        a.legend(loc = 2)
        return f, a, cal_data, cal_func

    def apply_on(self, dist, limit_to_cal_range = True):
        dist_t = dist.copy()
        bins_no = np.arange(dist_t.bins.shape[0])
        cal_f = self.get_calibrationFunctionSpline()

        new_d = cal_f(bins_no)
        df = pd.DataFrame(np.array([bins_no, new_d]).transpose(), columns = ['bin_no','d'])

        dist_t.bins = new_d

        start_d = self.data.d.iloc[0]
        end_d = self.data.d.iloc[-1]

        if limit_to_cal_range:
            dist_t = dist_t.zoom_diameter(start = start_d, end=end_d)
        return dist_t