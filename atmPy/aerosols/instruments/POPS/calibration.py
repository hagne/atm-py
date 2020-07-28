#from POPS_lib.fileIO import read_Calibration_fromFile,read_Calibration_fromString,save_Calibration
#import fileIO
from scipy import interpolate, optimize
import numpy as np
import pylab as plt
from io import StringIO as io
import pandas as pd
import warnings
from atmPy.aerosols.instruments.POPS import mie


#read_fromFile = fileIO.read_Calibration_fromFile
#read_fromString = fileIO.read_Calibration_fromString

def _msg(txt, save, out_file, verbose):
    if verbose:
        print(txt)
    if save:
        out_file.write(str(txt) + '\n')


def generate_calibration(single_pnt_cali_d=508,
                         single_pnt_cali_ior=1.6,
                         single_pnt_cali_int=1000,
                         noise_level = 12,
                         ior=1.5,
                         dr=[100, 5000],
                         no_pts=5000,
                         no_cal_pts=50,
                         plot=True,
                         raise_error=True,
                         test=False
                         ):
    """
    This function generates a calibration function for the POPS instrument based on its theoretical responds.

    Args:
        single_pnt_cali_d: float [508]
            Diameter of single point calibration in nm.
        single_pnt_cali_ior: float [1.6]
            Refractive index of material used in single point calibration.
        single_pnt_cali_int: float [1000]
            Raw intensity (digitizer bins) measured in single point calibration
        noise_level: int
            The way POPS detectes peak height it is effected by the noise which results in a positive bias for sizing
            close to the lower detection limit ... this corrects for it. If you don't want this set noise_level to None
        ior: float [1.5]
            Refractive index of the anticipated aerosol material.
        dr: array-like [[110, 3400]]
            Diameter range of the calibration. The calibration range will actually be a bit smaller than this, so make
            this range a little bit larger than you want it.
        no_pts: int [600]
            Number of points used in the Mie calculations... quite unimportant value.
        no_cal_pts: int
            Number of points in the generated calibration. Usually a mie curve is not bijective, this number defines how
            much a respondscurve is smoothened. This is merely a starting number. If bejectivity is not given, the number
            is reduced until bijectivity is achieved.
        plot: bool [False]
            If the plotting of the result is desired.
        raise_error: bool [True]
            If an error is raised in case the resulting calibration function is not bijective.
        test: bool [False]
            If True the calibration diameters are returned, so one can check if they are in the desired range.
    Returns:
        Calibration instance
        if plot: (Calibration instance, Axes instance)
        if test: Series instance
    """
    drum = np.array(dr)/1e3
    d, amp = mie.makeMie_diameter(diameterRangeInMikroMeter=drum,
                                 noOfdiameters=no_pts,
                                 IOR=ior)

    df = pd.DataFrame({'d': d, 'amp': amp})
    # a calibration function is created with no_cal_pts of calibration points. no_cal_pts is steadyly decreased until
    # the calibration function is bijective
    valid = False
    while not valid:
        no_cal_pts -= 1
        binedgs = np.logspace(np.log10(df.amp.min()), np.log10(df.amp.max()), no_cal_pts)
        binedgs[0] -= (binedgs[1] - binedgs[
            0]) * 0.01  # this is to ensure the first point is not onthe edge ... for cut function used later

        mie_cal = df.groupby(pd.cut(df.amp, binedgs)).median()
        dfstd = df.groupby(pd.cut(df.amp, binedgs)).mad()

        mie_cal.index = mie_cal.d
        dfstd.index = mie_cal.d
        mie_cal['sigma_d'] = dfstd.d
        mie_cal.index.name = None
        mie_cal.sort_values('amp', axis=0, inplace=True)

        # check for bijectivity
        if not ((mie_cal.d.values[1:] - mie_cal.d.values[:-1]) < 0).sum():
            valid = True

    # final conditioning: um2nm, indexname, drop d
    # mie_cal.drop('d', axis=1, inplace=True)
    mie_cal.index.name = 'd_nm'
    mie_cal.index *= 1e3
    mie_cal.d *= 1e3
    mie_cal.sigma_d *= 1.e3
    cali_inst_pre = Calibration(mie_cal)

    # single point calibration
    ## solve calibration function ot get amp at calibration diameter
    ### first guess
    dt = mie_cal.index[abs(mie_cal.index - single_pnt_cali_d).argmin()]
    at = mie_cal.loc[dt, 'amp']

    # cali_inst_at_single_pnt_calid_d = cali_inst_pre.calibrationFunction(single_pnt_cali_d)

    ### solve
    cali_inst_at_single_pnt_calid_d = optimize.fsolve(lambda x: cali_inst_pre.calibrationFunction(x) - single_pnt_cali_d, at)

    ## scale for ior mismatch
    if ior == single_pnt_cali_ior:
        scale_ioradj = 1
        single_pnt_cali_int_pre = cali_inst_at_single_pnt_calid_d
    else:
        # single_pnt_cali_d = 500
        single_pnt_cali_d *= 1e-3
        dt, mt = mie.makeMie_diameter(diameterRangeInMikroMeter=[single_pnt_cali_d, single_pnt_cali_d + 1e-3],
                                      IOR=single_pnt_cali_ior,
                                      noOfdiameters=2)
        single_pnt_cali_int_pre = mt[0]

        scale_ioradj = cali_inst_at_single_pnt_calid_d / single_pnt_cali_int_pre

    ## scale for instrument calibration
    scale_instrument = single_pnt_cali_int / single_pnt_cali_int_pre

    ## total scale
    scale = scale_ioradj * scale_instrument

    mie_cal.loc[:,'amp'] *= scale

    if not isinstance(noise_level, type(None)):
        mie_cal.loc[:, 'amp'] += noise_level

    cali_inst = Calibration(mie_cal)

    if plot:
        f, a = plt.subplots()
        a.plot(df.d * 1e3, df.amp * scale, label='POPS resp.')
        cali_inst.plot(ax = a)
        # a.plot(ampm.index * 1e3, ampm.values * scale, label='POPS resp. smooth')
        # g, = a.plot(cali.index, cali.values, label='cali')
        # g.set_linestyle('')
        # g.set_marker('x')
        # g.set_markersize(10)
        # g.set_markeredgewidth(2)
        g, = a.plot(single_pnt_cali_d * 1e3, single_pnt_cali_int, label='single ptn cal')
        g.set_linestyle('')
        g.set_marker('x')
        g.set_markersize(10)
        g.set_markeredgewidth(2)
        g.set_label('single pt. cali.')
        # # st.plot(ax = a)
        # a.loglog()
        a.legend()
        # return dft
        # return cali_inst, a

    ##########
    # start of old mie calcultations
    old = False
    if old:
        dr = np.array(dr)/1000

        single_pnt_cali_d *= 1e-3
        # rr = dr / 2 / 1000

        d, amp = mie.makeMie_diameter(noOfdiameters=no_pts,
                                      diameterRangeInMikroMeter = dr,
                                      # radiusRangeInMikroMeter=rr,
                                      IOR=ior)
        ds = pd.Series(amp, d)

        cal_d = pd.Series(index=np.logspace(np.log10(dr[0]), np.log10(dr[1]), no_cal_pts + 2)[1:-1])

        if test:
            return cal_d

        # single point calibration
        if ior == single_pnt_cali_ior:
            ds_spc = ds
        else:
            d, amp = mie.makeMie_diameter(noOfdiameters=no_pts,
                                          diameterRangeInMikroMeter = dr,
                                          # radiusRangeInMikroMeter=rr,
                                          IOR=single_pnt_cali_ior)
            ds_spc = pd.Series(amp, d)

        # rolling, the issue here is that I am rolling over windows of diameter rather then windows of amp, which would be better
        ampm = ds.rolling(int(no_pts / no_cal_pts), center=True).mean()

        #
        cali = ampm.append(cal_d).sort_index().interpolate().reindex(cal_d.index)

        spc_point = ds_spc.append(pd.Series(index=[single_pnt_cali_d])).sort_index().interpolate().reindex(
            [single_pnt_cali_d])  # .values[0]
        scale = single_pnt_cali_int / spc_point.values[0]

        cali *= scale
        cali.index *= 1e3

        cali_inst_pre = pd.DataFrame(cali, columns=['amp'])
        cali_inst_pre['d'] = cali_inst_pre.index
        cali_inst_pre = Calibration(cali_inst_pre)

        if raise_error:
            ct = cali.values
            if (ct[1:] - ct[:-1]).min() < 0:
                raise ValueError(
                    'Clibration function is not bijective. usually decreasing the number of calibration points will help!')

            cal_fkt_test = cali_inst_pre.calibrationFunction(cali_inst_pre.data.amp.values)
            if not np.all(~np.isnan(cal_fkt_test)):
                raise ValueError(
                    'Clibration function is not bijective. usually decreasing the number of calibration points will help!')

        if plot:
            f, a = plt.subplots()
            a.plot(ds.index * 1e3, ds.values * scale, label='POPS resp.')
            a.plot(ampm.index * 1e3, ampm.values * scale, label='POPS resp. smooth')
            g, = a.plot(cali.index, cali.values, label='cali')
            g.set_linestyle('')
            g.set_marker('x')
            g.set_markersize(10)
            g.set_markeredgewidth(2)
            g, = a.plot(single_pnt_cali_d * 1e3, single_pnt_cali_int, label='single ptn cal')
            g.set_linestyle('')
            g.set_marker('o')
            g.set_markersize(10)
            g.set_markeredgewidth(2)
            # st.plot(ax = a)
            a.loglog()
            a.legend()
            return cali_inst_pre, a
    ###### end of old 2

    return cali_inst

def get_interface_bins(fname, n_bins, imin=1.4, imax=4.8, save=False, verbose = True):
    """Prints the bins assosiated with what is seen on the POPS user interface and the serial output, respectively.

    Parameters
    ----------
    fname: string or calibration instance
        name of file containing a calibration or a calibration instance it self
    n_bins: int
        number of bins
    imin: float [1.4], optional
        log10 of the minimum value considered (digitizer bins)
    imax: float [4.8], optional
        log10 of the maximum value considered (digitizer bins)
    save: bool or string.
        If result is saved into file given by string.


    Returns
    -------
    matplotlib axes instance
    pandas DataFrame instance
    """
    if isinstance(fname, str):
        cal = read_csv(fname)
    else:
        cal = fname

    bin_ed = np.linspace(imin, imax, n_bins + 1)
    bin_center_log = 10 ** ((bin_ed[:-1] + bin_ed[1:]) / 2.)
    bin_center_lin = ((10 ** bin_ed[:-1] + 10 ** bin_ed[1:]) / 2.)
    bin_ed = 10 ** bin_ed
    bin_ed_cal = cal.calibrationFunction(bin_ed)
    bin_center_lin_cal = cal.calibrationFunction(bin_center_lin)
    bin_center_log_cal = cal.calibrationFunction(bin_center_log)
    if save:
        save_file = open(save, 'w')
    else:
        save_file = False

    txt = '''
bin edges (digitizer bins)
--------------------------'''
    _msg(txt, save, save_file, verbose)

    for e, i in enumerate(bin_ed):
        _msg(i, save, save_file, verbose)
    # bin_center_cal = cal.calibrationFunction(bin_center)


    txt = '''
bin centers (digitizer bins)
----------------------------'''
    _msg(txt, save, save_file, verbose)
    for e, i in enumerate(bin_center_lin):
        _msg(i, save, save_file, verbose)

    txt = '''
bin centers of logarithms (digitizer bins)
----------------------------'''
    _msg(txt, save, save_file, verbose)
    for e, i in enumerate(bin_center_log):
        _msg(i, save, save_file, verbose)

    txt = '''

bin edges (nm)
--------------'''
    _msg(txt, save, save_file, verbose)
    for e, i in enumerate(bin_ed_cal):
        _msg(i, save, save_file, verbose)
    # bin_center_cal = cal.calibrationFunction(bin_center)


    txt = '''
bin centers (nm)
----------------'''
    _msg(txt, save, save_file, verbose)
    for e, i in enumerate(bin_center_lin_cal):
        _msg(i, save, save_file, verbose)

    txt = '''
bin centers of logarithms (nm)
----------------'''
    _msg(txt, save, save_file, verbose)
    for e, i in enumerate(bin_center_log_cal):
        _msg(i, save, save_file, verbose)

    out = {}

    df_bin_c = pd.DataFrame(bin_center_lin_cal, index=bin_center_log, columns=['Bin_centers'])
    df_bin_e = pd.DataFrame(bin_ed_cal, index = bin_ed, columns = ['Bin_edges'])
    # a = df.Bin_centers.plot()

    if verbose:
        f, a = plt.subplots()
        d = df_bin_c.Bin_centers.values[1:-1]
        g, = a.plot(np.arange(len(d)) + 2, d)
        g.set_linestyle('')
        g.set_marker('o')
        # g.set_label('')
        a.set_yscale('log')
        a.set_xlim((1, 16))
        a.set_ylim((100, 3000))
        a.set_ylabel('Bin center (nm)')
        a.grid(which='both')
        a.set_xlabel('POPS bin')
        out['axes'] = a
    else:
        out['axes'] = None

    # a.set_title('Bin')


    out['bincenters_v_int'] = df_bin_c
    out['binedges_v_int'] = df_bin_e
    return out


def _string2Dataframe(data, log=True):
    sb = io(data)
    dataFrame = pd.read_csv(sb, sep = ' ', names = ('d','amp')).sort_values('d')
    if log:
        dataFrame.amp = 10 ** dataFrame.amp
    return dataFrame


def read_str(data, log=True):
    '''Read a calibration table from string.

    Arguments
    ---------
    data: string.
        Multiline string with a diameter-intensity pair seperated by space. Diameter in nm, intensity in digitizer bin
        or log_10(digitizer bins).
    log: bool, optional.
        Set True if the intensity values are given in log_10(digitizer bins).

    Example
    -------
    data = """140 88
    150 102
    173 175
    200 295
    233 480
    270 740
    315 880
    365 1130
    420 1350
    490 1930
    570 3050
    660 4200
    770 5100
    890 6300
    1040 8000
    1200 8300
    1400 10000
    1600 11500
    1880 16000
    2180 21000
    2500 28000s
    3000 37000"""
    read_str(data, log = False)
    '''

    dataFrame = _string2Dataframe(data, log=log)
    calibrationInstance = Calibration(dataFrame)
    return calibrationInstance


def read_csv(fname):
    """ most likely found here"""
    calDataFrame = pd.read_csv(fname)
    calibrationInstance = Calibration(calDataFrame)
    return calibrationInstance

def save_Calibration(calibrationInstance, fname):
    """should be saved hier cd ~/data/POPS_calibrations/"""
    calibrationInstance.data.to_csv(fname, index = False)
    return

class Calibration:
    def __init__(self,dataTabel):
        self.data = dataTabel
        self.calibrationFunction = self.get_calibrationFunctionSpline()
        
    def get_interface_bins(self, n_bins, imin=1.4, imax=4.8, save=False, verbose = False):
        out = get_interface_bins(self, n_bins, imin=imin, imax=imax, save=save, verbose = verbose)
        return out

    def save_csv(self,fname):
        save_Calibration(self,fname)
        return
        
    def get_calibrationFunctionSpline(self, fitOrder = 1):# = 1, noOfPts = 500, plot = False):
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
        if (self.data.amp.values[1:]-self.data.amp.values[:-1]).min() < 0:
            warnings.warn('The data represent a non injective function! This will not work. plot the calibration to see what I meen')

        # OLD 1
        #
        # sf = UnivariateSpline(self.data.d.values, self.data.amp.values, s=fitOrder)
        # d = np.logspace(np.log10(self.data.d.values.min()), np.log10(self.data.d.values.max()), 500)
        # amp = sf(d)
        #
        # ##### second step
        # cal_function = UnivariateSpline(amp, d, s=fitOrder)

        # OLD 2
        # not sure why I did it this way, I think I wanted to avoid overswingers ... not sure why I did not just do a
        # quatratic interpolation ... below

        # sf = UnivariateSpline(np.log10(self.data.d.values), np.log10(self.data.amp.values), s=0)
        # d = np.linspace(np.log10(self.data.d.values.min()), np.log10(self.data.d.values.max()), 500)
        # amp = sf(d)
        #
        # # us = UnivariateSpline(np.log10(self.data.amp), np.log10(self.data.d), s=0)
        # us = UnivariateSpline(amp, d, s=0)
        # cal_function = lambda amp: 10**us(np.log10(amp))

        # Current
        # simply interpolate (quadratic)
        cal_function = interpolate.interp1d(self.data.amp, self.data.d,
                                            # kind='quadratic',
                                            kind='cubic'
                                            )
        return cal_function
        
    def plot(self, ax=None, nofpts = 500):
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
        amp = np.logspace(np.log10(self.data.amp.min()), np.log10(self.data.amp.max()), nofpts)

        # rounding in logspace can get us out of the domain of cal_function. Therefore the first and last values are
        # rounded up and down, respectively
        amp[0] = round(amp[0] + 0.005, 2) #equivalent with floor with a given decimal
        amp[-1] = np.floor(amp[-1])

        d = cal_function(amp)

        if type(ax).__name__ == 'AxesSubplot':
            a = ax
            f = a.get_figure()
        else:
            f,a = plt.subplots()
        
        cal_data, = a.plot(self.data.d,  self.data.amp, 'o',label = 'data',)
        cal_func, = a.plot(d,amp, label = 'function')

        if 'sigma_d' in self.data.columns:
            a.errorbar(self.data.d, self.data.amp, xerr=self.data.sigma_d,
                   # yerr=dfstd.amp,
                   # zorder = 10,
                   ls ='')
        
        a.loglog()
    
        a.set_xlim(0.9*self.data.d.min(), 1.1*self.data.d.max())
        a.set_xlabel('Diameter (nm)')#($\mu$m)')
    
        a.set_ylim(0.9*self.data.amp.min(), 1.1*self.data.amp.max()) 
        a.set_ylabel('Amplitude (digitizer bins)')
    
        a.set_title('Calibration curve')
        a.legend(loc = 2)
        return f,a,cal_data, cal_func