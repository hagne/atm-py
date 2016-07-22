# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 21:23:22 2015

@author: htelg
"""

import warnings

import numpy as np
import pandas as pd
import pylab as plt
from scipy import stats

from atmPy.tools import array_tools, plt_tools
from atmPy.tools import math_linear_algebra as mla

# from scipy import integrate
from atmPy.radiation import solar
from atmPy.radiation.rayleigh import bucholtz_rayleigh as bray
from atmPy.general import atmosphere_standards as atmstd, timeseries
from scipy import signal
from atmPy.tools import time_tools
from copy import deepcopy
from matplotlib import colors, cm
import os

year = '2015'
miniSASP_channels = [550.4, 460.3, 671.2, 860.7]


def read_csv(fname, version = 'current', verbose=False):
    """Creates a single ULR instance from one file or a list of files.

    Arguments
    ---------
    fname: string or list
    version: str
        0.1: files till 2016-07-18 ... this includes Svalbard data
        current: files since 2016-07-18
    """
    if type(fname).__name__ == 'list':
        first = True
        for file in fname:
            if os.path.split(file)[-1][0] != 'r':
                continue
            if verbose:
                print(file)
            ulrt = miniSASP(file, verbose=verbose)
            if first:
                ulr = ulrt
                first = False
            else:
                ulr.data = pd.concat((ulr.data, ulrt.data))
    else:
        ulr = miniSASP(fname, verbose=verbose)
    ulr.data = ulr.data.sort_index()
    return ulr


# todo: see warning below
def _recover_Negatives(series, verbose=True):
    series = series.values
    where = np.where(series>2**16)
    series[where] = (series[where]*2**16)/2**16
    if where[0].shape[0] > 0:
        if verbose:
            warnings.warn("""This has to be checked!!! Dont know if i implemented it correctly! Arduino negatives become very large positive (unsigned longs)
           ;  recover using deliberate overflow""")


def _extrapolate(x, y):
    """This is a very simple extrapolation. 
    Takes: two series of the same pandas DataTable. x is most likely the index of the DataTable
    assumption: 
        - x is are the sampling points while y is the dataset which is incomplete and needs to be extrapolated
        - the relation ship is very close to linear
    proceedure:
        - takes the fist two points of y and performs a linear fit. This fit is then used to calculate y at the very first
          x value
        - similar at the end just with the last two points.
    returns: nothing. everthing happens inplace
    """
    
    xAtYnotNan = x.values[~np.isnan(y.values)][:2]
    YnotNan = y.values[~np.isnan(y.values)][:2]
    slope, intercept, r_value, p_value, slope_std_error = stats.linregress(xAtYnotNan,YnotNan)

    fkt = lambda x: intercept + (slope * x)
    y.values[0] = fkt(x.values[0])

    xAtYnotNan = x.values[~np.isnan(y.values)][-2:]
    YnotNan = y.values[~np.isnan(y.values)][-2:]
    slope, intercept, r_value, p_value, slope_std_error = stats.linregress(xAtYnotNan,YnotNan)

    fkt = lambda x: intercept + (slope * x)
    y.values[-1] = fkt(x.values[-1])
    
    return


def simulate_from_size_dist_LS(dist_LS, airmassfct=True, rotations=2, sun_azimuth=True, pressure=True, temp=True):
    """Simulates AOD and skybrightness from sizeditribution

    Arguments
    ---------
    dist_LS: sizedistribution_LS.
        A sizedistribution layer series
    miniSASP: MiniSASP instance.
        Can be the length of the enitire fight... will be cut accordingly.
    airmassfct: bool.
        if the slant angle is considered.
        False: result is OD
    sun_azimuth: bool.
        if the sun_azimuth is considered. When False the sun is at always at 0 degrees.
        When True suttle changes in the position of the sun can be seen as shifts in the features.

    Retruns:
    --------
    tuple of dicts containing:
        AODs
        skyBrs
    """
    optPs = {}
    for wl in miniSASP_channels:
        optP = dist_LS.calculate_optical_properties(wl, 1.455)
        optPs[wl] = optP

    skyBrs = {}
    aods = {}

    for k in optPs.keys():
        sky, aod = simulate_from_size_dist_opt(optPs[k], airmassfct=airmassfct, rotations=rotations,
                                               sun_azimuth=sun_azimuth, pressure=pressure, temp=temp)
        skyBrs[k] = sky
        aods[k] = aod

    return aods, skyBrs


# ToDo: include actual TEMP and pressure
def simulate_from_size_dist_opt(opt_prop, airmassfct=True, rotations=2, sun_azimuth=True, pressure=True, temp=True):
    """ Simulates miniSASP signal from a size distribution layer series (in particular from the optical property class
    derived from the layer series.
    The simulation calculates the position of the sun at the instruments position during the experiment. Slant angles
    are considert. The atmosphere above the top layer is unkonwn, therefore the measurement from the miniSASP at the top
    most layer should be added to all results.

    Note
    ----
    Temperature and pressure are currently not considered in the underlying Rayleigh calculations. Instead, an
    international standard atmosphere was used.

    Arguments
    ---------
    OpticalProperties class which was created from a layer series (dist_LS) using the sizedistribution module.
    airmassfct: bool, optional.
        If True, results will be corrected for the airmassfactor (slant angle only)
    rotations: int.
        Number of rotations of the mSASP to be simulated.
    pressure: bool or array-like.
        If True the opt_prop.paretn_timeseries.Barometric_pressure timeseries.
        If False standard atmosphere is used.
        If array-like the this array is used.
    temp: bool or array-like.
        If True the opt_prop.paretn_timeseries.Temperature timeseries.
        If False standard atmosphere is used.
        If array-like the this array is used.

    Returns
    -------
    dict:
        containing three (aerosol, rayleigh, sum) pandas DataFrames each with the sky brightness as a function of mSASPS
        azimuth angle.
    pandas DataFrame:
        AOD as a function of elevaton"""

    time_series = opt_prop.parent_dist_LS.parent_timeseries
    dist_ls = opt_prop.parent_dist_LS
    layerthickness = np.apply_along_axis(lambda line: line[1] - line[0], 1, dist_ls.layerbounderies)
    time_series = solar.get_sun_position_TS(time_series)
    where = array_tools.find_closest(time_series.data.Altitude.values, dist_ls.layercenters)
    alts = time_series.data.Altitude.values[where]
    solar_elev = time_series.data.Solar_position_elevation.values[where]
    solar_az = time_series.data.Solar_position_azimuth.values[where]
    # time = time_series.data.index[where]

    what_mSASP_sees_aerosols = pd.DataFrame()
    what_mSASP_sees_AOD_aerosols = np.zeros(alts.shape)

    for altitude in range(dist_ls.layercenters.shape[0]):
        # get the sun position at the time when the plane was at the particular altitude,
        sol_el = solar_elev[altitude]
        sol_az = solar_az[altitude]


        # angles between mSASP positions and sun. This is used to pick the angle in the phase functions
        if sun_azimuth:
            sun_azimuth = sol_az
        else:
            sun_azimuth = 0
        mSASP2Sunangles = angle_MSASP_sun(sol_el,
                                          sun_azimuth=sun_azimuth,
                                          no_angles=int(opt_prop.angular_scatt_func.shape[0] / 2) * rotations,
                                          # pretty arbitrary number ... this is just to get a reasonal number of angles
                                          no_rotations=rotations)

        # pick relevant angles in phase function for each layer, this includes selecting the relavant layers (selected altitude to top).
        closest_phase2sun_azi = array_tools.find_closest(opt_prop.angular_scatt_func.index.values,
                                                         mSASP2Sunangles.mSASP_sun_angle.values)
        # minimize so values are calculated only once
        closest_phase2sun_azi = np.unique(closest_phase2sun_azi)
        phase_fct_rel = opt_prop.angular_scatt_func.iloc[closest_phase2sun_azi, altitude:]
        # Integrate ofer selected intensities in phase function along vertical line (from selected height to top)
        # x = phase_fct_rel.columns.values
        # do_integ = lambda y: integrate.simps(y, x)
        # phase_fct_rel_integ = pd.DataFrame(phase_fct_rel.apply(do_integ, axis=1),
        #                                    columns=[alts[altitude]],
        #                                    # columns=[dist_ls.layercenters[altitude]]
        #                                    )  # these are the integrated intensities of scattered light into the relavant angles. Integration is from current (arbitrary) to top layer
        # print(phase_fct_rel.shape, layerthickness[altitude:].shape)
        phth = phase_fct_rel * layerthickness[altitude:]
        phase_fct_rel_integ = pd.DataFrame(phth.apply(np.sum, 1))
        # return phase_fct_rel, phase_fct_rel_integ

        if airmassfct:
            slant_adjust = 1. / np.sin(solar_elev[altitude])
        else:
            slant_adjust = 1.
        # similar to above this selects the different angels of mSASP to the sun. However, it keeps all of them (no unique)
        closest_phase2sun_azi = array_tools.find_closest(phase_fct_rel_integ.index.values,
                                                         mSASP2Sunangles.mSASP_sun_angle.values)

        what_mSASP_sees_aerosols[dist_ls.layercenters[altitude]] = pd.Series(
            phase_fct_rel_integ.iloc[closest_phase2sun_azi].values.transpose()[0] * slant_adjust)
        # what_mSASP_sees_aerosols[dist_ls.layercenters[altitude]] = pd.Series(
        #     phase_fct_rel_integ.iloc[closest_phase2sun_azi].values.transpose()[0] * slant_adjust)

        # what_mSASP_sees_AOD_aerosols[altitude] = opt_prop.data_orig['AOD_layer'][altitude:].sum().values[0] * slant_adjust
        what_mSASP_sees_AOD_aerosols[altitude] = opt_prop.data_orig['AOD_cum'].values[altitude][0] * slant_adjust
    what_mSASP_sees_aerosols.index = mSASP2Sunangles.index
    # what_mSASP_sees_AOD_aerosols = pd.DataFrame(what_mSASP_sees_AOD_aerosols, index = alts, columns = ['AOD_aerosols'])
    what_mSASP_sees_AOD = pd.DataFrame(what_mSASP_sees_AOD_aerosols, columns=['aerosol'])
    what_mSASP_sees_sky = {'aerosol': what_mSASP_sees_aerosols}

    what_mSASP_sees_rayleigh, what_mSASP_sees_AOD_rayleigh = simulate_from_rayleigh(time_series,
                                                                                    dist_ls.layerbounderies,
                                                                                    pressure,
                                                                                    temp,
                                                                                    opt_prop.wavelength,
                                                                                    what_mSASP_sees_aerosols.shape[0],
                                                                                    rotations,
                                                                                    airmassfct,
                                                                                    sun_azimuth)
    what_mSASP_sees_rayleigh.columns = dist_ls.layercenters
    what_mSASP_sees_sky['rayleigh'] = what_mSASP_sees_rayleigh

    what_mSASP_sees_sum = what_mSASP_sees_aerosols + what_mSASP_sees_rayleigh
    what_mSASP_sees_sky['sum'] = what_mSASP_sees_sum
    # what_mSASP_sees_sky['aerosols'] = what_mSASP_sees_aerosols


    what_mSASP_sees_AOD_sum = what_mSASP_sees_AOD_aerosols + what_mSASP_sees_AOD_rayleigh.values.transpose()[0]
    # what_mSASP_sees_AOD['aerosols'] = what_mSASP_sees_AOD_aerosols
    what_mSASP_sees_AOD['rayleigh'] = what_mSASP_sees_AOD_rayleigh.values.transpose()[0]
    what_mSASP_sees_AOD['sum'] = what_mSASP_sees_AOD_sum
    # return what_mSASP_sees_AOD_aerosols  , what_mSASP_sees_AOD_rayleigh.values.transpose()[0]

    # what_mSASP_sees_AOD_aerosols = pd.DataFrame(what_mSASP_sees_AOD_aerosols, index = alts, columns = ['AOD'])
    what_mSASP_sees_AOD.index = alts
    return what_mSASP_sees_sky, what_mSASP_sees_AOD


def simulate_from_rayleigh(time_series,
                           layerbounderies,
                           # altitude,
                           # layerbounderies,
                           pressure,
                           temp,
                           wl,
                           no_angles,
                           rotations,
                           airmassfct,
                           sun_azimuth):
    """ Fix this documentation!


    Simulates miniSASP signal from a size distribution layer series
    Arguments
    ---------
    layerbounderies: array-like
    altitude: float or array-like.
        Altitude for which the mSASP signal is simulated for in meters
    pressure: array, bool
        Atmospheric pressure in mbar. If False, value is taken from international standard atmosphere.
    temp: array, bool in K
    wl: wavelength in nm

    no_angles: int
        total number of angles considered. This included the number in multiple roations. most likely this is
        int(opt_prop.angular_scatt_func.shape[0] / 2) * rotations

    rotations: int.
        number of rotations the of the mSASP.

    Returns
    -------
    pandas.DataFrame
        containing the sky brightness as a function of mSASPS azimuth angle"""
    layerbounderies = np.unique(layerbounderies.flatten())
    altitude = (layerbounderies[1:] + layerbounderies[:-1]) / 2.
    time_series = solar.get_sun_position_TS(time_series)
    where = array_tools.find_closest(time_series.data.Altitude.values, altitude)
    solar_elev = time_series.data.Solar_position_elevation.values[where]
    solar_az = time_series.data.Solar_position_azimuth.values[where]
    alts = time_series.data.Altitude.values[where]  # thats over acurate, cal simply use the layerbounderies

    if (type(pressure).__name__ == 'bool') or (type(temp).__name__ == 'bool'):
        if pressure and temp:
            # temp = time_series.data.Temperature
            # pressure = time_series.data.Barometric_pressure_Pa
            lb = pd.DataFrame(index=layerbounderies)
            select_list = ["Temperature", "Altitude", "Pressure_Pa"]

            bla = []
            for i in ["Temperature", "Altitude", "Pressure_Pa"]:
                if i not in time_series.data.columns:
                    bla.append(i)

            if len(bla) != 0:
                txt='The underlying housekeeping data has to have the following attributes for this operation to work: %s'%(["Temperature", "Altitude", "Pressure_Pa"])
                txt+='\nmissing:'
                for i in bla:
                    txt += '\n \t' + i
                # print(txt)
                raise AttributeError(txt)

            hkt = time_series.data.loc[:, select_list]

            hkt.index = hkt.Altitude
            hkt = hkt.sort_index()

            hkt_lb = pd.concat([hkt, lb]).sort_index().interpolate()
            hkt_lb = hkt_lb.groupby(hkt_lb.index).mean().reindex(lb.index)
            temp = hkt_lb.Temperature.values + 273.15
            pressure = hkt_lb.Pressure_Pa.values

        else:
            p, t = atmstd.standard_atmosphere(layerbounderies)
            if type(pressure).__name__ == 'bool':
                if pressure == False:
                    pressure = p
            if type(temp).__name__ == 'bool':
                if temp == False:
                    temp = t
    # print(pressure, temp)
    if (layerbounderies.shape != pressure.shape) or (layerbounderies.shape != temp.shape):
        raise ValueError('altitude, pressure and tmp have to have same shape')

    # time = time_series.data.index[where]

    what_mSASP_sees_rayleigh = pd.DataFrame()
    what_mSASP_sees_AOD_rayleigh = np.zeros(altitude.shape)

    for alt in range(altitude.shape[0]):
        # get the sun position at the time when the plane was at the particular altitude,

        sol_el = solar_elev[alt]
        sol_az = solar_az[alt]
        # print(alts[alt:])

        # return ray_scatt_fct

        # angles between mSASP positions and sun. This is used to pick the angle in the phase functions
        if sun_azimuth:
            sun_azimuth = sol_az
        else:
            sun_azimuth = 0
        mSASP2Sunangles = angle_MSASP_sun(sol_el,
                                          sun_azimuth=sun_azimuth,
                                          no_angles=no_angles,
                                          # pretty arbitrary number ... this is just to get a reasonal number of angles
                                          no_rotations=rotations)

        ray_scatt_fct = bray.rayleigh_angular_scattering_intensity(layerbounderies[alt:], pressure[alt:], temp[alt:],
                                                                   wl, mSASP2Sunangles.values.transpose())
        ray_scatt_fct = pd.DataFrame(ray_scatt_fct, index=mSASP2Sunangles.index)
        # return layerbounderies[alt:], pressure[alt:],temp[alt:], wl, ray_scatt_fct
        if airmassfct:
            slant_adjust = 1. / np.sin(solar_elev[alt])
        else:
            slant_adjust = 1.
        # closest_phase2sun_azi = array_tools.find_closest(ray_scatt_fct.index.values,
        #                                                  mSASP2Sunangles.mSASP_sun_angle.values)
        what_mSASP_sees_rayleigh[alts[alt]] = pd.Series(ray_scatt_fct.values.transpose()[0] * slant_adjust)
        # what_mSASP_sees_rayleigh.index = mSASP2Sunangles.index.values

        what_mSASP_sees_AOD_rayleigh[alt] = bray.rayleigh_optical_depth(layerbounderies[alt:], pressure[alt:],
                                                                        temp[alt:], wl) * slant_adjust
        # return layerbounderies[alt:],pressure[alt:],temp[alt:],wl, what_mSASP_sees_AOD_rayleigh[alt], slant_adjust

    what_mSASP_sees_rayleigh.index = mSASP2Sunangles.index
    what_mSASP_sees_AOD_rayleigh = pd.DataFrame(what_mSASP_sees_AOD_rayleigh, index=altitude, columns=['AOD_ray'])
    return what_mSASP_sees_rayleigh, what_mSASP_sees_AOD_rayleigh


def angle_MSASP_sun(sun_elevation, sun_azimuth=0., no_angles=1000, no_rotations=1):
    """Calculates the angle between sun and mini-SASP orientation for one full rotation

    Arguments
    ---------
    sun_elevation: float
        elevation angle of the sun in radians.
    sun_azimuth: float,otional.
        azimuth angle of the sun in radians.
    no_angles:  int, optional.
        number of angles.
    no_rotations: int, otional.
        number of rotations.

    Returns
    -------
    pandas.DataFrame instance
    """
    sunPos = np.array([1, np.pi / 2. - sun_elevation, sun_azimuth])  # just an arbitrary example
    r = np.repeat(1, no_angles)
    theta = np.repeat(sunPos[1], no_angles)  # MSASP will always allign to the polar angle of the sun
    rho = np.linspace(0, no_rotations * 2 * np.pi, no_angles + 1)[:-1]  # no_rotations rotations around its axes
    mSASPpos = np.array([r, theta, rho]).transpose()
    ### trans to cartesian coordinates
    sunPos = mla.spheric2cart(np.repeat(np.array([sunPos]), no_angles, axis=0))
    mSASPpos = mla.spheric2cart(mSASPpos)
    angles = mla.angleBetweenVectors(sunPos, mSASPpos)
    angles = pd.DataFrame(angles, index=rho, columns=['mSASP_sun_angle'])
    angles.index.name = 'mSASP_azimuth'
    return angles

def _simplefill(series):
    """Very simple function to fill missing values. Should only be used for 
    values which basically do not change like month and day. 
    Will most likely give strange results when the day does change.
    Returns: nothing everything happens inplace"""
    
    series.values[0] = series.dropna().values[0]
    series.values[-1] = series.dropna().values[-1]
    series.fillna(method='ffill', inplace = True)
    return


# Todo: wouldn't it be better if that would be a subclass of timeseries?
class miniSASP(object):
    def __init__(self, fname, verbose=True):
        self.verbose = verbose
        self.read_file(fname)
        # self.assureFloat()
        self.recover_negative_values()
        self.normalizeToIntegrationTime()
        self.set_dateTime()
        self.remove_data_withoutGPS()
        self.remove_unused_columns()
        self.channels = miniSASP_channels

    def _time2vertical_profile(self, stack, ts, key='Altitude'):
        """Converts the time series of mSASP revolutions to a height profile by merging
        with a TimeSeries that has height information"""
        picco_t = ts.copy()
        picco_t.data = picco_t.data[[key]]
        out = stack.merge(picco_t)
        out.data.index = out.data[key]
        out.data = out.data.drop(key, axis=1)
        #     out = out.data.transpose()
        return out

    # Todo: inherit docstring
    def get_sun_position(self):
        """read docstring of solar.get_sun_position_TS"""
        out = solar.get_sun_position_TS(self)
        return out

    def merge(self, ts):
        """ Merges current with other timeseries. The returned timeseries has the same time-axes as the current
        one (as opposed to the one merged into it). Missing or offset data points are linearly interpolated.

        Argument
        --------
        ts: timeseries or one of its subclasses.
            List of TimeSeries objects.

        Returns
        -------
        TimeSeries object or one of its subclasses

        """
        ts_this = self.copy()
        ts_data_list = [ts_this.data, ts.data]
        catsortinterp = pd.concat(ts_data_list).sort_index().interpolate()
        merged = catsortinterp.groupby(catsortinterp.index).mean().reindex(ts_data_list[0].index)
        ts_this.data = merged
        return ts_this

    def _smoothen(self, stack, window=20, which='vertical'):
        """Smoothen a vertical profile of mSASP revolutions.
        "which" is in case there is the need for smootingin something else, like the timeseries."""
        stack = stack.copy()
        out = stack.data
        heigh = round(out.index.values.max())
        low = round(out.index.values.min())
        no = heigh - low
        try:
            out = out.reindex(np.linspace(low, heigh, no + 1), method='nearest')
        except ValueError:
            warnings.warn('The height was not monotonic ant needed to be sorted!!!! This can be the origin of errors.')
            out = out.sort_index()
            out = out.reindex(np.linspace(low, heigh, no + 1), method='nearest')
        out = pd.rolling_mean(out, window, min_periods=1, center=True)
        out = out.reindex(np.arange(low + (window / 2), heigh, window), method='nearest')
        stack.data = out
        return stack

    def split_revolutions(self, peaks='l', time_delta=(5, 20), revolution_period=26.):
        """This function reorganizes the miniSASP data in a way that all sun transits are stacked on top of eachother
        and the time is translated to an angle"""

        ulr = self.copy()
        # star = 10000
        # till = 20000
        # ulr.data = ulr.data[star:till]


        if peaks == 's':
            peaks_s = ulr.find_peaks()
        elif peaks == 'l':
            peaks_s = ulr.find_peaks(which='long')

        time_delta_back = time_delta[0]
        time_delta_forward = time_delta[1]

        #     wls = ['460.3', '550.4',  '671.2', '860.7']
        photos = [ulr.data.PhotoA, ulr.data.PhotoB, ulr.data.PhotoC, ulr.data.PhotoD]
        out_dict = {}
        for u, i in enumerate(ulr.channels):

            centers = peaks_s.data[str(i)].dropna().index.values

            #     res = []
            df = pd.DataFrame()
            PAl = photos[u]
            for e, center in enumerate(centers):
                # center = peaks_s.data['460.3'].dropna().index.values[1]
                start = center - np.timedelta64(time_delta_back, 's')
                end = center + np.timedelta64(time_delta_forward, 's')
                PAlt = PAl.truncate(before=start, after=end, copy=True)
                PAlt.index = PAlt.index - center
                PAlt = PAlt[
                    PAlt != 0]  # For some reasons there are values equal to 0 which would screw up the averaging I intend to do
                #         res.append(PAlt)
                df[center] = PAlt.resample('50ms')
            df.index = (df.index.values - np.datetime64('1970-01-01T00:00:00.000000000Z')) / np.timedelta64(1, 's')
            df.index = df.index.values / revolution_period * 2 * np.pi
            out = timeseries.TimeSeries(df.transpose())
            out_dict[i] = out
        return out_dict

    def create_sky_brightness_altitude(self, picco, peaks='l', time_delta=(5, 20), revolution_period=26., key='Altitude',
                                       window=20, which='vertical'):
        """Creates a smoothened vertical profile of the sky brightness."""
        strevol = self.split_revolutions(peaks=peaks, time_delta=time_delta, revolution_period=revolution_period)
        smoothened_vertical_profs = SkyBrightDict()
        for i in strevol.keys():
            strevol_t = strevol[i]
            vprof_t = self._time2vertical_profile(strevol_t, picco, key=key)
            smoothened_vertical_profs[i] = self._smoothen(vprof_t, window=window, which=which)
        return smoothened_vertical_profs

    def copy(self):
        return deepcopy(self)

    def zoom_time(self, start=None, end=None, copy=True):
        """ Selects a strech of time from a housekeeping instance.

        Arguments
        ---------
        start (optional):   string - Timestamp of format '%Y-%m-%d %H:%M:%S.%f' or '%Y-%m-%d %H:%M:%S'
        end (optional):     string ... as start
        copy (optional):    bool - if False the instance will be changed. Else, a copy is returned

        Returns
        -------
        If copy is True:  housekeeping instance
        else:             nothing (instance is changed in place)


        Example
        -------
        >>> from atmPy.aerosols.instruments.piccolo import piccolo
        >>> launch = '2015-04-19 08:20:22'
        >>> landing = '2015-04-19 10:29:22'
        >>> hk = piccolo.read_file(filename) # create housekeeping instance
        >>> hk_zoom = zoom_time(hk, start = launch, end= landing)
        """

        if copy:
            housek = self.copy()
        else:
            housek = self

        if start:
            start = time_tools.string2timestamp(start)
        if end:
            end = time_tools.string2timestamp(end)

        housek.data = housek.data.truncate(before=start, after=end)
        if copy:
            return housek
        else:
            return

    def find_peaks(self, which='short', min_snr=10, moving_max_window=23):
        """ Finds the peaks in all four photo channels (short exposure). It also returns a moving maximum as guide to
         the eye for the more "real" values.

         Parameters
         ----------
         which: 'long' or 'short'.
            If e.g. PhotoA (long) or PhotoAsh(short) are used.
         min_snr: int, optional.
            Minimum signal to noise ratio.
         moving_max_window: in, optionl.
            Window width for the moving maximum.


         Returns
         -------
         TimeSeries instance (AtmPy)
        """
        moving_max_window = int(moving_max_window / 2.)
        # till = 10000
        # photos = [self.data.PhotoAsh[:till], self.data.PhotoBsh[:till], self.data.PhotoCsh[:till], self.data.PhotoDsh[:till]]
        if which == 'short':
            photos = [self.data.PhotoAsh, self.data.PhotoBsh, self.data.PhotoCsh, self.data.PhotoDsh]
        elif which == 'long':
            photos = [self.data.PhotoA, self.data.PhotoB, self.data.PhotoC, self.data.PhotoD]
        else:
            raise ValueError('which must be "long" or "short" and not %s' % which)
        # channels = [ 550.4, 460.3, 860.7, 671.2]
        df_list = []
        for e, data in enumerate(photos):
            pos_indx = signal.find_peaks_cwt(data, np.array([10]), min_snr=min_snr)
            out = pd.DataFrame()
            # return out
            # break
            out[str(self.channels[e])] = pd.Series(data.values[pos_indx], index=data.index.values[pos_indx])

            out['%s max' % self.channels[e]] = self._moving_max(out[str(self.channels[e])], window=moving_max_window)
            df_list.append(out)

        out = pd.concat(df_list).sort_index()
        out = out.groupby(out.index).mean()
        out = Sun_Intensities_TS(out)
        return out

    def _moving_max(self, ds, window=3):
        # x = ds.dropna().index.values
        # y = ds.dropna().values
        # out = []
        # out_x = []
        # i = 0
        # while True:
        #     out.append(y[i:i+window].max())
        #     out_x.append(x[int(i+window/2)])
        #     i = i+window
        #     if (i+window/2) >= len(x):
        #         break
        # out = np.array(out)
        # out_x = np.array(out_x)
        # return pd.Series(out, index = out_x)

        out = pd.DataFrame(ds, index=ds.index)
        out = pd.rolling_max(out, window)
        out = pd.rolling_mean(out, int(window / 5), center=True)
        return out

    def remove_data_withoutGPS(self, day='08', month='01'):
        """ Removes data from before the GPS is fully initiallized. At that time the Date should be the 8th of January.
        This is an arbitray value, which might change

        Arguments
        ---------
        day (optional): string of 2 digit integer
        month (optional): string of 2 digit integer
        """
        self.data = self.data[((self.data.Day != day) & (self.data.Month != month))]
    
    def read_file(self,fname):
        df = pd.read_csv(fname,
                         encoding="ISO-8859-1",
                         skiprows=16,
                         header=None,
                         error_bad_lines=False,
                         warn_bad_lines= False
                        )

        #### set column labels
        collabels = ['PhotoAsh', 
                     'PhotoBsh',
                     'PhotoCsh',
                     'PhotoDsh',
                     'PhotoA',
                     'PhotoB',
                     'PhotoC',
                     'PhotoD',
                     'Seconds',
                     'caseflag',
                     'var1',
                     'var2',
                     'var3']

        df.columns = collabels

        #### Drop all lines which lead to errors
        df = df.convert_objects(convert_numeric='force')
        df = df.dropna(subset=['Seconds'])
        # self.data = df
        # return
        df = df.astype(float)
        # self.data = df
        # return

        #### add extra columns
        df['time'] = np.nan
        df['azimuth'] = np.nan
        df['homep'] = np.nan
        df['MicroUsed'] = np.nan
        df['lat'] = np.nan
        df['lon'] = np.nan
        df['Te'] = np.nan
        df['GPSHr'] = np.nan
        df['MonthDay'] = np.nan
        df['Month'] = np.nan
        df['Day'] = np.nan
        df['GPSReadSeconds'] = np.nan
        df['HKSeconds'] = np.nan
        df['Yaw'] = np.nan
        df['Pitch'] = np.nan
        df['Roll'] = np.nan
        df['BaromPr'] = np.nan
        df['BaromTe'] = np.nan
        df['Modeflag'] = np.nan
        df['GateLgArr'] = np.nan
        df['GateShArr'] = np.nan
        df['PhotoOffArr'] = np.nan


        ##### Case 0
        case = np.where(df.caseflag.values == 0)
        df.azimuth.values[case] = df.var1.values[case]
        df.homep.values[case] = df.var2.values[case]
        df.MicroUsed.values[case] = df.var3.values[case]

        ##### Case 1
        case = np.where(df.caseflag.values == 1)
        df.lat.values[case] = df.var1.values[case]
        df.lon.values[case] = df.var2.values[case]
        df.Te.values[case] = df.var3.values[case]

        ##### Case 2
        case = np.where(df.caseflag.values == 2)
        df.GPSHr.values[case] = df.var1.values[case]
        df.MonthDay.values[case] = df.var2.values[case]
        df.GPSReadSeconds.values[case] = df.var3.values[case].astype(float) / 100.
        df.HKSeconds.values[case] = df.Seconds.values[case]

        ##### Case 3
        case = np.where(df.caseflag.values == 3)
        df.Yaw.values[case] = df.var1.values[case]
        df.Pitch.values[case] = df.var2.values[case]
        df.Roll.values[case] = df.var3.values[case]

        ##### Case 4
        case = np.where(df.caseflag.values == 4)
        df.BaromPr.values[case] = df.var1.values[case]
        df.BaromTe.values[case] = df.var2.values[case]
        df.Modeflag.values[case] = df.var3.values[case].astype(float) + 0.5

        ##### Case 5
        case = np.where(df.caseflag.values == 5)
        df.GateLgArr.values[case] = df.var1.values[case]
        df.GateShArr.values[case] = df.var2.values[case]
        df.PhotoOffArr.values[case] = df.var3.values[case]
        _simplefill(df.GateLgArr)
        _simplefill(df.GateShArr)
        self.data = df

    def assureFloat(self):
        """ OLD!!!! This function is currently not used
        Note
        ----
        Sometimes there is a line, which start with 'GPRMC'. This Line causes trouble. Therefore """

        where = np.where(self.data.PhotoAsh == 'GPRMC')
        self.data.PhotoAsh.values[where] = np.nan
        self.data.PhotoBsh.values[where] = np.nan
        self.data.PhotoCsh.values[where] = np.nan
        self.data.PhotoDsh.values[where] = np.nan

        self.data.PhotoA.values[where] = np.nan
        self.data.PhotoB.values[where] = np.nan
        self.data.PhotoC.values[where] = np.nan
        self.data.PhotoD.values[where] = np.nan

        self.data.PhotoAsh = self.data.PhotoAsh.values.astype(float)
        self.data.PhotoBsh = self.data.PhotoBsh.values.astype(float)
        self.data.PhotoCsh = self.data.PhotoCsh.values.astype(float)
        self.data.PhotoDsh = self.data.PhotoDsh.values.astype(float)

        self.data.PhotoA = self.data.PhotoA.values.astype(float)
        self.data.PhotoB = self.data.PhotoB.values.astype(float)
        self.data.PhotoC = self.data.PhotoC.values.astype(float)
        self.data.PhotoD = self.data.PhotoD.values.astype(float)
        
    def normalizeToIntegrationTime(self):
        ##### normalize to integration time
        self.data.PhotoAsh = self.data.PhotoAsh.values.astype(float)/self.data.GateShArr.values.astype(float)
        self.data.PhotoA = self.data.PhotoA.values.astype(float)/self.data.GateLgArr.values.astype(float)

        self.data.PhotoBsh = self.data.PhotoBsh.values.astype(float)/self.data.GateShArr.values.astype(float)
        self.data.PhotoB = self.data.PhotoB.values.astype(float)/self.data.GateLgArr.values.astype(float)

        self.data.PhotoCsh = self.data.PhotoCsh.values.astype(float)/self.data.GateShArr.values.astype(float)
        self.data.PhotoC = self.data.PhotoC.values.astype(float)/self.data.GateLgArr.values.astype(float)

        self.data.PhotoDsh = self.data.PhotoDsh.values.astype(float)/self.data.GateShArr.values.astype(float)
        self.data.PhotoD = self.data.PhotoD.values.astype(float)/self.data.GateLgArr.values.astype(float)



    def set_dateTime(self, millisscale = 10):
        self.data.Seconds *= (millisscale/1000.)
        self.data.index = self.data.Seconds

        self.data.Month = np.floor(self.data.MonthDay/100.)
        _simplefill(self.data.Month)

        self.data.Day = self.data.MonthDay - self.data.Month*100
        _simplefill(self.data.Day)
        self.data.Month = self.data.Month.astype(int).apply(lambda x: '{0:0>2}'.format(x))
        self.data.Day = self.data.Day.astype(int).apply(lambda x: '{0:0>2}'.format(x))
        self.GPSHr_P_3 = self.data.GPSHr.copy()
        self.Month_P_1 = self.data.Month.copy()
        self.Day_P_1 = self.data.Day.copy()
        ## time from GPS
        # get rid of stepfunktion
#         GPSunique = df.GPSReadSeconds.dropna().unique()
#         for e,i in enumerate(GPSunique):
#             where = np.where(df.GPSReadSeconds == i)[0][1:]
#             self.data.GPSReadSeconds.values[where] = np.nan

        GPSunique = self.data.GPSHr.dropna().unique()
        for e,i in enumerate(GPSunique):
            where = np.where(self.data.GPSHr == i)[0][1:]
            self.data.GPSHr.values[where] = np.nan

        self.GPSHr_P_1 = self.data.GPSHr.copy()
        # extrapolate and interpolate the time    
        _extrapolate(self.data.index, self.data.GPSHr)
        self.data.GPSHr.interpolate(method='index', inplace= True)
        self.data.GPSHr.dropna(inplace=True)
        self.GPSHr_P_2 = self.data.GPSHr.copy()
        self.data.GPSHr = self.data.GPSHr.apply(lambda x: '%02i:%02i:%09.6f'%(x,60 * (x % 1), 60* ((60 * (x % 1)) %1)))

        ###### DateTime!!
        dateTime = year + '-' +  self.data.Month + '-' + self.data.Day +' ' + self.data.GPSHr
        self.data.index = pd.Series(pd.to_datetime(dateTime, format="%Y-%m-%d %H:%M:%S.%f"), name='Time_UTC')

        self.data = self.data[pd.notnull(self.data.index)]  # gets rid of NaT
    
    def recover_negative_values(self):
        """ this is most likely not working!"""
        _recover_Negatives(self.data.PhotoA, verbose=self.verbose)
        _recover_Negatives(self.data.PhotoB, verbose=self.verbose)
        _recover_Negatives(self.data.PhotoC, verbose=self.verbose)
        _recover_Negatives(self.data.PhotoD, verbose=self.verbose)
        _recover_Negatives(self.data.PhotoAsh, verbose=self.verbose)
        _recover_Negatives(self.data.PhotoBsh, verbose=self.verbose)
        _recover_Negatives(self.data.PhotoCsh, verbose=self.verbose)
        _recover_Negatives(self.data.PhotoDsh, verbose=self.verbose)
    
    def remove_unused_columns(self):
        self.data.drop('var1', axis=1, inplace= True)
        self.data.drop('var2', axis=1, inplace= True)
        self.data.drop('var3', axis=1, inplace= True)
        self.data.drop('time', axis=1, inplace= True)
        self.data.drop('GPSHr', axis=1, inplace= True)
        self.data.drop('MonthDay', axis=1, inplace= True)
        self.data.drop('Month', axis=1, inplace= True)
        self.data.drop('Day', axis=1, inplace= True)
        self.data.drop('GPSReadSeconds', axis=1, inplace= True)
        self.data.drop('GateLgArr', axis=1, inplace= True)
        self.data.drop('GateShArr', axis=1, inplace= True)


def load_sunintensities_TS(fname):
    data = pd.read_csv(fname, index_col=0)
    data.index = pd.to_datetime(data.index)
    return Sun_Intensities_TS(data)


class Sun_Intensities_TS(timeseries.TimeSeries):
    def plot(self, offset=[0, 0, 0, 0], airmassfct=True, move_max=True, legend=True, all_on_one_axis = False,
             additional_axes=False,
             errors = False,
             rayleigh=True):
        """plots ... sorry, but this is a messi function. Things should have been done different, e.g too much data
         processing whith the data not put out ... need fixn
        Arguments
        ---------
        offset: list
        airmassfct: bool.
            If the airmass factor is included or not.
            True: naturally the air-mass factor is included in the data, so this does nothing.
            False: data is corrected to correct for the slant angle
        rayleigh: bool or the aod part of the output of miniSASP.simulate_from_size_dist_LS.
            make sure there is no airmassfkt included in this!!
        all_on_one_axis: bool or axes instance
            if True all is plotted in one axes. If axes instances this axis is used.
        """

        m_size = 5
        m_ewidht = 1.5
        l_width = 2
        gridspec_kw = {'wspace': 0.05}
        no_axes = 4
        if all_on_one_axis:
            no_axes = 1
        if additional_axes:
            no_axes = no_axes + additional_axes

        if type(all_on_one_axis).__name__ == 'AxesSubplot':
            a = all_on_one_axis
            f = a.get_figure()
        else:
            f, a = plt.subplots(1, no_axes, gridspec_kw=gridspec_kw)
        columns = ['460.3', '460.3 max', '550.4', '550.4 max', '671.2', '671.2 max', '860.7', '860.7 max']
        # peaks_max = [460.3, '460.3 max', 550.4, '550.4 max', 860.7, '860.7 max', 671.2,
        #        '671.2 max']
        if not all_on_one_axis:
            f.set_figwidth(15)
        #################
        for i in range(int(len(columns) / 2)):
            col = plt_tools.wavelength_to_rgb(columns[i * 2]) * 0.8
            intens = self.data[columns[i * 2]].dropna()  # .plot(ax = a, style = 'o', label = '%s nm'%colums[i*2])
            x = intens.index.get_level_values(1)
            if type(rayleigh) == bool:
                if rayleigh:
                    rayleigh_corr = 0
            else:
                # print('mach ick')
                aodt = rayleigh[float(columns[i * 2])].loc[:, ['rayleigh']]
                intenst = intens.copy()
                intenst.index = intenst.index.droplevel(['Time', 'Sunelevation'])
                aodt_sit = pd.concat([aodt, intenst]).sort_index().interpolate()
                aodt_sit = aodt_sit.groupby(aodt_sit.index).mean().reindex(intenst.index)
                rayleigh_corr = aodt_sit.rayleigh.values / np.sin(intens.index.get_level_values(2))
                # return aodt

            if not airmassfct:
                amf_corr = np.sin(intens.index.get_level_values(2))
            else:
                amf_corr = 1
            if not all_on_one_axis:
                atmp = a[i]
            else:
                atmp = a


            y = (offset[i] - np.log(intens) - rayleigh_corr) * amf_corr
            g, = atmp.plot(y, x)
            g.set_label('%s nm' % columns[i * 2])
            g.set_linestyle('')
            g.set_marker('o')
            #         g = a.get_lines()[-1]
            g.set_markersize(m_size)
            g.set_markeredgewidth(m_ewidht)
            g.set_markerfacecolor('None')
            g.set_markeredgecolor(col)

            if move_max:
                #             sun_intensities.data.iloc[:,i*2+1].dropna().plot(ax = a)
                intens = self.data[
                    columns[i * 2 + 1]].dropna()  # .plot(ax = a, style = 'o', label = '%s nm'%colums[i*2])
                x = intens.index.values

                g, = a[i].plot(offset[i] - np.log(intens), x)
                #             g = a.get_lines()[-1]
                g.set_color(col)
                #             g.set_solid_joinstyle('round')
                g.set_linewidth(l_width)
                g.set_label(None)

            if i != 0 and not all_on_one_axis:
                atmp.set_yticklabels([])

            if i == 4:
                break
        if all_on_one_axis:
            a.legend()
        else:
            if legend:
                for aa in a:
                    aa.legend()
        if not airmassfct:
            txt = 'OD'
        else:
            txt = 'OD * (air-mass factor)'
        if all_on_one_axis:
            atmp = a
        else:
            atmp = a[0]
        atmp.set_xlabel(txt)
        if not all_on_one_axis:
            atmp.xaxis.set_label_coords(2.05, -0.07)
        atmp.set_ylabel('Altitude (m)')
        return a

    def add_sun_elevetion(self, picco):
        """
        doc is not correct!!!

        This function uses telemetry data from the airplain (any timeseries including Lat and Lon) to calculate
        the sun's elevation. Based on the sun's elevation an airmass factor is calculated which the data is corrected for.

        Arguments
        ---------
        sun_intensities: Sun_Intensities_TS instance
        picco: any timeseries instance containing Lat and Lon
        """

        picco_t = timeseries.TimeSeries(picco.data.loc[:, ['Lat', 'Lon', 'Altitude']])  # only Altitude, Lat and Lon
        sun_int_su = self.merge(picco_t)
        out = sun_int_su.get_sun_position()
        #     sun_int_su = sun_int_su.zoom_time(spiral_up_start, spiral_up_end)
        arrays = np.array([sun_int_su.data.index, sun_int_su.data.Altitude, sun_int_su.data.Solar_position_elevation])
        tuples = list(zip(*arrays))
        index = pd.MultiIndex.from_tuples(tuples, names=['Time', 'Altitude', 'Sunelevation'])
        sun_int_su.data.index = index
        sun_int_su.data = sun_int_su.data.drop(
            ['Altitude', 'Solar_position_elevation', 'Solar_position_azimuth', 'Lon', 'Lat'], axis=1)
        return sun_int_su


def load_skybrighness(fname):
    keys = [460.3, 860.7, 550.4, 671.2]
    outt = SkyBrightDict()
    for k in keys:
        fn = fname + '_' + str(k) + '.csv'
        df = pd.read_csv(fn, index_col=0)
        df.columns = df.columns.astype(float)
        outt[float(k)] = timeseries.TimeSeries(df)
    return outt


class SkyBrightDict(dict):
    def save(self, fname):
        for k in self.keys():
            fn = fname + '_' + str(k) + '.csv'
            df = self[k]
            df.save(fn)

    def plot(self):
        """plot the output of mSASP.create_sky_brightness_altitude.
        Returns
        -------
        - axes
        - legend"""
        f, ax = plt.subplots(2, 2)
        f.set_size_inches((15, 15))
        ax = ax.flatten()
        sl = list(self.keys())
        sl.sort()
        for e, i in enumerate(sl):
            a = self[i].data.transpose().plot(ax=ax[e])  # scalarMap.to_rgba(e))
            lines = a.get_lines()
            scalarMap = cm.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=len(lines)), cmap=plt.cm.brg)
            for u, l in enumerate(lines):
                l.set_color(scalarMap.to_rgba(u))
            a.set_title(i)
            a.set_ylim((0.003, 0.025))
            a.set_xlabel('azimuth (rad)')
            l = a.legend(prop={'size': 8})
        return ax,l
