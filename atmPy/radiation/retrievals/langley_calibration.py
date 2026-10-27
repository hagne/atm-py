#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 21:01:37 2022

@author: hagen
"""

# -*- coding: utf-8 -*-
import pandas as pd
# from scipy import stats
import numpy as np
import pathlib
from typing import Union
# from scipy.fftpack import dst
import xarray as xr
import datetime
from atmPy.opt_imports import plt
# import pptx
import io
# import statsmodels.api as sm
# import sp02.products.raw_nc
import copy
import pathlib as pl
# from pygam import f as pgf
import scipy.interpolate as scint
import scipy.stats
# import matplotlib.dates as pltdates

import atmPy.radiation.retrievals.langley_db as atmlldb
from atmPy.opt_imports import pptx
from atmPy.opt_imports import statsmodels
import atmPy.general.timeseries as atmts
import atmPy.opt_imports as opt_imports
pltdates = opt_imports.OptionalImport('matplotlib.dates')


def weighted_mean_dt(dft, center_datetime, order_stderr = 1, oder_dt = 1, minperiods = 1, stderr_threshold = 0.02, dt_window = 60, approach = 1):
    """I tried several approaches and to me the more suffisticated approaches did not yield the result that I would have hoped.
    
    approach 1 (default)
    --------------------
    simply 1/V0_stderr**order_stderr with order_stderr = 2 and stderr_threshold = 0.02
    
    approach 2
    -----------
    This considers in addition to the above the temporal distance (dt) of each V0 to the current date. 
    Here dt is decreasing the weight of V0 that are further away. Note, A more realistic approach is the next 
    one though which considers an increase in V0_stderr with dt. 
    
    approach 3
    -----------
    Temporal distance (dt) of each V0 to the current date increases V0_stderr. While this approach is probably the most reasist, 
    it just does not seem to make a big difference, which is why I will probably with the first approach
    """
    dft = dft.copy()
    dft[dft.intercept_stderr > stderr_threshold] = np.nan
    out = {}
    if dft.dropna().shape[0] < minperiods:
        wm = np.nan
        wmst = np.nan
        
    elif approach == 1: # approach 1
        weights = (1/dft.intercept_stderr)**order_stderr
        
        wm = (np.e**(dft.intercept) * weights).sum() / weights.sum()
        wmst = ((np.e**(dft.intercept) - wm)**2 * weights).sum() / weights.sum()
        # print(wmst)
        wmst = np.sqrt(wmst)
    elif approach == 2: # approach 2
        # Here I tried to 
        weights_err = (1 - (dft.intercept_stderr/stderr_threshold))**order_stderr
        td = abs(dft.index - center_datetime) / pd.to_timedelta(1, 'd')
        weights_dt = pd.Series((1 - (td / (dt_window/2)))**oder_dt, index=dft.index)
        weights = weights_err * weights_dt
        
        wm = (np.e**(dft.intercept) * weights).sum() / weights.sum()
        wmst = ((np.e**(dft.intercept) - wm)**2 * weights).sum() / weights.sum()
        wmst = np.sqrt(wmst)
    elif approach == 3: # approach 3
        td = abs(dft.index - center_datetime) / pd.to_timedelta(1, 'd')
        weights_dt = pd.Series((td / (dt_window/2))**(1/(oder_dt+0.0001)), index=dft.index)
        intererr = dft.intercept_stderr + (weights_dt * stderr_threshold * 0.25) # this will cause the uncertatinty to be increased by up to the stderr_threshold at the edge of the window
        weights_err = 1/(intererr**order_stderr)
        weights_err[weights_err < 0] = np.nan
        weights = weights_err# * weights_dt
        
        wm = (np.e**(dft.intercept) * weights).sum() / weights.sum()
        wmst = ((np.e**(dft.intercept) - wm)**2 * weights).sum() / weights.sum()
        wmst = np.sqrt(wmst)
    else:
        assert(False), f'{approach} not a valid value for kwarg "approach".'

    # out['weights_err'] = weights_err
    # out['weights_dt'] = weights_dt
    # out['weights'] = weights
    # out['dft'] = dft
        
    out['weighted_mean'] = wm
    out['weighted_std'] = wmst
    return out

def read_langley_params(p2f = '/home/grad/surfrad/aod/tbl_mfrhead', verbose = False):
    """
    This reads the fit parameters that John is generating. These parameters can be used to 
    reconstruct the V0s.

    Parameters
    ----------
    p2f : TYPE, optional
        DESCRIPTION. The default is '/home/grad/surfrad/aod/tbl_mfrhead'.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    p2f = pl.Path(p2f)

    with open(p2f, 'r') as rein:
        lines = rein.readlines()



    blocksize = 15
    # blocks = []
    dsses = []
    for i in range(100):
        if verbose:
            print(f'block: {i}')
        block = lines[i * blocksize: (i+1) * blocksize]
        if len(block) == 0:
            if verbose:
                print('done')
            break
        # blocks.append(block)
        # break
        dt = block[0].split('!')[0].strip()
        dt = pd.to_datetime(dt, format = '%y%j%H%M')

        mfrserial = int(block[1].strip())

        V0 = pd.DataFrame([[float(i) for i in l.split('!')[0].strip().split(' ')] for l in block[3::2]], columns=['const', 'lin', 'sin', 'cos'], index = [l.split('!')[1].strip().split(' ')[0] for l in block[3::2]])

        V0_err = pd.DataFrame([[float(i) for i in l.split('!')[0].strip().split(' ')] for l in block[4::2]], columns=['intercept', 'slope'], index = [l.split('!')[1].strip().split(' ')[0] for l in block[4::2]])
        
        if '615' in V0.index:
            if verbose:
                print('Has 615 nm channel -> continue')
            continue
        
        #### rename stats channels to nominal channels
        wl_nominal = np.array([415,500,1625, 615, 670, 870, 940])
        V0.rename({col: wl_nominal[abs(wl_nominal - int(col)).argmin()] for col in V0.index}, axis = 0, inplace = True)  
        V0_err.rename({col: wl_nominal[abs(wl_nominal - int(col)).argmin()] for col in V0_err.index}, axis = 0, inplace = True)  
        # V0.rename({'162': '1625',
        #            '670': '673'}, inplace=True)
        # V0_err.rename({'162': '1625',
        #                '670': '673'}, inplace=True)
        # return V0
        assert(np.all([wl in wl_nominal for wl in V0.index])), f'something is strange with the channel {V0.index}'
        assert(np.all([wl in wl_nominal for wl in V0_err.index ])), f'something is strange with the channel {V0_err.index}'
        # assert(np.all([wl in V0.index for wl in ['415', '500', '1625', '670', '870', '940']])), f'something is strange with the channel {V0.index}'
        # assert(np.all([wl in V0_err.index for wl in ['415', '500', '1625', '670', '870', '940']])), f'something is strange with the channel {V0_err.index}'
        V0.index = [int(col) for col in V0.index]
        V0.index.name = 'channel'
        V0.columns.name= 'V0_params'
        
        V0_err.index = [int(col) for col in V0_err.index]
        V0_err.index.name = 'channel'
        V0_err.columns.name= 'V0_err_params'
        
        ds = xr.Dataset()
        ds['V0'] = V0
        ds['V0_err'] = V0_err
        # ds = ds.expand_dims({'serial': [mfrserial]})
        ds['serial'] = mfrserial
        ds = ds.expand_dims({'datetime': [dt]})
        if verbose:
            print(dt)
        dsses.append(ds)
    return xr.concat(dsses, dim = 'datetime')

def open_langleys(p2fld: Union[str, pl.Path, list], pattern = '*.nc'):
    """This function is designed to open all langleys within a folder. Here the langleys 
    are those generated and saved in the Langley class (with Langley.save2netcdf()
    Parameters
    ----------
    p2fld : Union[str, pl.Path, list]
        This can be either a path to a folder containing the langley files or a list of paths 
        to the langley files or a list of already opened langley datasets.
    """

    if  isinstance(p2fld, (str, pl.Path)):
        p2fld = pl.Path(p2fld)
        p2flist = list(p2fld.glob(pattern))
    elif isinstance(p2fld, list):
        if all([isinstance(p, xr.Dataset) for p in p2fld]):
            p2flist = p2fld
        else:
            p2flist = [pl.Path(p) for p in p2fld]

    p2flist.sort()
    dslist = []
    for p2f in p2flist:
        if isinstance(p2f, xr.Dataset):
            ds = p2f
        else:
            ds = xr.open_dataset(p2f)
        if 'date' in ds.attrs:
            dt = pd.to_datetime(ds.attrs['date'])
        else:
            dt = pd.to_datetime(p2f.name.split('_')[-1].split('.')[0])
        if 'when' in ds.attrs:
            ampm = ds.attrs['when']
        else:
            ampm = p2f.name.split('_')[1]
        if ampm == 'pm':
            dt += pd.to_timedelta(6, 'h')
        ds = ds.expand_dims(datetime = [dt])
        # ds = ds.expand_dims(ampm = [ampm])
        ds['ampm'] =  ({'datetime':[dt]}, [ampm])
        dslist.append(ds)
    ds = xr.concat(dslist, 'datetime', join='outer')
    return Langley_Timeseries(ds)

def read_langleys_from_database(path2db, start = None, end = None, days = None):
    db = atmlldb.LangleyResultDatabase(path2db)
    start = pd.to_datetime(start) if start is not None else None
    end = pd.to_datetime(end) if end is not None else None
    days = pd.to_timedelta(days, 'd') if days is not None else None

    if start is not None and end is not None and days is not None:
        raise ValueError('one of start, end, and days must be None')
    if start is None and end is None and days is None:
        raise ValueError('at least one of start, end, and days must be given')
    if start is not None and end is not None:
        pass
    elif end is not None and days is not None:
        start = end - days
    elif start is not None and days is not None:
        end = start + days
    elif days is not None:
        end = pd.Timestamp.now()
        start = end - days
    elif start is not None:
        end = pd.Timestamp.now()
    else:
        raise ValueError('end alone cannot determine the start date')

    ds = db.dump_langley_results(local_day_range=(start, end),channel = None)

    return Langley_Timeseries(ds)

def open_langley_dailys(start = None, #'20160101',  
                        end = None, #'20220101',
                        p2fld = '/home/grad/htelg/data/grad/surfrad/mfrsr/langleys_concat/tbl/',
                        drop_variables = ['langley_residual_correlation_prop', 'valid_points', 'residual_stats', 'cleaning_iterations', 'status'], # they are not needed under normal conditions
                       ):
    """This is something used within SURFRAD, this should not be here or made more general"""
    p2fld = pl.Path(p2fld)
    df = pd.DataFrame(p2fld.glob('*.nc'), columns=['p2f'])
    df.index = df.apply(lambda row: pd.to_datetime(row.p2f.name.split('_')[4] + '01'), axis =1)
    df.sort_index(inplace=True)
    df = df.truncate(start, end)
    assert(df.shape[0] != 0), f'no files found in {p2fld}.'
    df['serial_no'] = df.apply(lambda row: row.p2f.name.split('_')[3], axis = 1)
    assert(df.serial_no.unique().shape[0] != 0), f'No serial_no found in {df.p2f}.'
    assert(df.serial_no.unique().shape[0] == 1), f'Files indicate more then one serial number (found {df.serial_no.unique()}), which should not be the case unless you updated the langley processing ... did you? ... programmin grequired.'   
    # return df
    # ds_langs = xr.open_mfdataset(df.p2f, drop_variables = drop_variables)
    with xr.open_mfdataset(df.p2f, drop_variables = drop_variables) as ds_langs:
        ds_langs.load()
    return Langley_Timeseries(ds_langs)

class CalibrationPrediction(object):
    def __init__(self, dataset):
        self.dataset = dataset
    
    def plot(self, wl=500, 
                 fitres = False, 
                 ax = None,
                 show_pm5 = True,
                 show_uncertainty = True, 
                 **kwargs):
        
        if isinstance(ax, type(None)):
            a = plt.subplot()
        else:
            a = ax
            
        dst = self.dataset.sel(channel = wl)
        if show_pm5:
            perc = 0.05
            a.fill_between(dst.datetime,  dst.V0 - ( dst.V0 * perc),  dst.V0 + ( dst.V0 * perc), 
                           color = '0.9', zorder = 0)
        
        if 'label' in kwargs:
            label = kwargs.pop('label')
        else:
            label = 'wl'
        g, = a.plot(dst.datetime, dst.V0, label = label, zorder = 10, **kwargs)
        if show_uncertainty:
            col = g.get_color()
            a.fill_between(dst.datetime, dst.confidence_interval.sel(interval_boundary = 'low'), dst.confidence_interval.sel(interval_boundary = 'high'), color = col, alpha = 0.3)
            
        a.set_ylim((dst.V0 - ( dst.V0 * perc)).min(), (dst.V0 + ( dst.V0 * perc)).max())
        return a
        

class Langley_Timeseries(object):
    def __init__(self, dataset):
        if isinstance(dataset, xr.Dataset):
            self.dataset = dataset
        elif isinstance(dataset, str):
            self.dataset = xr.open_dataset(dataset)
        self.predict_until = self.dataset.datetime.values[-1] + pd.to_timedelta(60, 'days')
        # self.daterange2predict = pd.date_range(start = self.dataset.datetime.values[0], end = end, freq='D')
        self._v0rol = None
        self._v0gam = None
        self._v0 = None
        self._ranked = None    
        self._threshold_mask = None
        self._prediction = None

    # def rank_by(self, wl = 500, fit_results = 'intercept_stderr'):
    #     """Ranks the current dataset by the given values"""
    #     ds = self.dataset
    #     dssel = ds.langley_fitres.sel(fit_results = fit_results, channel = wl, drop = True)
    #     ds['ranked'] = dssel.rank('datetime')
    #     ds['ranked'].attrs['description'] = f'Ranked by langley_fitres {fit_results} at {wl}nm'
    #     return
    
    @property
    def daterange2predict(self):
        return pd.date_range(start = self.dataset.datetime.values[0], end = self.predict_until, freq='D')
    
    def sort_by(self, wl = 500, fit_results = 'intercept_stderr'):
        """Sorts the current dataset by the given values"""
        ds = self.dataset
        ds = ds.sortby(ds.langley_fitres.sel(channel = wl, fit_results = fit_results))
        ds['ranked'] = ('datetime', (range(ds.datetime.shape[0])))
        self.dataset = ds
        return

    def plot_sorted(self, wl = 500, wlsort: int | None = None, sort_by = 'intercept_stderr', aa = None):
        """Plots the langley fit results ranked by their standard error. Note, this will sort the dataset!!
        Parameters
        ----------
        wl: int
            Channel to plot
        wlsort: int 
            Channel to sort by (if None, use wl)
        """

        # if 'ranked' not in self.dataset:
        #     self.rank_by(wl = wl)
        # ds = self.dataset.swap_dims({'datetime':'ranked'})
        # ds = ds.sortby(ds.ranked)
        if isinstance(wlsort, type(None)):
            wlsort = wl

        self.sort_by(wl = wlsort, fit_results = sort_by)
        ds = self.dataset

        # ds = ds.sortby(ds.langley_fitres.sel(channel = wlsort, fit_results = sort_by))
        # ds['ranked'] = ('datetime', (range(ds.datetime.shape[0])))
        # self.dataset = ds
        if isinstance(aa, type(None)):
            f,aa = plt.subplots(2, sharex=True, gridspec_kw={'hspace':0})
        else:
            f = aa[0].figure

        a = aa[0]
        # a.plot(ds.ranked, ds.langley_fitres.sel(fit_results = 'intercept', channel = wl, drop = True))
        a.errorbar(ds.ranked, ds.langley_fitres.sel(fit_results = 'intercept', channel = wl), ds.langley_fitres.sel(fit_results = 'intercept_stderr', channel = wl))
        if hasattr(a, 'at'):
            at = a.at
        else:
            at = a.twinx()
            a.at = at
        # next(a._get_lines.prop_cycler)
        at._get_lines.get_next_color()
        at.plot(ds.ranked, ds.langley_fitres.sel(fit_results = 'intercept_stderr', channel = wl), marker = '.',
                # ls = ''
                # color = at._get_lines.get_next_color()
            )
        a.set_ylabel('ln(V0)')
        at.set_ylabel('std err ln(V0)')
        ##################
        a = aa[1]
        # a.plot(ds.ranked, ds.langley_fitres.sel(fit_results = 'intercept', channel = wl, drop = True))
        a.errorbar(ds.ranked, ds.langley_fitres.sel(fit_results = 'slope', channel = wl), ds.langley_fitres.sel(fit_results = 'slope_stderr', channel = wl))
        if hasattr(a, 'at'):
            at = a.at
        else:
            at = a.twinx()
            a.at = at
        # next(a._get_lines.prop_cycler)
        at._get_lines.get_next_color()
        at.plot(ds.ranked, ds.langley_fitres.sel(fit_results = 'slope_stderr', channel = wl), marker = '.',
                # ls = ''
                # color = at._get_lines.get_next_color()
            )
        a.set_ylabel('slope')
        at.set_ylabel('std err slope')
        return f,aa

    @property
    def threshold_mask(self):
        if self._threshold_mask is None:
            fitresults = self.dataset.fitresults
            coverage_metrics = self.dataset.coverage_metrics

            slope = fitresults.sel(fitresult='slope', drop=True)
            intercept = fitresults.sel(fitresult='intercept', drop=True)
            slope_stderr = fitresults.sel(fitresult='slope_stderr', drop=True)
            intercept_stderr = fitresults.sel(
                fitresult='intercept_stderr', drop=True
            )
            resid_curvature = fitresults.sel(
                fitresult='resid_curvature', drop=True
            )

            effective_airmass_offset = coverage_metrics.sel(
                coverage_metric='effective_airmass_offset', drop=True
            )
            maximum_cdf_gap_ratio = coverage_metrics.sel(
                coverage_metric='maximum_cdf_gap_ratio', drop=True
            )
            retained_points = coverage_metrics.sel(
                coverage_metric='retained_points', drop=True
            )
            retained_fraction = coverage_metrics.sel(
                coverage_metric='retained_fraction', drop=True
            )
            cdf_spacing_cv = coverage_metrics.sel(
                coverage_metric='cdf_spacing_cv', drop=True
            )

            slope_lower_bound = slope.quantile(0.5,'datetime') * 1.3
            intercept_stderr_upper_bound = 7 * intercept.median('datetime') * 1e-4
            self.tp_intercept_stderr_upper_bound = intercept_stderr_upper_bound
            self._threshold_mask = xr.Dataset({
                'slope': (slope_lower_bound < slope) & (slope < 0),
                'intercept': intercept > 0,
                'intercept_stderr': (
                    intercept_stderr < intercept_stderr_upper_bound
                ),
                'slope_stderr': slope_stderr < 0.001,
                'resid_curvature': resid_curvature < 0.005,
                'effective_airmass_offset': abs(effective_airmass_offset) < 0.2,
                'maximum_cdf_gap_ratio': maximum_cdf_gap_ratio < 50,
                'retained_points': retained_points > 100,
                'retained_fraction': retained_fraction > 0.8,
                'cdf_spacing_cv': cdf_spacing_cv < 1,
            })
        return self._threshold_mask

    def plot(self, wl = 500, th = 0.2, order_stderr = 2, ampm = 'am', ax = None, **kwargs):
        """Plots all langleys as a scatterplot where the color is a measure of the standard error.
        
        Parameters
        ----------
        wl: int
        ...
        """

        if isinstance(ax, type(None)):
            a = plt.subplot()
        else:
            a = ax
        intercept = np.exp(self.dataset.sel(channel = wl, fit_results = 'intercept', ampm = ampm).langley_fitres)
        intercept_stdr = self.dataset.sel(channel = wl, fit_results = 'intercept_stderr', ampm = ampm).langley_fitres
        weights = 1/intercept_stdr**order_stderr
        intercept = intercept.where(intercept_stdr < th)
        # return intercept.to_pandas()
        df = pd.DataFrame({'intercept':intercept.to_pandas(),
                      'intercept_stderr': intercept_stdr.to_pandas(),
                      'weights': weights.to_pandas()})#, columns=['basdf', 'asd'])
    
        df.sort_values('intercept_stderr', ascending=False, inplace=True)
        a.scatter(df.index, df.intercept, s = 8, c = df.weights, cmap =plt.cm.Greens, **kwargs)
        return a 

    
   
    def _get_v0_gam_old(self,
                    th = 0.02, # this is about half the annual variability (winter vs summer) of ~0.035
                    order_stderr = 2, # determines the weights
                    lam_overtime = 2.5e4,
                    ns_overtime = 2, # 4 per year -> 
                    lam_season = 1e4,
                    ns_season = 6,
                    ):
        # importing here so not everyone has to install pygam
        from pygam import LinearGAM, s #, l, GammaGAM
        
        lam_overtime *= self.dataset.datetime.shape[0]
        ds_v0_list = []
        
        for wl in self.dataset.channel:
            wl = int(wl)
            
            ds_v0 = xr.Dataset()
        
            dfwl =self.dataset.langley_fitres.sel(channel = wl).to_pandas().copy() #not sure if the copy is needed, but have to make sure that I don't overwrite the dataset
            no_of_years = (dfwl.iloc[-1].name - dfwl.iloc[0].name) / pd.to_timedelta(365, 'days')
            
            dfwl[dfwl.intercept_stderr > th] = np.nan
            dfwl.dropna(inplace=True)
        
            dfwl['doy']= dfwl.apply(lambda row: row.name.day_of_year, axis = 1)
            dfwl['todmonth'] = dfwl.apply(lambda row: (row.name - dfwl.iloc[0].name) / pd.to_timedelta(1, 'days'), axis = 1) #total number of month
        
            # new_daterange
            new_X = pd.DataFrame(index = self.daterange2predict, columns=['doy', 'todmonth'])
            new_X['doy']= new_X.apply(lambda row: row.name.day_of_year, axis = 1)
            new_X['todmonth'] = new_X.apply(lambda row: (row.name - new_X.iloc[0].name) / pd.to_timedelta(1, 'days'), axis = 1) #total number of month
                
            # gam parameters
            lam = lam_overtime
            n_splines = round(ns_overtime*no_of_years)
            todmonth = s(0, lam =  lam,
                        n_splines=n_splines,)
            
            lam = lam_season
            n_splines = ns_season
            doy = s(1, basis = 'cp', lam = lam, 
                    n_splines=n_splines,)
        
            y = np.exp(dfwl.intercept).values
            X = dfwl.loc[:, ['todmonth','doy',]].values
        
            weights = 1/dfwl.intercept_stderr**order_stderr
        
            # gam fit
            gammod = LinearGAM(todmonth + doy)
            gam = gammod.fit(X, y, 
                             weights=weights, 
                            )
        
            # prediction onto gap-less timeseries and extrapolate from last v0
            X = new_X.loc[:, ['todmonth','doy',]].values
            pred = gam.predict(X)
            new_X['pred_w2'] = pred
            conf = gam.confidence_intervals(X,width=0.99)
            new_X['conv_w2_low']= conf[:,0]
            new_X['conv_w2_high']= conf[:,1]
        
            new_X.index.name = 'datetime'
        
            ds_v0['V0'] = new_X.pred_w2
        
            ci = new_X.loc[:,['conv_w2_low', 'conv_w2_high']]
            ci.columns = ['low', 'high']
            ci.columns.name = 'interval_boundary'
            ds_v0['confidence_interval'] = ci
        
            ds_v0 = ds_v0.expand_dims(channel = [wl,])
        
            ds_v0_list.append(ds_v0)
        
        ds_v0_all = xr.concat(ds_v0_list, 'channel')
        self._v0gam = CalibrationPrediction(ds_v0_all)
        return self._v0gam
        
    def _get_v0_rolling(self,
                        th = 0.02, # this is about half the annual variability (winter vs summer) of ~0.035
                        dt_window = 60,
                        order_stderr = 2,
                        oder_dt = None,
                        approach = 1,
                        smoothned = False, # smoothening does not work well, maybe for further development?
                        ):

        ds_list = []
        wm_list = []
        # new_daterange = pd.date_range(start = self.dataset.datetime.values[0], end = self.dataset.datetime.values[-1], freq='D')
    
        for wl in self.dataset.channel:
            wl = int(wl)
            # break
        
            # wl = 500
        
            dfwl =self.dataset.langley_fitres.sel(channel = wl).to_pandas()
            rollw = dfwl.rolling(f'{dt_window}d', 
                                center = True,
                                   # min_periods =1,
                                  )
        
            mean_wl = pd.DataFrame([weighted_mean_dt(ro, idx, approach = approach, order_stderr = order_stderr, oder_dt = order_stderr, stderr_threshold = th, dt_window = dt_window) for idx,ro in zip(rollw.max().index, rollw)], index = rollw.mean().index, 
                                  columns = ['weighted_mean', 'weighted_std'],
                                 )
        
            wm = mean_wl.weighted_mean
            std = mean_wl.weighted_std
            confint = pd.DataFrame({'low':wm-std, 'high':wm+std})
            wm.index.name = 'datetime'
            confint.index.name = 'datetime'
            confint.columns.name = 'interval_boundary'
            
            ds = xr.Dataset()
            ds['V0']  = wm
            ds['confidence_interval'] = confint
            ds = ds.assign_coords({'channel': [wl,]})
            wm_list.append(ds)
            
            if smoothned:
                #### TODO: do i use the smoothening? if not make optional!
                assert(wm.isna().sum() == 0)
                # new_daterange = pd.date_range(start = wm.index[0], end = wm.index[-1], 
                #               # periods=1, 
                #               freq='D')
                spline = scint.UnivariateSpline(pltdates.date2num(wm.index), wm.values, 
                                              # s = wm.shape[0] * 20, 
                                              s = self.daterange2predict.shape[0] * 15,
                                             )
            
                wmsmooth = spline(pltdates.date2num(self.daterange2predict))
                wmsmooth = pd.Series(wmsmooth, index = self.daterange2predict)
                    
                spline = scint.UnivariateSpline(pltdates.date2num(std.index), std.values, 
                                              # s = wm.shape[0] * 20, 
                                              s = self.daterange2predict.shape[0] * 15,
                                             )
            
                stdsmooth = spline(pltdates.date2num(self.daterange2predict))
                stdsmooth = pd.Series(stdsmooth, index = self.daterange2predict)
            
            
                confint = pd.DataFrame({'low':wmsmooth-stdsmooth, 'high':wmsmooth+stdsmooth})
            
                wmsmooth.index.name = 'datetime'
                confint.index.name = 'datetime'
                confint.columns.name = 'interval_boundary'
            
                ds = xr.Dataset()
            
                ds['V0']  = wmsmooth
            
                ds['confidence_interval'] = confint
            
                ds = ds.assign_coords({'channel': [wl,]})
            
                ds_list.append(ds)
        

        if smoothned:
            ds_all = xr.concat(ds_list, 'channel')
            self._v0rol =   CalibrationPrediction(ds_all)
        else:
            wm_all = xr.concat(wm_list, 'channel')
            #### interpolate and extrapolate
            wm_all = wm_all.interp(datetime = self.daterange2predict,
                                     method='nearest', 
                                    kwargs={"fill_value": "extrapolate"})
           
            self._v0rol =   CalibrationPrediction(wm_all)
        
        return self._v0rol
    
    @property
    def V0_simple(self):
        #todo: expand this documentation
        """This simply returns the V0 based on all langley result in this object. Therefore, kick out what you don't want!"""
        dsout = xr.Dataset()
        v0 = np.exp(self.dataset.langley_fitres.sel(fit_results = 'intercept', drop = True))
        dsout['V0'] = v0.mean('datetime')
        dsout['V0_std'] = v0.std('datetime', ddof = 1) 
        dsout.V0_std.attrs['description'] = 'nbiased standard deviation, ddof = 1'
        ns = []
        osubs = []
        alpha = 0.05
        for wl in self.dataset.channel:
            n = self.dataset.langley_fitres.sel(fit_results = 'intercept', channel = wl).dropna('datetime').shape[0]
            ns.append(n)
            chi2val = scipy.stats.chi2.ppf(alpha, df=n-1)
            osub = np.sqrt((n-1)/chi2val) # this is the one-sided upper bound factor with a 95% confidence
            osubs.append(osub)
            
        dsout['number_of_cals'] = ('channel', ns)
        dsout.number_of_cals.attrs['description'] = 'Number of langleys used.'
        dsout['one_sided_upper_bound_factor_95conf'] = ('channel', osubs)
        dsout.one_sided_upper_bound_factor_95conf.attrs['description'] = r'one-sided upper bound factor with a 95% confidence, mostly relevant for small number of langleys.'


        additional_uncertainty = 0.005 # ARM introduced this additional uncertainty to cover remaining uncertainties. In Michalsky 2001 this is reflected in the remaining uncertainty even in the MLO calibrations.
        dsout['OD_uncertainty'] = (dsout['one_sided_upper_bound_factor_95conf'] * dsout['V0_std'] / dsout['V0']) + additional_uncertainty
        dsout.OD_uncertainty.attrs['description'] = r'(V0_std / V0 * osub) + 0.005. osub: one-sided upper bound factor with a 95% confidence, mostly relevant for small number of langleys. This still needs to be multiplide by 1/AMF to get the uncertainty in (A)OD. Michalsky 2001.'
        dsout['V0_stderr'] = self.dataset.langley_fitres.sel(fit_results = 'intercept_stderr', drop = True).mean('datetime')
        # dsout = dsout.drop_vars('quantile')
        return dsout

    def makev0prediction4date(self, date):

        where_ds = self.threshold_mask
        # todo: the below (and a few lines below that) is only when am and pm (time_of_day) are combined, there is probably interested in separating that out at some point?
        where_ds = where_ds.stack(obs = ("datetime", "time_of_day"), create_index=False).drop_vars("time_of_day").swap_dims({"obs": "datetime"})
        where = ~xr.zeros_like(where_ds.intercept, dtype = bool)
        for thd in where_ds:
            where &= where_ds[thd]

        intercept = self.dataset.fitresults.sel({'fitresult' : 'intercept',
                                    # 'channel': channel,
                                    # 'time_of_day': time_of_day
                                            }, drop = True)

        intercept = intercept.stack(obs = ("datetime", "time_of_day"), create_index=False).drop_vars("time_of_day").swap_dims({"obs": "datetime"})

        channel_l = []
        for ch in intercept.channel.values:
            # break
            interc = np.exp(intercept.where(where).sel(channel = ch, drop = True))
            # df = pd.DataFrame(interc).dropna()
            df = interc.drop_vars('quantile', errors='ignore').to_dataframe().dropna()
            prediction_date = pd.to_datetime(date)
            valid_dates = df.index.sort_values()
            before_index = valid_dates.searchsorted(prediction_date, side='right') - 1
            after_index = valid_dates.searchsorted(prediction_date, side='left')
            df.loc[prediction_date] = np.nan
            # df.loc[pd.to_datetime('2029-02-01')] = np.nan
            df = df.sort_index()
            days_to_next_valid_langley = np.full(2, np.nan)
            if before_index >= 0:
                days_to_next_valid_langley[0] = (
                    prediction_date - valid_dates[before_index]
                ) / pd.to_timedelta(1, 'day')
            if after_index < valid_dates.size:
                days_to_next_valid_langley[1] = (
                    valid_dates[after_index] - prediction_date
                ) / pd.to_timedelta(1, 'day')
        
            ts = atmts.TimeSeries(df)
            gam = ts.statistics.gamcl_r
            gam.robust = True
            # gam.
            # gam.rscript_path = "/Users/htelg/miniconda3/envs/py13/bin/Rscript"
            # gam.ncv = True
            # gam.ncv_window_days = 1
            # gam.ar1 = True
            # gam.newey_west = True
            gamall = gam
            ch_prediction = gam.prediction.sel(datetime_full = date, drop = True)
            ch_prediction['days_to_next_valid_langley'] = xr.DataArray(
                days_to_next_valid_langley,
                coords={'direction': ['before', 'after']},
                dims=('direction',),
                attrs={
                    'long_name': 'Days to nearest valid Langley calibration',
                    'description': (
                        'Distance from the prediction date to the nearest valid '
                        'Langley calibration on or before and on or after it. '
                        'The value is NaN when the prediction is extrapolated '
                        'in that direction.'
                    ),
                    'units': 'days',
                },
            )
            channel_l.append(dict(channel = ch, gam = gam, prediction = ch_prediction))
        self._prediction = channel_l
        return channel_l

    @property
    def v0prediction(self):
        """
        The idea is to have the programm to decide which of the two (gam or roll) to use

        Returns
        -------
        None.

        """
        if isinstance(self._v0, type(None)):
            if (self.dataset.datetime.values[-1] - self.dataset.datetime.values[0]) / pd.to_timedelta(365, 'days') < 2:
                self._v0 = self.v0prediction_rolling
            else:
                self._v0 = self.v0prediction_gam
        return self._v0
        
    @property
    def v0prediction_gam(self):
        if isinstance(self._v0gam, type(None)):
            self._get_v0_gam()
        return self._v0gam
            
    @property
    def v0prediction_rolling(self):
        if isinstance(self._v0rol, type(None)):
            self._get_v0_rolling()
        return self._v0rol

    def plot_threshold_effects(self):
        where_ds = self.threshold_mask
        # todo: the below (and a few lines below that) is only when am and pm (time_of_day) are combined, there is probably interested in separating that out at some point?
        where_ds = where_ds.stack(obs = ("datetime", "time_of_day"), create_index=False).drop_vars("time_of_day").swap_dims({"obs": "datetime"})
        intercept = self.dataset.fitresults.sel({'fitresult' : 'intercept',
                                    # 'channel': channel,
                                    # 'time_of_day': time_of_day
                                            }, drop = True)

        intercept = intercept.stack(obs = ("datetime", "time_of_day"), create_index=False).drop_vars("time_of_day").swap_dims({"obs": "datetime"})
        for ch in where_ds.channel.values:
            intercept_sel = intercept.sel(channel = ch, drop = True)
            where_ds_sel = where_ds.sel(channel = ch, drop = True)
            lim = (0,15)
            # lim = (6,9)
            # lim = (7,8)
            # center = intercept_sel
            # width = 0.5
            # lim = (center-width, center+width)
            
            
            f,aa = plt.subplots(2)
            f.set_figheight(f.get_figheight() * 1.5)
            bins = np.linspace(lim[0], lim[1], 30)
            aa[0].hist(intercept_sel, bins = bins)
            intercept_sel.plot(marker = '.', ls = '', ax = aa[1])
            
            
            where = ~xr.zeros_like(where_ds_sel.datetime, dtype = bool)
            for thd in where_ds_sel:
                where &= where_ds_sel[thd]
                aa[0].hist(intercept_sel.where(where), bins = bins, label = f'{thd}')
                intercept_sel.where(where).plot(marker = '.', ls = '', ax = aa[1])
            
                
            a = aa[0]
            a.set_xlim(*lim)
            a.legend()
            a.set_title(f'channel = {ch} nm')
            
            a = aa[1]
            a.set_ylim(*lim)
            
            g = a.get_lines()[-1]
            g.set_color('black')
            
            no_month = ((intercept_sel.datetime.max() - intercept_sel.datetime.min()) / pd.to_timedelta(30, 'D'))
            points_per_month = where.sum()/ no_month # points per month
            a.text(0.05, 0.85, f'{points_per_month:0.2f}', transform = a.transAxes)

    def plot_gam_results(self):
        if self._prediction is None:
            raise ValueError('Make a prediction first. Use makev0prediction4date')
        for ch in self._prediction:
            out = ch['gam'].plot_overview()
            f,aa = out
            aa[0].set_title(ch['channel'])
        
    
def fit_langley(langley, robust = True,
                # airmasslimits = (2.5,4),
                weighted = True, error_handling = 'skip', parent = None,
                weights = None):
    """
    Performs a linear fit to the langleys. 
    Parameters
    ----------
    langley : TYPE
        DESCRIPTION.
    weighted : bool, optional
        If the fit is to be weighted. This is a good idea due to the change in density of data points as the AMF decreases. The default is ....
    weights : pandas.Series or xarray.DataArray, optional
        Fixed fit weights indexed by airmass. When omitted, weights are
        calculated from the spacing of the points supplied to this function.

    Returns
    -------
    None.

    """
    langley = langley.copy()
    if weighted and weights is not None:
        if isinstance(weights, xr.DataArray):
            weights = weights.to_series()
        elif not isinstance(weights, pd.Series):
            weights = pd.Series(np.asarray(weights, dtype=float), index=langley.index)
    langley = langley.replace([np.inf, -np.inf], np.nan).dropna()    
    # langley = langley.dropna() # this is necessary as otherwise the fit fails
    langley.sort_index(ascending=False, inplace=True) # otherwise the weigted fit will fail for pm
    # langley = langley.truncate(*airmasslimits)
    if not (langley.shape[0] > 1):
        # f'Not enough data points to fit langley. Only {langley.shape[0]} points available.'
        lser = {'slope': np.nan, 'intercept': np.nan, 'slope_stderr': np.nan, 'intercept_stderr': np.nan, 'resid_curvature': np.nan}
        fitfailed = True
        res = False
        out = {}
        out['res_series'] = pd.Series(lser)
        out['res'] = res
        return out

    
    y = langley
    x = langley.index
    x = statsmodels.api.add_constant(x)
    if weighted:
        if weights is None:
            idx = langley.index.to_numpy(dtype=float)
            w = np.empty_like(idx)
            w[0] = (idx[0] - idx[1]) / 2
            w[-1] = (idx[-2] - idx[-1]) / 2
            w[1:-1] = (idx[:-2] - idx[2:]) / 2
        else:
            fit_weights = weights.reindex(langley.index)
            w = fit_weights.to_numpy(dtype=float, copy=True)
            if np.any(~np.isfinite(w)) or np.any(w <= 0):
                raise ValueError('Fit weights must be finite and greater than zero for every retained point.')
            w /= w.mean()
        if parent is not None:
            parent.tp_w = w
            parent.tp_langley = langley

    try:
        if robust and langley.shape[0] > x.shape[1]:
            if weighted:
                y_values = y.to_numpy(dtype=float)
                combined_weights = w
                xtw = x.T * combined_weights
                params = np.linalg.solve(xtw @ x, xtw @ y_values)
                huber = statsmodels.api.robust.norms.HuberT()

                for _ in range(50):
                    resid = y_values - x @ params
                    scale = statsmodels.api.robust.scale.mad(resid)
                    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
                        break
                    combined_weights = w * huber.weights(resid / scale)
                    xtw = x.T * combined_weights
                    next_params = np.linalg.solve(xtw @ x, xtw @ y_values)
                    if np.max(np.abs(next_params - params)) <= 1e-8 * max(1.0, np.max(np.abs(params))):
                        params = next_params
                        break
                    params = next_params

                mod = statsmodels.api.WLS(y, x, weights=combined_weights)
            else:
                mod = statsmodels.api.RLM(y, x)
        else:
            if weighted:
                mod = statsmodels.api.WLS(y, x, weights=w)
            else:
                mod = statsmodels.api.OLS(y, x)

        if robust and weighted and langley.shape[0] > x.shape[1]:
            res = mod.fit(cov_type='HC1')
        else:
            res = mod.fit()
        assert(~np.isnan(res.params['x1'])), 'The fit failed and returned a nan. This should not happen.'
        lser = {'slope': res.params['x1'], 'intercept': res.params['const'], 'slope_stderr': res.bse['x1'], 'intercept_stderr': res.bse['const']}
        fitfailed = False
        
    except np.linalg.LinAlgError:
        if error_handling == 'raise':
            raise
        lser = {'slope': np.nan, 'intercept': np.nan, 'slope_stderr': np.nan, 'intercept_stderr': np.nan}
        fitfailed = True
        res = False
        
    if not fitfailed:
        ### test for curved residual
        x = res.resid.index.values
        y = res.resid.values
        pfres = np.polyfit(x,y,2)    
        lser['resid_curvature'] =  pfres[0]
    else:
        lser['resid_curvature'] =  np.nan

    lser=pd.Series(lser)
    out = {}
    out['res_series'] = lser
    out['res'] = res
    return out

def plot_langley(lang_ds, wls = 'all', date = 'Test'):
    if wls == 'all':
        wls = lang_ds.channel.values
    else:
        pass
#     wl = 500
    aas = []
    if 'datetime' in lang_ds.coords:
        lang_ds = lang_ds.isel(datetime = 0)
    for wl in wls:
        lan = lang_ds.langleys.to_pandas()[wl]
        resid = lang_ds.langley_fit_residual.to_pandas()[wl]
        f,aa = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0})
        aas.append(aa)
        alan, ares = aa
        
        fitres = lang_ds.langley_fitres.to_pandas().loc[wl]
        fit = fitres.intercept + fitres.slope * lang_ds.airmass.to_pandas()
        fit.plot(ax = alan)
        lan.plot(ax = alan, ls = '', marker = '.')
        g = alan.get_lines()[-1]
        g.set_markersize(2)
        
        fitres['corrmatrixdet'] = lang_ds.langley_residual_correlation_prop.values
        
        alan.text(0.01,0.05, fitres.__str__(), transform=alan.transAxes, fontsize = 'x-small')

        resid.plot(ax = ares)
        ares.grid()
        
        alan.set_title(f'{date} - channel {wl}nm')
    return f,aas

class Langley(object):
    def __init__(self, parent, langleys, airmass_limits = (2.2,4.7),
                 langley_fit_settings = None, when = None):
        self.parent = parent
        amf_min, amf_max = airmass_limits
        ds = langleys.where((langleys.airmass > amf_min) & (langleys.airmass < amf_max), drop = True)
        self.dataset = ds
        langleys = ds.direct_normal.to_pandas()
        langleys.columns.name = 'channel'
        # langleys = langleys.truncate(*airmass_limits)
        # langleys = langleys.copy()
        # langleys = langleys.where(aimass_limits[0] < langleys.airmass).where(aimass_limits[1] > langleys.airmass).dropna('airmass')
        self.langleys = langleys
        self.langley_fit_settings = langley_fit_settings
        self.refresh()
        self.langley_pre_clean = False
        self.when = when
    
    def refresh(self):
        self._langley_fitres = None
        self._langley_residual_correlation_prop = None
        self._langley_fit_residual = None

    @property
    def expected_airmass_cdf(self):
        """Expected airmass CDF for uniform temporal sampling."""
        if 'expected_airmass_cdf' not in self.dataset:
            ds = self.dataset.sortby('airmass')
            count = ds.sizes['airmass']
            if count == 0:
                cdf = xr.DataArray(
                    np.empty(0, dtype=float),
                    coords={'airmass': ds.airmass},
                    dims=('airmass',),
                )
            elif count == 1:
                cdf = xr.zeros_like(ds.airmass, dtype=float)
            else:
                if bool(ds.datetime.isnull().any()):
                    raise ValueError('Langley datetime values must not contain missing values.')
                timestamp = ds.datetime.astype('datetime64[ns]').astype('int64')
                time_span = timestamp.isel(airmass=-1) - timestamp.isel(airmass=0)
                if time_span.item() == 0:
                    raise ValueError('Langley datetime values must span more than one instant.')
                cdf = (timestamp - timestamp.isel(airmass=0)) / time_span
                if bool((cdf.diff('airmass') <= 0).any()):
                    raise ValueError('Langley datetime must vary monotonically with airmass.')

            cdf = cdf.rename('expected_airmass_cdf')
            cdf.attrs = {
                'long_name': 'expected cumulative distribution of airmass',
                'description': 'Normalized elapsed time, oriented from low to high airmass.',
                'units': '1',
            }
            self.dataset['expected_airmass_cdf'] = cdf
        return self.dataset.expected_airmass_cdf.sortby('airmass')

    @property
    def airmass_weights(self):
        """Fixed inverse-density weights on the original airmass sequence."""
        if 'airmass_weights' not in self.dataset:
            cdf = self.expected_airmass_cdf
            airmass = np.asarray(cdf.airmass.data, dtype=float)
            expected_quantile = np.asarray(cdf.data, dtype=float)
            count = airmass.size

            weights = np.ones(count, dtype=float)
            if count > 1:
                weights[0] = ((airmass[1] - airmass[0]) /
                              (expected_quantile[1] - expected_quantile[0]))
                weights[-1] = ((airmass[-1] - airmass[-2]) /
                               (expected_quantile[-1] - expected_quantile[-2]))
                if count > 2:
                    weights[1:-1] = ((airmass[2:] - airmass[:-2]) /
                                     (expected_quantile[2:] - expected_quantile[:-2]))

            if np.any(~np.isfinite(weights)) or np.any(weights <= 0):
                raise ValueError('Expected airmass weights must be finite and greater than zero.')
            if count:
                weights /= weights.mean()

            self.dataset['airmass_weights'] = xr.DataArray(
                weights,
                coords={'airmass': cdf.airmass},
                dims=('airmass',),
                name='airmass_weights',
                attrs={
                    'long_name': 'fixed Langley airmass fit weights',
                    'description': 'Inverse expected sampling density, dm/du, normalized to mean one.',
                    'units': '1',
                },
            )
        return self.dataset.airmass_weights.sortby('airmass')

    @property
    def airmass_coverage_metrics(self):
        """Coverage of retained observations in expected-CDF coordinates."""
        expected_cdf = self.expected_airmass_cdf.to_pandas()
        expected_weights = self.airmass_weights.to_pandas()
        channels = self.langleys.columns.to_numpy()
        channel_count = channels.size
        retained_points = np.zeros(channel_count, dtype=int)
        retained_fraction = np.full(channel_count, np.nan)
        cdf_span = np.full(channel_count, np.nan)
        maximum_cdf_gap = np.full(channel_count, np.nan)
        maximum_cdf_gap_ratio = np.full(channel_count, np.nan)
        cdf_spacing_cv = np.full(channel_count, np.nan)
        effective_airmass = np.full(channel_count, np.nan)
        effective_airmass_offset = np.full(channel_count, np.nan)

        if expected_cdf.size > 1:
            reference_spacing = np.median(np.diff(expected_cdf.to_numpy(dtype=float)))
        else:
            reference_spacing = np.nan
        if expected_cdf.size:
            airmass_min = float(expected_cdf.index.min())
            airmass_max = float(expected_cdf.index.max())
            airmass_midpoint = (airmass_min + airmass_max) / 2
            airmass_half_range = (airmass_max - airmass_min) / 2
        else:
            airmass_midpoint = np.nan
            airmass_half_range = np.nan

        for position, channel in enumerate(channels):
            retained_airmass = self.langleys[channel].dropna().index
            retained_cdf = expected_cdf.reindex(retained_airmass).dropna().sort_values().to_numpy(dtype=float)
            retained_weights = expected_weights.reindex(retained_airmass).dropna()
            retained_points[position] = retained_cdf.size
            if expected_cdf.size:
                retained_fraction[position] = retained_cdf.size / expected_cdf.size
            if retained_weights.size:
                effective_airmass[position] = np.average(
                    retained_weights.index.to_numpy(dtype=float),
                    weights=retained_weights.to_numpy(dtype=float),
                )
                if airmass_half_range > 0:
                    effective_airmass_offset[position] = (
                        (effective_airmass[position] - airmass_midpoint) / airmass_half_range
                    )
            if retained_cdf.size:
                cdf_span[position] = retained_cdf[-1] - retained_cdf[0]
                coverage_spacing = np.diff(np.concatenate(([0.0], retained_cdf, [1.0])))
                maximum_cdf_gap[position] = coverage_spacing.max()
                maximum_cdf_gap_ratio[position] = coverage_spacing.max() / reference_spacing
            elif expected_cdf.size:
                maximum_cdf_gap[position] = 1.0
                maximum_cdf_gap_ratio[position] = 1.0 / reference_spacing
            if retained_cdf.size > 1:
                spacing = np.diff(retained_cdf)
                cdf_spacing_cv[position] = spacing.std() / spacing.mean()

        metrics = xr.Dataset(
            data_vars={
                'expected_points': (
                    'channel',
                    np.full(channel_count, expected_cdf.size, dtype=int),
                ),
                'airmass_midpoint': (
                    'channel',
                    np.full(channel_count, airmass_midpoint),
                ),
                'airmass_half_range': (
                    'channel',
                    np.full(channel_count, airmass_half_range),
                ),
                'retained_points': ('channel', retained_points),
                'retained_fraction': ('channel', retained_fraction),
                'cdf_span': ('channel', cdf_span),
                'maximum_cdf_gap': ('channel', maximum_cdf_gap),
                'maximum_cdf_gap_ratio': ('channel', maximum_cdf_gap_ratio),
                'cdf_spacing_cv': ('channel', cdf_spacing_cv),
                'effective_airmass': ('channel', effective_airmass),
                'effective_airmass_offset': ('channel', effective_airmass_offset),
            },
            coords={'channel': self.dataset.channel},
        )
        metrics.expected_points.attrs = {'long_name': 'number of expected points', 'units': '1'}
        metrics.airmass_midpoint.attrs = {'long_name': 'midpoint of the expected airmass range', 'units': '1'}
        metrics.airmass_half_range.attrs = {'long_name': 'half-width of the expected airmass range', 'units': '1'}
        metrics.retained_points.attrs = {'long_name': 'number of retained points', 'units': '1'}
        metrics.retained_fraction.attrs = {'long_name': 'fraction of expected points retained', 'units': '1'}
        metrics.cdf_span.attrs = {'long_name': 'span of retained expected-airmass CDF', 'units': '1'}
        metrics.maximum_cdf_gap.attrs = {'long_name': 'maximum uncovered expected-airmass CDF interval', 'units': '1'}
        metrics.maximum_cdf_gap_ratio.attrs = {'long_name': 'maximum CDF gap relative to expected spacing', 'units': '1'}
        metrics.cdf_spacing_cv.attrs = {'long_name': 'coefficient of variation of retained CDF spacing', 'units': '1'}
        metrics.effective_airmass.attrs = {'long_name': 'weighted mean airmass of retained points', 'units': '1'}
        metrics.effective_airmass_offset.attrs = {
            'long_name': 'effective airmass offset from the range midpoint, normalized by the half-range',
            'units': '1',
        }
        return metrics

    
    def to_xarray_dataset(self):
        if "sn_mfrsr" in self.parent.dataset.attrs:
            serialnoattr = 'sn_mfrsr'
        elif "serial_no" in self.parent.dataset.attrs:
            serialnoattr = 'serial_no'
        else:
            raise KeyError('No serial number found in dataset attributes. Expected either "sn_mfrsr" or "serial_no".')
        serialno = self.parent.dataset.attrs[serialnoattr]
        try:
            serialno = serialno.values[0]
        except IndexError:
            serialno = serialno.values
        except AttributeError:
            serialno = serialno
        ds = xr.Dataset({'langleys': self.langleys,
                         'langley_fit_residual': self.langley_fit_residual,
                         'langley_fitres': self.langley_fitres,
                         'langley_residual_correlation_prop': self.langley_residual_correlation_prop['determinant'],
                         'sp02_serial_no': serialno},)
        ds.attrs = dict(when = self.when,
                        date= pd.to_datetime(self.parent.dataset.datetime.values[0]).strftime('%Y%m%d'))

        ds['airmass_coverage_metrics'] = self.airmass_coverage_metrics.to_dataarray()        

        return ds
    
    def save2netcdf(self, path2file):
        ds = self.to_xarray_dataset()
        ds.to_netcdf(path2file)
        return ds

    def save2database(self, langley_results_db, path2file = None, version = None):
        """
        Parameters
        ----------
        langley_results_db: str or pathlib.Path
            Path to database file. If it does not exist, it will be created.
            """
        database = atmlldb.LangleyResultDatabase(
            langley_results_db,
            create_if_not_exist=True,
            version=version
        )
        database.write_langley(
            date=self.parent.dataset.datetime.values[0],
            path2file=path2file,
            langley=self,
        )
        return database
    
    @property
    def langley_fit_residual(self):
        if isinstance(self._langley_fit_residual, type(None)):
            self._fit_langles()
        return self._langley_fit_residual
    
    @property
    def langley_fitres(self):
        if isinstance(self._langley_fitres, type(None)):
            self._fit_langles()
        return self._langley_fitres
    
    @property
    def langley_residual_correlation_prop(self):
        if isinstance(self._langley_residual_correlation_prop, type(None)):
            self._fit_langles()
        return self._langley_residual_correlation_prop
    
    def _fit_langles(self, verbose = False, error_handling = 'raise'):
        langleys = self.langleys
        airmass_weights = self.airmass_weights.to_pandas()
        # df_langres = pd.DataFrame(index=['slope', 'intercept', 'stderr'])
        fit_res_dict = {}
        # resid_dict = {}
        df_resid = pd.DataFrame(index=langleys.index)
        for wl in langleys:
            # lrg_res = stats.linregress(langleys.index, langleys[wl])
            # df_langres[wl] = pd.Series({'slope': lrg_res.slope, 'intercept': lrg_res.intercept, 'stderr': lrg_res.stderr})
            out =  fit_langley(langleys[wl], error_handling=error_handling, parent = self,
                                   weights=airmass_weights)
            # df_langres[wl] = out['res_series']
            fit_res_dict[wl] = out['res_series']
            if out['res']:
                df_resid[wl] = out['res'].resid
                # resid_dict[wl] = out['res'].resid
            else:
                # resid_dict[wl] = np.nan
                df_resid[wl] = np.nan
        
        df_langres = pd.DataFrame(fit_res_dict)
        df_langres = df_langres.transpose()
        df_langres.index.name = 'channel'
        df_langres.columns.name = 'fit_results'
        self._langley_fitres = df_langres
        
        # df_resid = pd.DataFrame(resid_dict)
        df_resid.columns.name = 'channel'
        
        resid_corr_matrix = df_resid.corr()
        lrcp = {}
        lrcp['determinant'] = np.linalg.det(resid_corr_matrix)
        
        self._langley_fit_residual = df_resid
        self._langley_residual_correlation_prop = lrcp
        return df_langres
    
    def clean(self, threshold = 3, use_channels = None, clean_channels_individually = True, verbose = False,
              *, expected_instrument_noise = 0.01, min_channels = 2, max_iterations = 3):
        """Remove points using channel-scaled, leverage-adjusted residuals.

        Parameters
        ----------
        threshold : float or mapping
            Dimensionless residual threshold. A mapping may specify a different
            threshold for each selected channel.
        use_channels : sequence, optional
            Channels used for screening. By default all channels except 940 nm
            are used.
        clean_channels_individually : bool
            Mask only the flagged channel when true. Otherwise, remove a complete
            airmass row when at least ``min_channels`` channels flag it. While a 
            nice idea, this is not recommended because if a single channel
            malfunctions, it will effect all other channels and might not be detected
            as a problem.
        verbose : bool
            Print the number of newly flagged points after each fit.
        expected_instrument_noise : float or mapping, optional
            Lower bound for each channel's residual scale. Values are in log-signal
            units and may be specified independently for each channel.
        min_channels : int
            Number of selected channels that must flag a point before its complete
            airmass row is removed. Ignored when channels are cleaned individually.
        max_iterations : int
            Maximum number of removal and refit passes.

        Returns
        -------
        dict
            The cleaned Langley and diagnostics. ``removed_points`` is a boolean
            mask whose true values identify removed rows or channel values.
        """
        if isinstance(max_iterations, bool) or not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
            raise ValueError('max_iterations must be a positive integer.')

        langley = copy.deepcopy(self)
        langley.langley_pre_clean = self

        if use_channels is None:
            channels = [channel for channel in langley.langleys.columns if channel != 940]
        elif np.isscalar(use_channels):
            channels = [use_channels]
        else:
            channels = list(use_channels)
        missing_channels = [channel for channel in channels if channel not in langley.langleys.columns]
        if missing_channels:
            raise KeyError(f'Channels not found in Langley data: {missing_channels}')
        if not channels:
            raise ValueError('No channels available for cleaning.')
        if not clean_channels_individually:
            if isinstance(min_channels, bool) or not isinstance(min_channels, (int, np.integer)) or min_channels < 1:
                raise ValueError('min_channels must be a positive integer.')
            if len(channels) == 1:
                min_channels = 1
            elif min_channels > len(channels):
                raise ValueError('min_channels cannot exceed the number of selected channels.')

        if np.isscalar(threshold):
            thresholds = pd.Series(float(threshold), index=channels)
        else:
            thresholds = pd.Series(threshold, dtype=float).reindex(channels)
            if thresholds.isna().any():
                missing = thresholds.index[thresholds.isna()].to_list()
                raise KeyError(f'Missing thresholds for channels: {missing}')
        if np.any(~np.isfinite(thresholds)) or np.any(thresholds <= 0):
            raise ValueError('All thresholds must be finite and greater than zero.')

        if expected_instrument_noise is None:
            instrument_noise = pd.Series(0.0, index=channels)
        elif np.isscalar(expected_instrument_noise):
            instrument_noise = pd.Series(float(expected_instrument_noise), index=channels)
        else:
            instrument_noise = pd.Series(expected_instrument_noise, dtype=float).reindex(channels)
            if instrument_noise.isna().any():
                missing = instrument_noise.index[instrument_noise.isna()].to_list()
                raise KeyError(f'Missing expected instrument noise for channels: {missing}')
        if np.any(~np.isfinite(instrument_noise)) or np.any(instrument_noise < 0):
            raise ValueError('Expected instrument noise must be finite and non-negative.')

        if clean_channels_individually:
            removed_points = pd.DataFrame(False, index=self.langleys.index, columns=self.langleys.columns)
        else:
            removed_points = pd.Series(False, index=self.langleys.index, name='removed')

        status = 'maximum iterations reached'
        removal_iterations = 0
        standardized_residuals = pd.DataFrame(index=self.langleys.index, columns=channels, dtype=float)
        residual_scale = pd.Series(np.nan, index=channels, dtype=float)
        reference_weights = langley.airmass_weights.to_pandas()

        for _ in range(max_iterations):
            if verbose:
                print(f'cleaning iteration {removal_iterations}')

            lfr = langley.langley_fit_residual.loc[:, channels]
            if lfr.dropna(how='all').empty:
                status = 'fit failed'
                break
            langley.tp_lfr = lfr.copy()

            center = lfr.median()
            mad = lfr.subtract(center, axis='columns').abs().median() * 1.4826
            numerical_scale = pd.Series(np.sqrt(np.finfo(float).eps), index=channels) # the only reason for this is to avoid division by zero. It provides machine precision, but is not a realistic scale. 
            residual_scale = pd.concat([mad, instrument_noise, numerical_scale], axis=1).max(axis=1)

            leverage = pd.DataFrame(np.nan, index=lfr.index, columns=channels)
            for channel in channels:
                channel_residual = lfr[channel].dropna().sort_index(ascending=False)
                if channel_residual.shape[0] <= 2:
                    continue
                idx = channel_residual.index.to_numpy(dtype=float)
                weights = reference_weights.reindex(channel_residual.index).to_numpy(dtype=float, copy=True)
                weights /= weights.mean()
                design = np.column_stack([np.ones_like(idx), idx])
                covariance = np.linalg.pinv(design.T @ (weights[:, np.newaxis] * design))
                channel_leverage = weights * np.einsum('ij,jk,ik->i', design, covariance, design)
                leverage.loc[channel_residual.index, channel] = np.clip(channel_leverage, 0, 1)

            leverage_adjustment = np.sqrt((1 - leverage).clip(lower=np.finfo(float).eps))
            denominator = leverage_adjustment.multiply(residual_scale.replace(0, np.nan), axis='columns')
            iteration_standardized_residuals = lfr.subtract(center, axis='columns').abs().divide(denominator)
            standardized_residuals.update(iteration_standardized_residuals)
            threshold_ratio = iteration_standardized_residuals.divide(thresholds, axis='columns')

            if clean_channels_individually:
                new_removed = threshold_ratio.gt(1) & lfr.notna()
            else:
                outlier_channels = threshold_ratio.gt(1) & lfr.notna()
                new_removed = outlier_channels.sum(axis=1).ge(min_channels)

            if verbose:
                print(new_removed.sum().to_list() if clean_channels_individually else int(new_removed.sum()))

            if not np.any(new_removed):
                status = 'converged'
                break

            # self.tp_cutoff = thresholds * residual_scale
            # self.tp_where = ~new_removed
            # self.tp_langley = langley.langleys.copy()

            if clean_channels_individually:
                langley.langleys.loc[:, channels] = langley.langleys.loc[:, channels].mask(new_removed)
                removed_points.loc[:, channels] = removed_points.loc[:, channels] | new_removed.reindex(
                    removed_points.index, fill_value=False
                )
            else:
                langley.langleys = langley.langleys.loc[~new_removed]
                removed_points = removed_points | new_removed.reindex(removed_points.index, fill_value=False)

            removal_iterations += 1
            langley.refresh()

        out = {'langley': langley}
        out['iterations'] = removal_iterations
        out['status'] = status
        out['removed_points'] = removed_points
        out['standardized_residuals'] = standardized_residuals
        out['residual_scale'] = residual_scale
        return out
    
    def plot(self, channel = None, show_pre_clean = True, plot_residual = True, textpos = [0.1, 0.1], ax = None):
        def plot_fit(wl, ax = None):
            res = self.langley_fitres.loc[wl]
            fit_index = self.langleys.index.append(pd.Index([0])).sort_values()
            fit = pd.DataFrame(res.intercept + (res.slope * fit_index), index = fit_index)#, columns=['fit'])
            fit.columns = ['fit',]
            if not isinstance(self.langley_pre_clean, (type(None), bool)) and show_pre_clean: # if there was a cleaning step, plot the original data as well
                self.langley_pre_clean.langleys[wl].plot(ax = ax, marker = '.', ls = '', markersize = 3, 
                                                         color = 'red', 
                                                        #  alpha = 0.5
                                                         )
            self.langleys[wl].plot(ax = ax, marker = '.', ls = '', markersize = 5)

            fit.plot(ax = ax)
            fr = self.langley_fitres.loc[wl]
            txt = f'channel: {wl} nm\n'
            txt += f'slope: {fr.slope:0.2g} ± {fr.slope_stderr:0.1g}\n'
            txt += f'intercept: {fr.intercept:0.4f} ± {fr.intercept_stderr:0.4f}\n'
            txt += f'no points: {self.langleys[wl].dropna().shape[0]}'
            ax.text(*textpos, txt, transform = ax.transAxes)
            return None
        def plot_residual(wl, ax = None):
            res = self.langley_fit_residual[wl]
            res.plot(ax = ax, marker = '.', ls = '', markersize = 5, zorder = 10, label = 'cleaned')
            g = ax.get_lines()[-1]
            if not isinstance(self.langley_pre_clean, (type(None), bool)) and show_pre_clean: # if there was a cleaning step, plot the original data as well
                self.langley_pre_clean.langley_fit_residual[wl].plot(ax = ax, marker = '.', ls = '',
                                                                      markersize = g.get_markersize()*0.8, zorder = 5, label = 'pre-cleaned',
                                                        #  color = g.get_color(), 
                                                        #  alpha = 0.5
                                                         )
                gr = ax.get_lines()[-1]
                # gr.set_color(g.get_color())
                # gr.set_markerfacecolor('red')
                gr.set_color('red')
                gr.set_markerfacecolor(g.get_color())
                gr.set_markeredgewidth(0.5)

            ax.axhline(0, color = 'black', ls = '--')
            ax.set_ylabel('residual')
            ax.set_xlabel('airmass')
            gs = ax.get_lines()
            ymax = max([np.nanmax(g._y) for g in gs if g.get_label() == 'cleaned'])
            ymin = min([np.nanmin(g._y) for g in gs if g.get_label() == 'cleaned'])
            pad = abs(ymax - ymin) * 0.2
            ax.set_ylim(ymin - pad, ymax + pad)
            return None
        if isinstance(channel, type(None)):
            channel = self.langleys.columns
            # for wl in self.langleys.columns:
            #     f,a = plt.subplots()
            #     plot_fit(wl, ax)
        elif np.isscalar(channel):
            channel = [channel,]
        axout = []
        for e,wl in enumerate(channel):
            if isinstance(ax, type(None)):
                if plot_residual:
                    f,a = plt.subplots(2, sharex=True, height_ratios=[2,1], gridspec_kw={'hspace': 0})
                    f.set_figheight(f.get_figheight() * 1.5)
                else:
                    f,a = plt.subplots()
                axout.append(a)
            else:
                a = ax[e] if isinstance(ax, (list, np.ndarray)) else ax
                # f = ax.figure
            if plot_residual:
                plot_fit(wl, ax = a[0])
                plot_residual(wl, ax = a[1])
            else:
                plot_fit(wl, ax = a)
        return axout
    

class Calibration_SP02(object):
    def __init__(self, raw2lang_output):
        """I believe this is an old SP02 calibration object. not sure if used anymore"""
        self.raw2lang_output = raw2lang_output
        self.langley_fitres = raw2lang_output['ds_results']
        
    
#     sp02.calibration.get_best_results#, save_as_result = save)

    def plot_od(self, sn, ax = None, show_molecular = False):
        best_10 = self
        fitres_mean = best_10.langley_fitres.langley_fitres.mean(dim = 'datetime').to_pandas()
        wl = best_10.langley_fitres.channle_wavelengths
        if len(wl.dims) == 1:
            wl = wl.to_pandas()
        else:
            wl = wl[:,0].to_pandas()
        fitres_mean.rename(wl, inplace=True)
        fitres_mean['slope'] = fitres_mean.slope.abs()
    
        a = fitres_mean.slope.plot(ls = '', marker = 'o', ax = ax, label = f'sn{sn}')
        a.set_yscale('log')
        a.set_ylabel('OD')
        a.set_xscale('log')

        if show_molecular:
            x = np.linspace(368, 1050, 5)
            y = 1e10 * x**-4
            g, = plt.plot(x, y)
            g.set_zorder(0)
            g.set_linestyle('--')
            g.set_color('0.6')
            g.set_dash_capstyle('round')
            g.set_label('molecular')
        a.legend()  
        return a
        
    def plot_langley_results(self, intensity4ts = 'resid_curvature', add2pptx = False,test = False, tollerance = [2.5, 2.5]):
        """
        slope 	intercept 	slope_stderr 	intercept_stderr 	resid_curvature

        Parameters
        ----------
        result_ds : TYPE
            DESCRIPTION.
        intensity4ts : TYPE, optional
            DESCRIPTION. The default is ''.
        add2pptx : TYPE, optional
            DESCRIPTION. The default is False.
        test : TYPE, optional
            DESCRIPTION. The default is False.
        tollerance : TYPE, optional
            DESCRIPTION. The default is [2.5, 2.5].

        Returns
        -------
        out : TYPE
            DESCRIPTION.

        """
        ds = self.langley_fitres
        out = {}
        if add2pptx:
            date = pd.to_datetime(datetime.datetime.now()).date()
            date = str(date)
            date = date.replace('-','')
            path2ppt = f'/mnt/telg/projects/sp02/calibration/langley_calibration_summary_{date}.ppt'
            path2ppt
            prs = pptx.Presentation()
            title_slide_layout = prs.slide_layouts[0]
            slide = prs.slides.add_slide(title_slide_layout)
            title = slide.shapes.title
            subtitle = slide.placeholders[1]

            title.text = "SP02 revival"
            now = pd.to_datetime(datetime.datetime.now())

            subtitle.text = f"created on {now.date()} by Hagen"



        plt.rcParams['hatch.linewidth'] = 4
        for wl in ds.channel.values:
            f,a = plt.subplots()
            res = ds.langley_fitres.sel(channel = wl).to_pandas()

            rest = res.copy()

            restdna = rest.dropna()

            cm = plt.cm.plasma_r
            pc = plt.scatter(restdna.slope,restdna.intercept_stderr, c = np.e**restdna.intercept, s = 40, cmap = cm, 
                       )
            cb = f.colorbar(pc)
            cb.set_label('Langley fit - intercept')
            zm = (np.e**restdna.intercept).median()
            zstd  = (np.e**restdna.intercept).std()
            cb.ax.axhline(zm, color = 'black')
            hs = cb.ax.axhspan(zm-zstd, zm+zstd)
            hs.set_facecolor([0,0,0,0])
            hs.set_edgecolor([0,0,0,0.5])
            hatch = hs.set_hatch('//')

            a.set_title(f'channel: {wl}nm')
            a.set_xlabel('Langley fit - slope (mV)')
            a.set_ylabel('Langley fit - std of residual')

        # langley over time
            f_ivt,a = plt.subplots()
            out['tp_restdna'] = restdna
            restdna.intercept.plot(ax = a)
            pc = a.scatter(restdna.index, np.e**restdna.intercept, c = restdna[intensity4ts], cmap=plt.cm.gnuplot)
            f_ivt.colorbar(pc)
            a.set_xlabel('')
            a.set_ylabel('Langley fit - intercept')
            a.set_title('Langley fit intercept as a function of time')

            if add2pptx:

                blank_slide_layout = prs.slide_layouts[6]
                slide = prs.slides.add_slide(blank_slide_layout)
                add_fig2slide(f, slide)
                add_fig2slide(f_ivt, slide, top = 7.5/2)
#             if test:
            break
        if add2pptx:
            prs.save(path2ppt)
        return out
    
    def get_best_results(self, top = 10):
        ds = self.langley_fitres
        def sortbyparam(ds, param, top):
            df = (ds.langley_fitres.sel(linreg_param = f'{param}_stderr') / ds.langley_fitres.sel(linreg_param = param)).mean(dim = 'channel').to_pandas()
            df = df[df>0]
            df.sort_values(inplace=True)
            out = {}
            df = df.iloc[:top]
            out = ds.sel(datetime =  df.index.values)
            return out
        out = {}
        slopestdrmin = sortbyparam(ds,'slope', top)
        interceptstdrmin = sortbyparam(ds,'intercept',top)
        curvature_date = ds.langley_fitres.sel(linreg_param = 'resid_curvature').mean(dim = 'channel').to_pandas().abs().sort_values().index[:top].values
        curvature = ds.sel(datetime = curvature_date)
        df_cmp = ds.correlation_matrix_properties.to_pandas()
        determimin_dates = df_cmp.sort_values('determinant', ascending=False).index[:top].values
        determimin = ds.sel(datetime = determimin_dates)
    #     uberlap = ds.datetime.values
    #     for other in [interceptstdrmin, 
    #                  curvature, 
    #     #              determimin,
    #                  ]:
    #         uberlap = np.intersect1d(uberlap, other)
    #     out['crossection'] = uberlap

        out['slope_stderr'] = slopestdrmin
        out['intercept_stderr'] = interceptstdrmin
        out['curvature'] = curvature
        out['corr_matrix_determinant_max'] = determimin
        out['ds_results'] = interceptstdrmin
        return Calibration(out)
    
    def plot_good_enough(self, best=None, col = 'A', std_error_goal = [2, 2e-4, 2e-3, 1e-1]):
        # col = 'A'
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        df = self.langley_fitres.langley_fitres.sel(channel = col).to_pandas()
        df_corr = self.langley_fitres.correlation_matrix_properties.to_pandas()
        # ['intercept_stderr'].langley_fitres.sel(channel = col).to_pandas()

        df.sort_values('intercept_stderr', inplace=True)
        # df
        df_corr = df_corr.loc[df.index]
        # return df, df_corr  
        
        df.index = range(df.shape[0])
        df_corr.index = range(df.shape[0])
        inter_v0 = np.e**df.intercept
        df['intercept_stderr'] = np.e**(df['intercept'] + df['intercept_stderr']) - inter_v0
        df['intercept'] = inter_v0
        
        

        if not isinstance(best, type(None)):
            df = df.iloc[:best]
            df_corr = df_corr.iloc[:best]

        f, aa = plt.subplots(4, sharex=True, gridspec_kw={'hspace': 0})
        f.set_figheight(f.get_figheight() * 1.5)
        a_int = aa[0]
        a_slope = aa[1]
        a_curv = aa[2]
        a_corr = aa[3]
        
        ### intercept

        df.intercept.plot(ax = a_int, marker = 'o')
        at_int = a_int.twinx()
        df.intercept_stderr.plot(ax = at_int, color = colors[1], marker = 'o')
        at_int.axhline(y = std_error_goal[0], ls = '--', color = 'black')
        

        ### slope
        df.slope.plot(ax = a_slope, marker = 'o')
        at_slope = a_slope.twinx()
        df.slope_stderr.plot(ax = at_slope, color = colors[1], marker = 'o')
        at_slope.axhline(y = std_error_goal[1], ls = '--', color = 'black')

        ### curvature
        df.resid_curvature.plot(ax = a_curv, marker = 'o')
        a_curv.axhline(y = std_error_goal[2], ls = '--', color = 'black')


        ### correlation
        df_corr.plot(ax = a_corr, marker = 'o')        
        a_corr.axhline(y = std_error_goal[3], ls = '--', color = 'black')

        ### Labels 
        a_int.set_ylabel('V0 (mV)')
        at_int.set_ylabel('stderr', color = colors[1])

        a_slope.set_ylabel('slope')
        at_slope.set_ylabel('stderr', color = colors[1])

        a_curv.set_ylabel('curvature')
        # return a_corr,df_corr
        return aa
    
    def plot_data_with_fit(self, top = 0, wls = 'all'):
        which = top

        path2ncfld = pathlib.Path('./langleys/')
        best10 = self
#         which = 0

        best10sel = best10.langley_fitres
        dtt = best10sel.langley_fitres.sel(linreg_param = 'intercept_stderr', channel = 'A').to_pandas().sort_values().index[which]
        date = pd.to_datetime(dtt)
        date = f'{date.year:04d}{date.month:02d}{date.day:02d}'

        # sn = int(best10sel.sp02_serial_no.values[0])
        sn = best10sel.sp02_serial_no.values
        if len(sn.shape) == 0:
            sn = int(sn)
        else:
            sn = int(sn[0])
        
        fn = f'sn{sn}_{date}_am.nc'

        ds = xr.open_dataset(path2ncfld.joinpath(fn))

        f, aa = plot_langley(ds, wls = wls, date = f'{which} - {date}')
        return aa

class CalibrationsOverTime_SP02(object):
    def __init__(self, path2historic = '/export/htelg/projects/sp02/calibration/BRWsp02calhist20062018SEDcorr.dat', list0path2modern = [], additional_result_instance = None):
        "SP02 related"
        self.paths = dict(path2historic = pathlib.Path(path2historic),
                          list0path2modern = [pathlib.Path(p) for p in list0path2modern])
        self._historic = None
        self._results = None
        self._modern  = None
        self._additional = None
        self.additional = additional_result_instance
#         self._additional = additional_result_instance
        
    @property
    def additional(self):
        return self._additional

    @additional.setter
    def additional(self,value):
        if not isinstance(self._additional, type(None)):
            self._remove_additional()
        self._additional = value
        self._add_additional()
        
    def _add_additional(self):
        if isinstance(self._additional, type(None)):
            return
        if not isinstance(self._results, type(None)):
            meana, stda = self._ds2resutltdf(self.additional.langley_fitres)
            adt = pd.to_datetime(self.additional.langley_fitres.datetime.mean().values)
            self._results['mean'].loc[adt] = meana
            self._results['std'].loc[adt] = stda
    
    def _remove_additional(self):
        mean = self.results['mean']
        adt = pd.to_datetime(self.additional.langley_fitres.datetime.mean().values)
        mean.drop([adt], inplace = True)
        return 
    
    @property
    def modern(self):
        if isinstance(self._modern, type(None)):
            self._load_modern()
        return self._modern
    
    def _load_modern(self):
        return [] 
    

    @property
    def results(self):
        if isinstance(self._results, type(None)):
            mean = self.historic_cals.copy()
            std = self._historic_std

            self._results = dict(mean = mean,
                                 std = std)
            
            
            for p2m in self.paths['list0path2modern']:
                langley_fitres = xr.open_dataset(p2m)
                meana, stda = self._ds2resutltdf(langley_fitres)
                adt = pd.to_datetime(langley_fitres.datetime.mean().values)
                self._results['mean'].loc[adt] = meana
                self._results['std'].loc[adt] = stda
            
            self._add_additional()
            
            # hist
            # modern
            #additional
        return self._results
            
                      
    @property
    def historic_cals(self):
        if isinstance(self._historic, type(None)):
            self._read_historic()
        return self._historic
    
    def _read_historic(self):
        path2file = self.paths['path2historic']
        df = pd.read_csv(path2file, delim_whitespace=True, skiprows=2,
                    names = ['year', 'doy', 'ch', 'cal'],
                   )

        df.index = pd.to_datetime(df.year, format='%Y') + pd.to_timedelta(df.doy, 'd')

        channels = pd.Series([int(i) for i in '368 1050 610 778 412 500 675 862'.split()])
        channels.index += 1

        df['wl'] = df.ch.apply(lambda x: channels.loc[x])
        df = df.drop(['year', 'doy', 'ch'], axis=1)

        trans_list = []
        for dt in df.index.unique():
            df_tt = df[df.index == dt]
        #     print(df_tt.shape)
            # df_trans.columns = df_trans.loc['wl']
            df_tt = pd.DataFrame([df_tt.cal.values.transpose()], columns=df_tt.wl, index = [dt])
            df_tt.index.name = 'date'
            trans_list.append(df_tt)

        df_v0_ts = pd.concat(trans_list)
        self._historic = df_v0_ts
        self._historic_std = df_v0_ts.copy()
        self._historic_std[:] = np.nan

    def _ds2resutltdf(self,ds):
        out = {}
        dst = np.e**ds.langley_fitres.sel(linreg_param = 'intercept')
        
        cw = ds.channle_wavelengths
        if cw.dims == ('channel',):
            cw = cw.to_pandas()
        elif cw.dims == ('channel','datetime',):
            cw = cw.transpose('channel','datetime').to_pandas().iloc[:,0]
        else:
            raise ValueError(f'should not be possible ....ds.dims = {ds.dims}')

        mean = dst.mean(dim = 'datetime').to_pandas()
    #     print(mean)
    #     print(cw)
        mean.rename(cw, inplace=True)

        std = dst.std(dim = 'datetime').to_pandas()
        std.rename(cw, inplace=True)
        return mean, std
    
    def plot(self, slides = False, serial_no = None, show_change_of_last = True, show_desired_diff = 0.025):
        """
        Gray span shows the exceptable change (plus minus mean*sow_desired_diff)

        Parameters
        ----------
        slides : TYPE, optional
            DESCRIPTION. The default is False.
        serial_no : TYPE, optional
            DESCRIPTION. The default is None.
        show_change_of_last : TYPE, optional
            DESCRIPTION. The default is True.
        show_desired_diff : TYPE, optional
            DESCRIPTION. The default is 0.025.

        Returns
        -------
        out : TYPE
            DESCRIPTION.

        """
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        out = {}
        df_v0_ts = self.results['mean'].sort_index()
        df_v0_ts_error = self.results['std'].sort_index()
        newcols = df_v0_ts.iloc[-1].dropna().index
        df_v0_ts =df_v0_ts.reindex(newcols, axis=1)
        df_v0_ts_error = df_v0_ts_error.reindex(newcols, axis=1)

        if slides:
            date = pd.to_datetime(datetime.datetime.now()).date()
            date = str(date)
            date = date.replace('-','')
            # path2ppt = f'/mnt/telg/projects/sp02/calibration/mlo_v0_over_time_{date}.ppt'
            path2ppt = f'calibration_{date}.ppt'
            
            if type(slides).__name__ == 'Presentation':
                prs = slides
            else:
                prs = pptx.Presentation()
                if 1:
                    title_slide_layout = prs.slide_layouts[0]
                    slide = prs.slides.add_slide(title_slide_layout)
                    title = slide.shapes.title
                    subtitle = slide.placeholders[1]
        
                    title.text = "v0 over time"
                    now = pd.to_datetime(datetime.datetime.now())
        
                    subtitle.text = f"created on {now.date()} by Hagen"


            ttop = 0
            tleft = 0.1
            lefts = [0.1,5]*2
            tops = [0.2]*2 + [3.5]*2
            width = 4.5
        # average over everything before 18
        ### average over all but last

        oldcal = df_v0_ts.iloc[:-1, :].mean()
        
        f,aa = plt.subplots(2,2, constrained_layout=True)
        aa = [a for at in aa for a in at]
        f.set_size_inches(9, 7)
        for e,ch in enumerate(df_v0_ts):
            out_ch = {}
            if slides:
                if np.mod(e,4) == 0:
                    blank_slide_layout = prs.slide_layouts[6]
                    slide = prs.slides.add_slide(blank_slide_layout)
                    txBox = slide.shapes.add_textbox(tleft, ttop, width, 0.2)
                    tf = txBox.text_frame
                    tf.text = f'sn #{serial_no}'
            
        #     df_v0_ts[ch].plot(ax = a)
            a = aa[e]
            a.plot(df_v0_ts[ch].index, df_v0_ts[ch].values,marker = 'o', markersize = 5, ls = '--')

            a.errorbar(df_v0_ts[ch].index, df_v0_ts[ch].values, df_v0_ts_error[ch].values, capsize = 2, ls = '', ecolor = 'red', zorder = 10)
        #     a.errorbar()

            # this allows to look at differences of current (last) calibration to the previous ones
            if show_change_of_last:
                # average over everything but the last
                a.axhline(y = oldcal.loc[ch], color = colors[1])

                # annotate change from previous average
                anhere = df_v0_ts.index[-1] - pd.to_timedelta(4*365, 'd')
                avg = oldcal.loc[ch]
                last = df_v0_ts[ch].values[-1]
                mid = (avg + last)/2
                diff = (last - avg)/avg *100
    
                a.annotate('',
                        xy=(anhere, avg), xycoords='data',
                        xytext=(anhere, last), textcoords='data',
                        arrowprops=dict(arrowstyle="<->",
                                        connectionstyle="arc3"),
                        ha = 'center',
                        va = 'center'
                        )
                a.text(anhere, mid, f" {diff:0.1f} %")
                
                
                # annotate change from last time
                anhere = df_v0_ts.index[-1] - pd.to_timedelta(1*365, 'd')
                avg = df_v0_ts[ch].values[-2]
                last = df_v0_ts[ch].values[-1]
                mid = (avg + last)/2
                diff = (last - avg)/avg *100
                dist = avg - last
    
                a.annotate('',
                        xy=(anhere, avg), xycoords='data',
                        xytext=(anhere, last), textcoords='data',
                        arrowprops=dict(arrowstyle="<->",
                                        connectionstyle="arc3"),
                        ha = 'center',
                        va = 'center'
                        )
                a.text(anhere, mid, f" {diff:0.1f} %")

            else:
                a.axhline(y = df_v0_ts.mean().loc[ch], color = colors[1])

            # for e,at in enumerate(aa):
            mean = df_v0_ts.mean().iloc[e]
            perc = show_desired_diff
            a.axhspan(mean - (mean * perc), mean + (mean * perc),
                      # color = colors[1], 
                      color = '0.8',
                      # alpha = 0.3
                      )

            a.set_title(f'channel: {ch} nm')
            a.grid()
            a.set_ylabel('V0')
            a.set_xlabel('Date of calibration')
            if slides:
                add_fig2slide(f, slide,left = lefts[np.mod(e,4)], top = tops[np.mod(e,4)], width = width)
        #     break
        if slides:
            prs.save(path2ppt)
            out['slides'] = prs
        out['f'] = f
        out['aa'] = aa
        out['df_v0_ts'] = df_v0_ts
        return out
    


# class SP02RawData(object):
#     def __init__(self, dataset, site, langley_fit_settings = None):
#         self.raw_data = dataset
#         self.site = site
#         self.langley_fit_settings = langley_fit_settings
#         self._sun_position = None
#         self._am = None
#         self._pm = None
#         # self._langleys_am = None
#         # self._langleys_pm = None
#         # self._langley_fitres_am = None
#         # self._langley_fitres_pm = None
        
#     @property
#     def sun_position(self):
#         if isinstance(self._sun_position, type(None)):
#             self._sun_position = self.site.get_sun_position(self.raw_data.datetime)
#         return self._sun_position
    
#     @property
#     def am(self):
#         if isinstance(self._am, type(None)):
#             self._get_langley_from_raw() 
#         return self._am
    
#     @property
#     def pm(self):
#         if isinstance(self._pm, type(None)):
#             self._get_langley_from_raw() 
#         return self._pm
    
#     def tp_get_rdl(self):
#         raw_df = self.raw_data.raw_data.to_pandas()
        
#         # changing to local time
#         raw_df_loc = raw_df.copy()
#         index_local = raw_df.index + pd.to_timedelta(self.site.time_zone[1], 'h')
#         raw_df_loc.index = index_local
#         self.raw_df_loc = raw_df_loc
        
    
#     def _get_langley_from_raw(self):
#         raw_df = self.raw_data.raw_data.to_pandas()
        
#         #### changing to local time
#         raw_df_loc = raw_df.copy()
#         index_local = raw_df.index + pd.to_timedelta(self.site.time_zone[1], 'h')
#         raw_df_loc.index = index_local
#         # self.tp_rdl = raw_df_loc.copy()
        
#         ##### getting the one day
#         sunpos = self.sun_position.copy()
#         start = raw_df_loc.index[0]
#         if sunpos.iloc[0].airmass > 0:
#             start = pd.to_datetime(f'{start.year}{start.month:02d}{start.day:02d}') + pd.to_timedelta(1,'d')
#         end = start + pd.to_timedelta(1, 'd')
#         raw_df_loc = raw_df_loc.truncate(start, end)

#         #### localize and cut day for sunposition
#         sunpos.index = index_local
#         sunpos = sunpos.truncate(start, end)

#         #### remove the night
#         sunpos[sunpos.airmass < 0] = np.nan

#         #### get the minimum airmass befor I start cutting it out
#         noon = sunpos.airmass.idxmin()

#         #### normalize to the sun_earth_distance
#         raw_df_loc = raw_df_loc.multiply(sunpos.sun_earth_distance**2, axis=0)
    
#         # langleys are the natural logarith of the voltage over the AMF ... -> log
#         # to avoid warnings and strange values do some cleaning before log
#         raw_df_loc[raw_df_loc <= 0] = np.nan
# #         self.tp_raw_df = raw_df.copy()
#         raw_df_loc = np.log(raw_df_loc)    
    
#         # keep only what is considered relevant airmasses
#         amf_min = 2.2 
#         amf_max = 4.7
#         sunpos[sunpos.airmass < amf_min] = np.nan
#         sunpos[sunpos.airmass > amf_max] = np.nan

#         sunpos_am = sunpos.copy()
#         sunpos_pm = sunpos.copy()

#         sunpos_am[sunpos.index > noon] = np.nan
#         sunpos_pm[sunpos.index < noon] = np.nan
        

#         langley_am = raw_df_loc.copy()
#         langley_pm = raw_df_loc.copy()

#         self.tp_sp_am = sunpos_am
#         self.tp_sp_pm = sunpos_pm
#         self.tp_df_am = langley_am[~sunpos_am.airmass.isna()].copy()
#         self.tp_df_pm = langley_am[~sunpos_pm.airmass.isna()].copy()

#         langley_am.index = sunpos_am.airmass
#         langley_pm.index = sunpos_pm.airmass

#         self._am = Langley(self,langley_am[~langley_am.index.isna()], langley_fit_settings = self.langley_fit_settings)
#         self._pm = Langley(self,langley_pm[~langley_pm.index.isna()], langley_fit_settings = self.langley_fit_settings)
#         return True

def load_calibration_history():
    l0m = ['/export/htelg/projects/sp02/calibration/2013_14_mlo_cal/calibration_result_1032.nc',
           '/export/htelg/projects/sp02/calibration/2020_summer_mlo_cal/calibration_result_1032.nc',
           # '/mnt/telg/projects/sp02/calibration/2020_mlo_cal/calibration_result_1032.nc',
           '/export/htelg/projects/sp02/calibration/2020_21_mlo_cal/calibration_result_1032.nc',
           '/export/htelg/projects/sp02/calibration/2021_22_mlo_cal/calibration_result_1032.nc',
          ]
    history_1032 = CalibrationsOverTime(list0path2modern=l0m)
    
    l0m = ['/export/htelg/projects/sp02/calibration/2020_summer_mlo_cal/calibration_result_1046.nc',
           # '/mnt/telg/projects/sp02/calibration/2020_mlo_cal/calibration_result_1046.nc',
           '/export/htelg/projects/sp02/calibration/2013_14_mlo_cal/calibration_result_1046.nc',
           '/export/htelg/projects/sp02/calibration/2020_21_mlo_cal/calibration_result_1046.nc',
           '/export/htelg/projects/sp02/calibration/2021_22_mlo_cal/calibration_result_1046.nc']
    history_1046 = CalibrationsOverTime(list0path2modern=l0m)
    
    calibrations = {1032: history_1032,
                    1046: history_1046}
    return calibrations
    
def raw2langleys(site,
                 start_date = '20200205',
                 end_date = '20200207',
                 path2netcdf = '/mnt/telg/data/baseline/mlo/2020/',
                 pattern = '*.nc',
                 pathout_fld = '.',
                 # local_day_in_one_file = False,
                 break_when_error = False,
                 test = False):
    """
    Note
    -----
    This is very similar ot the langley product that is generated by the script ...
    In fact thie script is based on this here file.
    
    Deprecation warning
    --------------------
    The script mentioned above should eventually replace 
    this here file.
    
    Takes the raw data (from the netcdf files) and generates the langleys

    Parameters
    ----------
    site : TYPE
        DESCRIPTION.
    start_date : TYPE, optional
        DESCRIPTION. The default is '20200205'.
    end_date : TYPE, optional
        DESCRIPTION. The default is '20200207'.
    path2netcdf : TYPE, optional
        DESCRIPTION. The default is '/mnt/telg/data/baseline/mlo/2020/'.
    pattern : TYPE, optional
        DESCRIPTION. The default is '*.nc'.
    pathout_fld : TYPE, optional
        DESCRIPTION. The default is '.'.
    break_when_error : TYPE, optional
        DESCRIPTION. The default is False.
    test : TYPE, optional
        DESCRIPTION. The default is False.

    Deprecated
    ----------
    local_day_in_one_file : bool, optional
        Usually files contain a single day in UTC, therefore, one need to 
        combine two files to get the entire (am and pm). Older versions from
        Longenecker took care of this; choose True in this case. The default
        is False.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # local_day_in_one_file = False
    out = {}
    mlo = site
    if not isinstance(path2netcdf, list):
        path2netcdf = [path2netcdf]
        
    path2netcdf = [pathlib.Path(p2n) for p2n in path2netcdf]
    # path2netcdf.mkdir()

    #### generate folder names and work plan
    runname = 'langleys'
    
    pathout_fld = pathlib.Path(pathout_fld)
    pathout_fld = pathout_fld.resolve()
    pathout_fld.mkdir(exist_ok=True)
    
    pathout_fld_run = pathout_fld.joinpath(runname)
    pathout_fld_run.mkdir(exist_ok=True)
    
    # df = pd.DataFrame(list(path2netcdf.glob(pattern)), columns=['path'])
    
    p2nglob = [p2f for p2n in path2netcdf for p2f in list(p2n.glob(pattern))]
    df = pd.DataFrame(p2nglob, columns=['path'])
    
    ### start time index
    df.index = df.apply(lambda row: pd.to_datetime(row.path.name.split('.')[4]), axis=1)
    df.sort_index(inplace=True)

    df_sel = df.truncate(before = start_date,
                         after = end_date
                        )
    
    #### open first file to get serial number for output names and testing
    ds = xr.open_dataset(df_sel.iloc[0,0])
    serial_no = int(ds.serial_no.values)
    # return ds
    ds.close()
    
    #### generate output filespaths and check if exists
    df_sel['path_out'] = df_sel.apply(lambda x: pathout_fld_run.joinpath(f'sn{serial_no}_{x.name.year:04d}{x.name.month:02d}{x.name.day:02d}_am.nc'), axis = 1)
    df_sel['output_file_exists'] = df_sel.path_out.apply(lambda x: x.is_file())    
    total_no = df_sel.shape[0] # -1 since I alsways process 2 at a time
    df_sel = df_sel[~df_sel.output_file_exists]
    work_no = df_sel.shape[0]
    

    print(f'processing {work_no} of a total of {total_no} files: ', end = '')
    # return df_sel, serial_no, ds, pathout_fld_run
    
    reslist_am = []
    # dstl = []
    
    #### go through workplan
    for i in range(len(df_sel)):
        print(i,end = '.')
        
        try:
            row_duo = df_sel.iloc[i:i+2]
        except KeyError:
            assert(False), 'this error was expected, take care of it now. '
        ds = xr.open_mfdataset(row_duo.path, combine='by_coords')
        try:
            serial_no_raw = ds.serial_no.values[0]
        except IndexError:
            serial_no_raw = ds.serial_no.values
        
        raw = sp02.products.raw_nc.SP02RawData(ds, mlo)        
        if serial_no_raw != serial_no:
            raise ValueError('serial no is {serial_no_raw}, but should be {serial_no}')
        try:
            raw.am
        except:
            print(f'problem with langley retieval at i = {i}')
            if break_when_error:
                break
            else:
                work_no -= 1
                continue


        try:
            raw.am.langley_fitres
        except:
            print(f'problem with fit at i = {i}')
            if break_when_error:
                break
            else:
                work_no -= 1
                continue

        reslist_am.append(dict(datetime = df_sel.index[i],
                               langley_fitres = raw.am.langley_fitres.copy(),
                               langleys = raw.am.langleys.copy(),
                               correlation_prop = raw.am.langley_residual_correlation_prop.copy(),
                               langley_fit_resid = raw.am.langley_fit_residual.copy(),
                              )
                         )
        # datestr = df_sel.index[i].date().__str__().replace('-', '')
        # sn = int(np.unique(raw.raw_data.serial_no.values))
        # pathout_langleyres = pathout_fld_run.joinpath(f'sn{sn}_{datestr}_am.nc')
        
        # raw.am['channle_wavelengths'] = raw.raw_data.channle_wavelengths.to_pandas().iloc[1]
        # return raw
        raw.am.save2netcdf(df_sel.iloc[i].path_out)
        ds.close()
        if test:
            break
    #     break
    print('done')
    
    pathout_file = pathout_fld_run.joinpath(f'sn{serial_no}_all_fitres.nc')
    
    if work_no > 0:
        ### integrate everthing in one Dataset
        fitresarray = np.array([res['langley_fitres'].values for res in reslist_am])
        coord_time = [res['datetime'] for res in reslist_am]
        correlation_props = [res['correlation_prop'] for res in reslist_am]    
        coord_linreg = raw.am.langley_fitres.columns.values
        # coord_amf = raw.am.langleys.index
        coord_channels = raw.am.langley_fitres.index.values
        langley_fitres_da = xr.DataArray(fitresarray, coords=[coord_time, coord_channels, coord_linreg], dims=['datetime', 'channel', 'linreg_param'])
        
        
        ds = xr.Dataset()
        ds['langley_fitres'] = langley_fitres_da
        # ds['langleys'] = langleys_da
        # ds['langley_fit_residual'] = langleys_resid_da
        
        cmp = pd.DataFrame(correlation_props, index = coord_time)
        cmp.index.name = 'datetime'
        cmp.columns.name = 'matrix_properties'
        ds['correlation_matrix_properties'] = cmp
        
        channel_wl = raw.raw_data.channle_wavelengths.to_pandas()#.iloc[1]
        print(f'channel_wl: {channel_wl.shape}')
        if len(channel_wl.shape)==2:
            channel_wl = channel_wl.iloc[1]
        ds['channle_wavelengths'] = channel_wl
        ds['sp02_serial_no'] = serial_no
    
    
        
        if pathout_file.exists():
            ds_others = xr.open_dataset(pathout_file) 
            ds_new = xr.concat([ds, ds_others], dim = 'datetime')
            ds_new = ds.sortby('datetime')
            ds_new['channle_wavelengths'] = ds_new['channle_wavelengths']
            ds_new['sp02_serial_no'] = ds['sp02_serial_no']
            ds = ds_new
            ds_others.close()
    
        ds.to_netcdf(pathout_file)
    else:
        if pathout_file.exists():
            ds = xr.open_dataset(pathout_file) 
        else:
            ds = None
        raw = None
        
        
    ds.close()
    out['ds_results'] = ds
    out['last_raw_instance'] = raw
    return Calibration(out)

def add_fig2slide(f, slide, left = 0, top = 0, width = 5.5, add_text_box = False):
    imstream = io.BytesIO()
    # f.patch.set_alpha(0)
    f.tight_layout()
    f.savefig(imstream, format='png', bbox_inche = 'tight')

#     blank_slide_layout = prs.slide_layouts[6]
#     slide = prs.slides.add_slide(blank_slide_layout)

    left =  pptx.util.Inches(left)
    top = pptx.util.Inches(top)
    width = pptx.util.Inches(width)
    pic = slide.shapes.add_picture(imstream, left, top, width = width)


#         slide = prs.slides.add_slide(blank_slide_layout)
    left = width # pptx.util.Inches(0.5)
#         width = pptx
    top = top + pptx.util.Inches(0.3)
    height = pptx.util.Inches(3)
    width = pptx.util.Inches(4)
    if add_text_box:
        txBox = slide.shapes.add_textbox(left, top, width, height)
        tf = txBox.text_frame
        tf.text = '...'
    return

def deprecated_plot_langley_results(result_ds, intensity4ts = 'resid_curvature', add2pptx = False,test = False, tollerance = [2.5, 2.5]):
    """
    slope 	intercept 	slope_stderr 	intercept_stderr 	resid_curvature

    Parameters
    ----------
    result_ds : TYPE
        DESCRIPTION.
    intensity4ts : TYPE, optional
        DESCRIPTION. The default is ''.
    add2pptx : TYPE, optional
        DESCRIPTION. The default is False.
    test : TYPE, optional
        DESCRIPTION. The default is False.
    tollerance : TYPE, optional
        DESCRIPTION. The default is [2.5, 2.5].

    Returns
    -------
    out : TYPE
        DESCRIPTION.

    """
    ds = result_ds
    out = {}
    if add2pptx:
        date = pd.to_datetime(datetime.datetime.now()).date()
        date = str(date)
        date = date.replace('-','')
        path2ppt = f'/mnt/telg/projects/sp02/calibration/langley_calibration_summary_{date}.ppt'
        path2ppt
        prs = pptx.Presentation()
        title_slide_layout = prs.slide_layouts[0]
        slide = prs.slides.add_slide(title_slide_layout)
        title = slide.shapes.title
        subtitle = slide.placeholders[1]

        title.text = "SP02 revival"
        now = pd.to_datetime(datetime.datetime.now())

        subtitle.text = f"created on {now.date()} by Hagen"



    plt.rcParams['hatch.linewidth'] = 4
    for wl in ds.channel.values:
        f,a = plt.subplots()
        res = ds.langley_fitres.sel(channel = wl).to_pandas()

    #     x,y,z = res.slope.copy(), res.stderr.copy(), res.intercept.copy()
        rest = res.copy()
        if 0:
            c = rest.slope.median()
            d = abs(c) * 0.5
            xlim =((c - d),(c + d))
            rest.slope[rest.slope > xlim[1]] = np.nan
            rest.slope[rest.slope < xlim[0]] = np.nan
        else:
            xlim = (None, None)
            
        
        if 1:
            bottom = res.intercept_stderr.min()
            top = tollerance[1] * bottom
            dlim = (top - bottom) * 0.1
            ylim = (bottom - dlim , top + dlim)
    #         print(ylim)
        #     ylim =(0,(c + d))
            rest.intercept_stderr[rest.intercept_stderr > top] = np.nan
            rest.intercept_stderr[rest.intercept_stderr < bottom] = np.nan
        else:
            ylim = (None, None)

        restdna = rest.dropna()

        noall = res.shape[0]
        inside = restdna.shape[0]
    #     inside = np.logical_and(~x.isna(), ~y.isna()).sum()
        outside = noall - inside
        txt = f'{outside} of {noall} data\npoints not within\nplotting limits'
    #     print(txt)
        a.text(0.95,0.95, txt,ha = 'right', va = 'top', transform=a.transAxes,)
    #     mi = 0
    #     ma = 10
    #     y[y > ma] = np.nan
    #     y[y < mi] = np.nan

        if 0:
            pc = plt.hexbin(res.slope, res.stderr, gridsize=100)
            cb = f.colorbar(pc)
        #     pc.set_linewidth(0.2)
            pc.set_edgecolor([0,0,0,0])
            cm = plt.cm.plasma_r
            pc.set_cmap(cm)
            cm.set_under('w')
            pc.set_clim(0.01)
            cb.set_label('# of datapoints on bin')

        # elif 0:
        #     g, = plt.plot(x, y)
        #     g.set_linestyle('')
        #     g.set_marker('.')
        #     g.set_markersize(1)
        #     a.set_xlim(xlim)
        #     a.set_ylim(ylim)
        else:
            cm = plt.cm.plasma_r
            pc = plt.scatter(restdna.slope,restdna.intercept_stderr, c = restdna.intercept, s = 40, cmap = cm, 
    #                     edgecolors='black',
                       )
            cb = f.colorbar(pc)
            cb.set_label('Langley fit - intercept')
            zm = restdna.intercept.median()
            zstd  = restdna.intercept.std()
            cb.ax.axhline(zm, color = 'black')
            hs = cb.ax.axhspan(zm-zstd, zm+zstd)
            hs.set_facecolor([0,0,0,0])
            hs.set_edgecolor([0,0,0,0.5])
    # #         hs.set_edgecolor([1,1,1,1])
            hatch = hs.set_hatch('//')
    #         hs.set_linewidth(10)
            a.set_xlim(xlim)
            a.set_ylim(ylim)
    #         print(ylim)
        a.set_title(f'channel: {wl}nm')
        a.set_xlabel('Langley fit - slope (mV)')
        a.set_ylabel('Langley fit - std of residual')

    # langley over time
        f_ivt,a = plt.subplots()
        out['tp_restdna'] = restdna
        # restdna.intercept.plot(ax = a)
        pc = a.scatter(restdna.index, restdna.intercept, c = restdna[intensity4ts], cmap=plt.cm.gnuplot)
        f_ivt.colorbar(pc)
        # g = a.get_lines()[0]
        # g.set_linestyle('')
        # g.set_marker('o')
        a.set_xlabel('')
        a.set_ylabel('Langley fit - intercept')
        a.set_title('Langley fit intercept as a function of time')

        if add2pptx:

            blank_slide_layout = prs.slide_layouts[6]
            slide = prs.slides.add_slide(blank_slide_layout)
            add_fig2slide(f, slide)
            add_fig2slide(f_ivt, slide, top = 7.5/2)
        if test:
            break
    if add2pptx:
        prs.save(path2ppt)
    return out

def deprecated_get_best_results(ds, top = 10):
    def sortbyparam(ds, param, top):
        df = (ds.langley_fitres.sel(linreg_param = f'{param}_stderr') / ds.langley_fitres.sel(linreg_param = param)).mean(dim = 'channel').to_pandas()
        df = df[df>0]
        df.sort_values(inplace=True)
        out = {}
        df = df.iloc[:top]
        out = ds.sel(datetime =  df.index.values)
        return out
    out = {}
    slopestdrmin = sortbyparam(ds,'slope', top)
    interceptstdrmin = sortbyparam(ds,'intercept',top)
    curvature_date = ds.langley_fitres.sel(linreg_param = 'resid_curvature').mean(dim = 'channel').to_pandas().abs().sort_values().index[:top].values
    curvature = ds.sel(datetime = curvature_date)
    df_cmp = ds.correlation_matrix_properties.to_pandas()
    determimin_dates = df_cmp.sort_values('determinant', ascending=False).index[:top].values
    determimin = ds.sel(datetime = determimin_dates)
#     uberlap = ds.datetime.values
#     for other in [interceptstdrmin, 
#                  curvature, 
#     #              determimin,
#                  ]:
#         uberlap = np.intersect1d(uberlap, other)
#     out['crossection'] = uberlap
    
    out['slope_stderr'] = slopestdrmin
    out['intercept_stderr'] = interceptstdrmin
    out['curvature'] = curvature
    out['corr_matrix_determinant_max'] = determimin
    return out


def future_add_fig2slide(f, slide, left = 0, top = 0, width = 5.5):
    imstream = io.BytesIO()
    f.patch.set_alpha(0)
    f.tight_layout()
    f.savefig(imstream, format='png', bbox_inche = 'tight')

#     blank_slide_layout = prs.slide_layouts[6]
#     slide = prs.slides.add_slide(blank_slide_layout)

    left =  pptx.util.Inches(left)
    top = pptx.util.Inches(top)
    width = pptx.util.Inches(width)
    pic = slide.shapes.add_picture(imstream, left, top, width = width)


#         slide = prs.slides.add_slide(blank_slide_layout)
#     left = width # pptx.util.Inches(0.5)
# #         width = pptx
#     top = top + pptx.util.Inches(0.3)
#     height = pptx.util.Inches(3)
#     width = pptx.util.Inches(4)
#     txBox = slide.shapes.add_textbox(left, top, width, height)
#     tf = txBox.text_frame
#     tf.text = '...'
    return
