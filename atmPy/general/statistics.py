import pandas as _pd
import matplotlib.pylab as _plt
import numpy as _np
# from matplotlib.dates import MonthLocator  as _MonthLocator
# from matplotlib.dates import DateFormatter as _DateFormatter
from matplotlib.ticker import FuncFormatter as _FuncFormatter
from matplotlib.ticker import MultipleLocator as _MultipleLocator
import matplotlib.gridspec as mplgrid
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerTuple
import matplotlib.colors as _mcolors
colors = _plt.rcParams['axes.prop_cycle'].by_key()['color']
#import plt_tools as _plt_tools
from atmPy.tools import array_tools as _array_tools
import atmPy.tools.plt_tool_kit.colors as _atmcols
# from pygam import LinearGAM, pygam_s, pygam_l

import datetime
import xarray as _xr 
import matplotlib.dates as mdates
import  matplotlib.lines as _mpllines
import matplotlib.dates as _mpldates

# from statsmodels.graphics.tsaplots import plot_acf as statsmod_plot_acf
from atmPy.opt_imports import statsmodels as sm
from atmPy.opt_imports import pygam

import scipy.stats as scistats


class Statistics(object):
    def __init__(self, parent_ts):
        self._parent_ts = parent_ts
        self.seasonality = Seasonality(self)
        self.diurnality = Diurnality(self)
        self.gamcl = GamClimatology(self)

    def define_custom(self,reference = 'index', frequency = 'M', bins = None):
        self.custom = Climatology(self, reference = reference, frequency = frequency, bins = bins)

class GamClimatology(object):
    def __init__(self, parent_stats = None):
        '''
        This class hosds a few tools for a GAM analysis. It is design to take
        the Statistics class as an argument to connect it to its parents, in
        particular, the Timeseries class. It can also invoced with an argument. 
        You can than overwrite properties, e.g. prediction, and can use the plot
        function

        Parameters
        ----------
        parent_stats : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        '''
        self.distribution = 'normal'# 'normal' ['binomial' 'poisson' 'gamma' 'inv_gauss']
        self.link = 'identity' # 'identity' ('logit' 'inverse' 'log' 'inverse-squared')
        self.trend_model = 'spline' # [spline, linear'] this is the overall trend model
        self.interactions = False
        # self.interactions_lam = None
        # self.interactions_nsplines = None
        self.interactions_only = False
        self._splines_per_year = 15
        self._parent_stats = parent_stats
        self._data = None
        self._seasonality_lam = None
        self._seasonality_nsplines = None
        self._trend_lam = None
        self._trend_nsplines = None
        self._prediction_grid_size = None
        self._prediction_confidence = None
        self._autocorrelation = None
        self._rerun()
        
    def _rerun(self):
        self._fit_res = None
        self._prediction = None
        self._gam = None
    
    @property
    def prediction_confidence(self):
        if isinstance(self._prediction_confidence, type(None)):
            self._prediction_confidence = 0.95
            
        return self._prediction_confidence
    
    @prediction_confidence.setter
    def prediction_confidence(self,value):
        assert(0<value<1)
        self._prediction_confidence = value
        self._prediction = None
        return
    
    @property
    def prediction_grid_size(self):
        if isinstance(self._prediction_grid_size, type(None)):
            self._prediction_grid_size = 100
        return self._prediction_grid_size
    
    @prediction_grid_size.setter
    def prediction_grid_size(self, value):
        self._prediction_grid_size = value
        self._prediction = None
        return
        
            
    @property
    def splines_per_year(self):
        return self._splines_per_year
    
    @splines_per_year.setter 
    def splines_per_year(self, value):
        self._splines_per_year =  value
        self._rerun()
        
        self._seasonality_nsplines = None
        self._trend_nsplines = None

    
    
    @property
    def trend_nsplines(self):
        if isinstance(self._trend_nsplines, type(None)):
            # self._trend_nsplines  = 1e1 #(self.data.index.max() - self.data.index.min()) / _pd.to_timedelta(1, 'd')
            n_splines = (self.data.index.max() - self.data.index.min())/ _pd.to_timedelta(1, 'd')/356 # no_of years
            n_splines*= self.splines_per_year # no of month
            # n_splines*= 4 # random
            self._trend_nsplines = int(n_splines)
        return self._trend_nsplines
    
    @trend_nsplines.setter
    def trend_nsplines(self, value):
        self._trend_nsplines = int(value)
        self._rerun()
        
        
    @property
    def trend_lambda(self):
        if isinstance(self._trend_lam, type(None)):
            self._trend_lam  = (self.data.index.max() - self.data.index.min()) / _pd.to_timedelta(1, 'd')
        return self._trend_lam
    
    @trend_lambda.setter
    def trend_lambda(self, value):
        self._trend_lam = value
        self._rerun()
        
        
    @property
    def seasonality_nsplines(self):
        if isinstance(self._seasonality_nsplines, type(None)):
            self._seasonality_nsplines  = self.splines_per_year
        return self._seasonality_nsplines
    
    @seasonality_nsplines.setter
    def seasonality_nsplines(self, value):
        self._seasonality_nsplines = value
        self._rerun()


    @property
    def seasonality_lambda(self):
        if isinstance(self._seasonality_lam, type(None)):
            self._seasonality_lam  = 1e1 #(self.data.index.max() - self.data.index.min()) / _pd.to_timedelta(1, 'd')
        return self._seasonality_lam
    
    @seasonality_lambda.setter
    def seasonality_lambda(self, value):
        self._seasonality_lam = value
        self._rerun()
    
    @property
    def data(self, data_column = None):#, linear = False, resolution = 1e5):
        
        """
        resolution: 1e5 ~ seasonal
                    1e7 ~ annual
        """
        if isinstance(self._data, type(None)):
            #### make_X_y_Xdf
            # dsn = xr.Dataset()
            # dsn['aod'] = da
            data = self._parent_stats._parent_ts.data
            
            Xdf = _pd.DataFrame()
            
            row_s = data.iloc[0,:]
            start_date = row_s.name
            self._start_date = start_date
            Xdf['dsincestart'] = data.apply(lambda row: (row.name - start_date )/ datetime.timedelta(days = 1), axis = 1)
            Xdf['doy'] = data.apply(lambda row: (row.name - _pd.to_datetime(row.name.year, format= '%Y')) / datetime.timedelta(days = 1), axis = 1)
            # if sun:
            #     Xdf['sunelev'] = data['sunelev']
            #     Xdf['sunaz'] = data['sunaz']
            # X = Xdf.values
            if isinstance(data_column, type(None)):
                y = data.iloc[:,0]
            else:
                y  = data.loc[:, data_column]
                
            # y.name = f'y_{y.name}'
            y.name = 'y'
            Xdf.columns = [f'x_{c}' for c in Xdf.columns]
            data_pretty = _pd.concat([y,Xdf], axis = 1)
            # data_pretty = 
            self._data = data_pretty#(data, Xdf)
        return self._data
    
    def optimize_lambda(self, elim = (-3,3), size = 100, lambda_grid = None):
        if isinstance(lambda_grid, type(None)):
            lams = create_lambda_set(no_col = self.data.shape[1] - 1, elim = elim, size = size)
        else: 
            lams = lambda_grid
            
        self.lambda_grid = lams
        
        X = self.data.iloc[:,1:].values
        y = self.data.iloc[:,0].values
        self.gam_inst.gridsearch(X, y, lam=lams)
        
        # update or reset some of the parameters!
        self._prediction = None 
        self._trend_lam = self.fit_res.terms[0].info['lam'][0]
        self._seasonality_lam = self.fit_res.terms[1].info['lam'][0]
        return None #self._fit_res # this is still the same instance ..
        
    @property
    def gam_inst(self):
        if isinstance(self._gam, type(None)):
            if self.trend_model == 'linear':
                year = pygam.l(0)#, lam = resolution, n_splines=int(n_splines))
            elif self.trend_model == 'spline':
                year = pygam.s(0, lam = self.trend_lambda, n_splines=self.trend_nsplines)#int(n_splines))
            else:
                raise ValueError(f'trend_model has to be either "linear" or "spline". It is {self.trend_model} though.')
                
            doy = pygam.s(1, basis = 'cp', lam = self.seasonality_lambda, n_splines= self.seasonality_nsplines)
            self.tp_term_doy = doy
            models = year + doy

            if self.interactions:
                # Alternatively you could define the terms with s(1, ....). 
                # I tried it with l too, but I could not find out how to add an intersect to that linear term
                # Also tried to just use the doy and year from above (also newly defined so I am not using the very same instances
                # but the problem is that I am using to many spline in those ... that just never finishes, might be worth trying again?
                
                te_term = pygam.te(0,1, 
                                   lam = self.interactions_lam,#(self.trend_lambda, self.seasonality_lambda), 
                                   n_splines = self.interactions_nsplines, #(self.trend_nsplines, self.seasonality_nsplines)
                                   basis = ['ps', 'cp'], # first is for trend, second for the cyclic day of year
                                  )

                
                if self.interactions_only:
                    models = te_term
                else:
                    models += te_term
            
            if self.distribution == 'normal':
                self._gam = pygam.LinearGAM(models)
            else:
                self._gam = pygam.GAM(models , distribution = self.distribution, link = self.link)
        return self._gam
    
    @property
    def fit_res(self, data_column = None, 
                # linear = False, 
                resolution = 1e5):
        
        """
        resolution: 1e5 ~ seasonal
                    1e7 ~ annual
        """
        if isinstance(self._fit_res, type(None)):
            X = self.data.iloc[:,1:].values
            y = self.data.iloc[:,0].values
            
            # if linear:
            #     year = pygam.l(0)#, lam = resolution, n_splines=int(n_splines))
            # else:
            #     year = pygam.s(0, lam = self.trend_lambda, n_splines=self.trend_nsplines)#int(n_splines))
            
            # # doy = pygam.s(1, basis = 'cp', lam = 1e1, n_splines= 12 * 4)
            # doy = pygam.s(1, basis = 'cp', lam = self.seasonality_lambda, n_splines= self.seasonality_nsplines)
            # # sune = s(2, 
            # #          # n_splines= 12 * 4,
            # #         )
            # # suna = s(3)
            # # hod = s(2, basis = 'cp', lam = 1e-3,# by =4, edge_knots = 2)
            
            # gam = pygam.LinearGAM(year + doy )
            # self.tp_gam = gam
            
            self._fit_res = self.gam_inst.fit(X, y) #this is just a nother handle to the very same gam instance from above!! If the lambda optimization is run it still connects to the same gam instance!!
        return self._fit_res
    
    def get_quantile(self, q = None, std = None):
        if isinstance(q, type(None)):
            if isinstance(std, type(None)):
                raise ValueError('Dude, at last one needs to be given!')
            else:
                q = scistats.norm.cdf(std)
        XX = self.data.iloc[:,1:].values
        pred = self.gam_inst.prediction_intervals(XX, quantiles=q)
        dindex = self._start_date + _pd.to_timedelta(XX[:,0], 'd')
        da = _xr.DataArray(pred[:,0], coords = {'datetime_full': dindex})
        return da
    
    @property
    def prediction(self):
        if isinstance(self._prediction, type(None)):
            ds = _xr.Dataset()
            Xdf = self.data.iloc[:,1:]
            Xdf.columns = [c.replace('x_', '') for c in Xdf.columns]
            gam = self.fit_res
            
            # if self.interactions_only:
            #     i = 0
            #     XX = gam.generate_X_grid(term=i, meshgrid=True)#n = self.prediction_grid_size)
            #     pdep = gam.partial_dependence(term=i, X=XX, meshgrid = True)#width=self.prediction_confidence)
            #     self.tp_pdep = pdep
            #     # self.tp_confi = confi
            #     self.tp_ds = ds
            #     self.tp_XX = XX
            #     return
            if not self.interactions_only:
                for i,col in enumerate(Xdf.columns):
                    term = gam.terms[i]
                    XX = gam.generate_X_grid(term=i, n = self.prediction_grid_size)
                    pdep, confi = gam.partial_dependence(term=i, X=XX, width=self.prediction_confidence)
                    self.tp_confi = confi
                    self.tp_pdep = pdep
                    dsincestart_feat = XX[:, term.feature]
                    colname = Xdf.columns[i]
                    if colname == 'dsincestart':
                        startd = Xdf.index[0]
                        index = startd + (dsincestart_feat * _pd.to_timedelta(1, 'day'))
                        colname = 'datetime'
                    else:
                        index = dsincestart_feat
                    # print(pdep.shape)
                    self.tp_colname = colname
                    self.tp_index = index
                    ds[f'partial_{colname}'] =  _xr.DataArray(pdep, coords={colname:index})
                    ds[f'partial_{colname}_confidence'] =  _xr.DataArray(confi, coords={colname:index,
                                                                                        'boundary': ['top', 'bottom']})
                    
            if self.interactions:
                if self.interactions_only:
                    i = 0
                else:
                    i = 2
                # XX_trend = gam.generate_X_grid(term=0, n = self.prediction_grid_size)
                # XX_season = gam.generate_X_grid(term=1, n = self.prediction_grid_size)
                # self.tp_XX_trend = XX_trend
                # self.tp_XX_season = XX_season
                XX = gam.generate_X_grid(term=i,n = self.prediction_grid_size, meshgrid=True)#n = self.prediction_grid_size)
                pred = gam.partial_dependence(term=i, X=XX, meshgrid = True)#width=self.prediction_confidence)
                self.tp_pdep = pred
                # self.tp_confi = confi
                self.tp_ds = ds
                self.tp_XX = XX
                # ds['interaction'] = xr.DataArray(self.tp_pdep, #dims = ('bla', 'blub'), 
                #                                  coords = {'dsincestart': self.tp_XX[0][:,0],
                #                                            'doy' : self.tp_XX[1][0]})

                dt = self.data.index[0] + (XX[0][:,0] * _pd.to_timedelta(1, 'day'))
                ds['interactions'] = _xr.DataArray(pred, #dims = ('bla', 'blub'), 
                                                   coords = {'datetime': dt,
                                                             'doy' : XX[1][0]})
            
            XX = self.data.iloc[:,1:].values
            pred = gam.predict(XX)  
            dindex = self._start_date + _pd.to_timedelta(XX[:,0], 'd')
            ds['prediction'] = _xr.DataArray(pred, coords = {'datetime_full': dindex})
            ds.attrs['intercept'] = gam.coef_[-1]
            if hasattr(gam, 'prediction_intervals'):
                pred = gam.prediction_intervals(XX, width=self.prediction_confidence)
                ds['prediction_confidence'] = _xr.DataArray(pred, coords = {'datetime_full': dindex,
                                                                            'boundary': ['bottom', 'top']})

            #### add residual
            res = self.data.iloc[:,0] - ds.prediction
            res.index.name = 'datetime_full'
            ds['residual'] = res
            
            self._prediction = ds
        return self._prediction
    
    @prediction.setter 
    def prediction(self, value):
        self._prediction = value

    @property
    def autocorrelation(self):
        if isinstance(self._autocorrelation, type(None)):
            self.autocorrelation_create_custom()
        return self._autocorrelation
        
    def autocorrelation_create_custom(self, lagmax = None, z_score = 1.96):
        """
        Two confidence intervals are given:
            Simple: based on number of datapoints in each correlation. Note, numbers are getting less with increasing lag.
            Bartlett’s formula: Brockwell and Davis, 1987. Time Series Theory and Methods. This takes varibiltiy of acf into account.
        """
        autocorr = [1]
        residuals = self.prediction.residual 
        acfconf = [_np.sqrt(1/ residuals.shape[0]) * z_score]
        lags = [0,]
        ci_bartlett = [_np.sqrt(1/ residuals.shape[0]) * z_score]
        # ci_lower_bartlett = []
        
        lag = 0
        # lagmax = None
        lagmax_above_ci = None
        while 1:
            lag += 1
            lags.append(lag)
            dt = _np.array([residuals[:-lag], residuals[lag:]])
            ac = _np.corrcoef(dt)[0, 1]
            autocorr.append(ac)
            
            ci = _np.sqrt(1/ (dt.shape[1])) * z_score
            acfconf.append(ci)
        
        
            # for k in range(len(acf_values)):
            n = residuals.shape[0]
            # n = dt.shape[1]
            bartlett_se = _np.sqrt((1 + 2 * _np.sum(_np.array(autocorr)[1:]**2)) / n)
            ci_b = 1.96 * bartlett_se
            ci_bartlett.append(ci_b)
        
            # max lag test parameter
            # if len(autocorr) >= 3:
            maxlagparam =  max(autocorr[-3:])
            # else:
                # maxlagparam = ac
            # print(acfconf[-3:])
            # print(maxlagparam)
            if isinstance(lagmax_above_ci, type(None)):
                if ci_b > maxlagparam:
                    lagmax_above_ci = lag - 1
                    if isinstance(lagmax, type(None)):
                        lagmax = int(lag * 1.5)
            if not isinstance(lagmax, type(None)):
              if lag >= lagmax:
                  break
                  
            # elif ci_b > maxlagparam:
            #     lagmax_above_ci = lag - 1
            #     lagmax = int(lag * 1.5)
    
        acfdf = _pd.DataFrame({'acf': autocorr, 'ci': acfconf, 'ci_bartlett': ci_bartlett}, index = lags)
        acfdf.index.name = 'lag'
        acfdf.columns.name = 'acf_params'
        if 'autocorrelation' in self.prediction.variables:
            self.prediction = self.prediction.drop_vars(['autocorrelation', 'lag'])
            
        self.prediction['autocorrelation'] = acfdf
        self.prediction.autocorrelation.attrs['lagmax'] = lagmax
        self.prediction.autocorrelation.attrs['lagmax_above_ci'] = lagmax_above_ci
        self._autocorrelation = self.prediction.autocorrelation
        return self.prediction.autocorrelation

    @property
    def standard_error_adjustment_autocorr(self):
        self.autocorrelation #make sure it is generated
        autocorrs = self.prediction.autocorrelation.sel(acf_params = 'acf')[:self.prediction.autocorrelation.lagmax_above_ci]
        rho_sum = autocorrs.sum()
        n = self.gam_inst.statistics_['n_samples']
        n_eff = n / (1 + 2 * rho_sum)
        adjustment = float(_np.sqrt(n/n_eff))
        self.prediction.autocorrelation.attrs['sandard_error_adjustment_autocorr'] = adjustment
        return adjustment

    @property
    def standard_error_adjustment_newey_west(self):
        """I think this only works for linear term, because of the self.gam_inst.statistics_['se'][0] below"""
        self.autocorrelation #make sure it is generated
        data = _pd.DataFrame({'residuals': self.prediction.residual.values, 'x':self.data['x_dsincestart'].values})
        data['const'] = 1
        
        model = sm.api.OLS(data['residuals'], data[['const', 'x']])
        maxlags = self.prediction.autocorrelation.lagmax_above_ci
        result = model.fit(cov_type='HAC', cov_kwds={'maxlags': maxlags})
    
        sandard_error_adjustment_newey_west = result.bse.x/self.gam_inst.statistics_['se'][0]
        self.prediction.autocorrelation.attrs['sandard_error_adjustment_newey_west'] = sandard_error_adjustment_newey_west
        return sandard_error_adjustment_newey_west

    @property
    def confidence_interval_adjusted_newey_wet(self):
        bla   = (self.prediction.partial_datetime_confidence - self.prediction.partial_datetime)
        bla  *=  self.standard_error_adjustment_newey_west
        bla  +=  self.prediction.partial_datetime
        self.prediction['partial_datetime_confidence_adjusted_nw'] = bla
        return bla
        
    def plot_seasonality(self, ax = None,
                         offset = 0,
                         xticklablesmonth = True,
                         show_confidence = True,
                         show_observations = True,
                         orientation='horizontal',
                         transform = None,
                         obs_gridsize = 50,
                         # vcenter = True,
                         **plot_kwargs):
        """"
        Parameters
        ===========
        offset: float or str ['center', 'intercept']
        transform: str ['log2lin', 'percent']
            """
        
        if isinstance(ax, type(None)):
            f, a = _plt.subplots()
        else:
            a = ax
            f = a.get_figure()
        
        # if 'label' not in plot_kwargs:
        #     plot_kwargs['label'] = 'seasonal'
        if orientation == 'vertical':
            y = 'doy'
        else:
            y = None


        # assert(not (vcenter and offset)), 'if vcenter is True, offset is set automatically' 
        if offset == 'center':
            offset = - float(self.prediction.partial_doy.mean())
        elif offset == 'intercept':
            offset = self.prediction.intercept
        
        da = self.prediction.partial_doy + offset
        if transform == 'percent':
            da = (10**da - 1) * 100
        if transform == 'log2lin':
            da = 10**da
            
        ##################################
        #### plot prediction
        ##############
        da.plot(ax = a, y = y, zorder = 10, **plot_kwargs)   
        g = a.get_lines()[-1]
        col_major = g.get_color()
        # if transform == 'percent':
        #     da = (10**self.prediction.partial_doy - 1) * 100
        # else: 
        #     da = self.prediction.partial_doy
        # (da + offset).plot(ax = a, y = y,**plot_kwargs)   


        ###################################
        #### conficence intervals
        #############
        if show_confidence:
            bot,up = self.prediction.partial_doy_confidence.values.transpose()
            bot = bot + offset
            up = up + offset
            if transform == 'percent': 
                bot = (10**bot - 1) * 100
                up = (10**up - 1) * 100
            if orientation == 'vertical':
                fill = a.fill_betweenx(self.prediction.doy,bot,up, 
                                       # color = [0,0,0,0.3],
                                       color = col_major, alpha = 0.4,
                                       zorder = 5 )
            else: 
                fill = a.fill_between(self.prediction.doy,bot,up, 
                                      color = col_major, 
                                      alpha = 0.4,
                                      zorder = 5 )
        else:
            fill = None

        ######################################
        #### observations
        ################
        
        hbin = None
        if show_observations:
            pdts = self.prediction.partial_datetime.to_pandas()
            data = self.data.copy()
            data['partial_dt'] = pdts
            data = data.interpolate()
            data['y_tc'] = data.y - data.partial_dt
            cm = _plt.cm.Oranges_r
            cm.set_under([0,0,0,0])
            y = data.y_tc + offset - self.prediction.intercept
            hbin = a.hexbin(data.x_doy, y, linewidths=0.2, gridsize = obs_gridsize, vmin = 0.1, cmap = cm, zorder = 1)
            hbin.set_clim(vmax = hbin.get_clim()[1] * 0.5)
            pcolors = _plt.cm.Oranges_r(_np.linspace(0, 1, 100))
            patch = [mpatches.Patch(facecolor=c,) for c in pcolors]
            handles=[patch]
            labels=['data density',]
            a.legend(handles, labels, ncol=len(labels), fontsize='small',
                       handler_map = {list: HandlerTuple(None)})
        #######################################
        #### Settings
        #############
        
        if xticklablesmonth:
            if orientation == 'vertical':
                axis = a.yaxis
            else: 
                axis = a.xaxis
            axis.set_major_locator(mdates.MonthLocator())  # Set the major ticks to be at the beginning of each month
            axis.set_minor_locator(mdates.WeekdayLocator())  # Set the minor ticks to be at the beginning of each week
            axis.set_major_formatter(mdates.DateFormatter('%b'))  # Format the major ticks with abbreviated month names
            axis.set_minor_locator(_plt.NullLocator())
            axis.set_label_text('')
        # a.set_xlabel('')
        # a.set_xlabel('Day of year')
        return f,a,fill, hbin
    
    def plot_prediction(self, ax = None, show_confidence=True, show_original_data=True, **plot_kwargs):
        if isinstance(ax, type(None)):
            f, a = _plt.subplots()
        else:
            a = ax
            f = a.get_figure()

    
        if 'label' not in plot_kwargs:
            plot_kwargs['label'] = 'prediction'
        

        self.prediction.prediction.plot(ax = a, zorder = 3, label = 'prediction')
        g = a.get_lines()[-1]
        col_major = g.get_color()
        if show_confidence:
            a.fill_between(self.prediction.datetime_full, 
                           self.prediction.prediction_confidence.sel(boundary='top'), 
                           self.prediction.prediction_confidence.sel(boundary='bottom'), 
                           alpha = 0.3, 
                           zorder = 2, 
                           # color = '0.5', 
                           color = col_major,
                           label = 'confidence')
        if show_original_data:
            data = self._parent_stats._parent_ts.data
            a.plot(data.index, data.iloc[:,0], ls = '', marker = '.', markersize = 1, zorder = 1, label = 'observation', color = '0.4')
        return f,a

    def plot_interactions(self, ax = None, xticklablesmonth = True, transform = None, cbkwargs = {}, offset = 'intercept'): 
        """
        Parameters
        ===========
        cbkwargs: dict or bool
            if False, no colorbar will be shown
        """
        if isinstance(ax, type(None)):
            f,a = _plt.subplots()
        else:
            a = ax 
            f = a.get_figure()

        if offset == 'intercept':
            offset = self.prediction.intercept
        
        da = self.prediction.interactions + offset
        
        if transform == 'percent':
            da = (10**da -1) * 100
        else:
            pass
        
            
        pc = da.plot(x = 'datetime', ax = a, add_colorbar = False)
        pc.set_cmap(_plt.cm.Spectral_r)

        cbar = None
        if cbkwargs:
            cbar = f.colorbar(pc, **cbkwargs)
                              # ax = axcb, 
                              # # location = 'left',
                              # anchor = (0,1.1), 
                              # # pad = 1.5,
                              # aspect = 10,
                              # shrink = 1.8) 
        
        
        if xticklablesmonth:
            axis = a.yaxis
            axis.set_major_locator(mdates.MonthLocator())  # Set the major ticks to be at the beginning of each month
            axis.set_minor_locator(mdates.WeekdayLocator())  # Set the minor ticks to be at the beginning of each week
            axis.set_major_formatter(mdates.DateFormatter('%b'))  # Format the major ticks with abbreviated month names
            axis.set_minor_locator(_plt.NullLocator())
            axis.set_label_text('')
        return f,a, pc, cbar
        
    def plot_overview(self, axis = None, show_confidence = True, show_original_data = True, 
                      shade_seasons = True, transform = None, offset_season = 0,
                      season_obs_gridsize = (80, 25)):
        """
        Parameters
        ===========
        offset_season: see plot_seasonality
        """
        
        out = []
        if isinstance(axis, type(None)):
            if self.interactions:
                f,aa = _plt.subplots(3,2, #sharex=True, sharey=True,
                    height_ratios=[1,1,3], 
                    width_ratios=[3,1], 
                    gridspec_kw={'hspace':0, 'wspace': 0})
                aaf = aa.flatten()
                
                for a in (aa[0,1], aa[1,1]):
                    a.spines['right'].set_visible(False)
                    # a.spines['left'].set_visible(False)
                    a.spines['top'].set_visible(False)
                    a.spines['bottom'].set_visible(False)
                    a.xaxis.set_ticks([])
                    a.yaxis.set_ticks([])
                # a = aa[1,1]
                # a.spines['right'].set_visible(False)
                # a.spines['top'].set_visible(False)
            else:
                aa = []
                f = _plt.figure()
                f.set_figheight(f.get_figheight() * 1.5)
                aa.append(f.add_subplot(3,1,1))
                aa.append(f.add_subplot(3,1,2, 
                                        # sharex = aa[0]
                                       ))
                aa.append(f.add_subplot(3,1,3))
                
            if shade_seasons:
                shade_seasons = {'color': '#02401A', 'alpha': 0.5}
        else:
            aa = axis
            f = aa[0].get_figure()
            shade_seasons = False

        out.append(f)
        out.append(aa)
            
        if self.interactions:
            a_pred = aa[0,0]
            a_trend = aa[1,0]
            a_season = aa[2,1]
            a_inter = aa[2,0]
            orientation='vertical'

        else: 
            a_pred = aa[0]
            a_trend = aa[1]
            a_season = aa[2]
            orientation='horizontal'
        
        xshift = 0
        
        self.plot_prediction(ax = a_pred, show_confidence=show_confidence, show_original_data=show_original_data)
        self.plot_trend(ax=a_trend, shade_seasons=shade_seasons, show_confidence=show_confidence, transform = transform)


        self.plot_seasonality(ax=a_season, 
                              show_confidence = show_confidence, 
                              orientation=orientation, transform = transform, 
                              offset = offset_season,
                              obs_gridsize = season_obs_gridsize)
        
        if self.interactions:
            # a = aa[2,0]
            a = a_inter
            cbkwargs = dict(ax = aa[0,1], 
                               anchor = (0,1.1), 
                               aspect = 10,
                               shrink = 1.8
                    ) 
            f,a,pc,cb = self.plot_interactions(ax = a, cbkwargs = cbkwargs, transform = transform)
            out.append(pc)
            out.append(cb)
            
        if self.interactions:
            for a in (aa[0,0], aa[1,0]):
                a.set_xlim(aa[2,0].get_xlim())
            
            a = aa[2,1]
            a.set_ylim(aa[2,0].get_ylim())
            a.yaxis.tick_right()
            
            cb.ax.yaxis.tick_left()

            a_inter.set_xlabel('')
            
        else:
            for e,a in enumerate(aa):
                if e == 0:
                    continue
                poslast = aa[e-1].get_position().bounds
                pos = list(a.get_position().bounds)
                if e == 1:
                    xshift += poslast[1] - (pos[1] + pos[3])
                pos[1] = pos[1] + xshift
                a.set_position(pos)

            # make sure the trend and full data have the same xrange
            a_trend.set_xlim(a_pred.get_xlim())
            a_pred.set_xticklabels([])
            
        a_trend.legend(fontsize = 'small')
        
        return out
    
    def plot_trend(self, ax = None, 
                   shade_seasons = False, 
                   show_obs = True, 
                   show_confidence = True,
                   show_confidence_adjusted = True,
                   offset= 0,
                   transform = None,
                   **plot_kwargs):
        """
        

        Parameters
        ----------
        ax : TYPE, optional
            DESCRIPTION. The default is None.
        shade_seasons : bool or dict, optional
            If to shade the seasons. If dict is provided they will be used as axvspan kwargs, e.g. color or alpha. The default is False.
        show_confidence : TYPE, optional
            DESCRIPTION. The default is True.
        offset : TYPE, optional
            DESCRIPTION. The default is 0.
        **plot_kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        f : TYPE
            DESCRIPTION.
        a : TYPE
            DESCRIPTION.

        """
        
        if isinstance(shade_seasons, dict):
            shade_kwargs = shade_seasons
            shade_seasons = True
        else:
            shade_kwargs = {}
            
        if not 'color' in shade_kwargs:
            shade_kwargs['color'] = '0.5'
         
        if 'label' not in plot_kwargs:
            plot_kwargs['label'] = 'seasonal'
            
        if 'zorder' not in plot_kwargs:
            plot_kwargs['zorder'] = 5
        
        if isinstance(ax, type(None)):
            f, a = _plt.subplots()
            
        else:
            a = ax
            f = a.get_figure()

        if transform == 'percent':
            da = (10**self.prediction.partial_datetime -1) * 100
        else:
            da = self.prediction.partial_datetime 

          
        
        ###################
        #### Plot trend
        (da + offset).plot(ax = a, **plot_kwargs)   
        g = a.get_lines()[-1]
        col_major = g.get_color()
        
        ##################
        #### Seasonal corrected data
        if show_obs:
            residual = self.prediction.residual
            dat = (da + offset)
            residual_plustrend = dat.interp(datetime = self.prediction.datetime_full) + residual
            if 1:
                g, = residual_plustrend.plot(ax = a, ls = '', marker = '.', markersize = 0.8, zorder = 1)
                g.set_color('0.4')
            else:
                scale = 1
                qm = a.hexbin(mdates.date2num(residual_plustrend.datetime), residual_plustrend, 
                               # gridsize = (scale*50,scale*100), 
                              
                               linewidth = 0.2)
                cm = plt.cm.inferno_r
                cm = plt.cm.gnuplot2
                cm.set_under([1,1,1,0])
                qm.set_cmap(cm)
                qm.set_clim(vmin = 0.1)
        
        ###################
        ### confidence inervals
        if show_confidence:
            if transform == 'percent':
                top = (10**self.prediction.partial_datetime_confidence.sel(boundary='top') -1) * 100
                bottom = (10**self.prediction.partial_datetime_confidence.sel(boundary='bottom') -1) * 100
            else:                
                top = self.prediction.partial_datetime_confidence.sel(boundary='top') + offset
                bottom = self.prediction.partial_datetime_confidence.sel(boundary='bottom') + offset                
            
            a.fill_between(self.prediction.datetime, 
                           top, 
                           bottom, 
                           alpha=0.3, 
                           zorder = 3, 
                           # color = '0.5',
                           color = col_major,
                           label = 'ci')

        if show_confidence_adjusted:
            # if transform == 'percent':
            #     top = (10**self.prediction.partial_datetime_confidence.sel(boundary='top') -1) * 100
            #     bottom = (10**self.prediction.partial_datetime_confidence.sel(boundary='bottom') -1) * 100
            # else:         
            cia = self.confidence_interval_adjusted_newey_wet
            top = cia.sel(boundary='top') + offset
            bottom = cia.sel(boundary='bottom') + offset                
            
            a.fill_between(self.prediction.datetime, 
                           top, 
                           bottom, 
                           alpha=0.3, 
                           zorder = 2, 
                           # color = '0.7', 
                           color = '#ff7f0e',
                           label = 'ci adjusted')
            
            
        if shade_seasons:
            seasons = [3,6,9,12]
            # colorlam = lambda x: x/12 - 0.1
            
            coll = []
            col = _atmcols.Color(shade_kwargs['color'],
                                # colors[0]
                                # model = 'hex'
                               )
            sats = _np.linspace(col.saturation, 0 + 0.0, 4)
            brits = _np.linspace(col.brightness, 1 - 0.0, 4)
            
            for brt,sat in zip(brits, sats):
                col.saturation = sat
                col.brightness = brt
                coll.append(col.rgb)
            
            shade_kwargs.pop('color')
            start_year = _pd.Timestamp(self.prediction.datetime.min().values).year
            end_year = _pd.Timestamp(self.prediction.datetime.max().values).year + 1
            for year in range(start_year, end_year): 
                for e, month in enumerate(seasons):
                    # start = _pd.Timestamp(year,month,1)
                    # end = start + _pd.DateOffset(month = 3)
                    # print(f'{start}, {end}')
                    # a.axvspan(start, end, color = 'black', alpha = 1/month, lw = 0)
                    start = _pd.Timestamp(year,month,1)
                    end = start + _pd.to_timedelta(31*3, 'd')
                    end = _pd.Timestamp(end.year, end.month, 1)
                    
                    # a.axvspan(start, end, color = f'{colorlam(month)}', lw = 0)
                    a.axvspan(start, end, color = coll[e], lw = 0, **shade_kwargs)
                    
                    
            a.set_xlim(self.prediction.datetime.min(), self.prediction.datetime.max())
            
            # custom_lines = [_mpllines.Line2D([0], [0], color=f'{colorlam(seasons[0])}', lw=4),
            #                 _mpllines.Line2D([0], [0], color=f'{colorlam(seasons[1])}', lw=4),
            #                 _mpllines.Line2D([0], [0], color=f'{pltcolorlam(seasons[2])}', lw=4),
            #                 _mpllines.Line2D([0], [0], color=f'{colorlam(seasons[3])}', lw=4),
            #                ]
            custom_lines = [_mpllines.Line2D([0], [0], color=coll[0], lw=4, **shade_kwargs),
                            _mpllines.Line2D([0], [0], color=coll[1], lw=4, **shade_kwargs),
                            _mpllines.Line2D([0], [0], color=coll[2], lw=4, **shade_kwargs),
                            _mpllines.Line2D([0], [0], color=coll[3], lw=4, **shade_kwargs),
                           ]
            
            a.legend(custom_lines, ['spring', 'summer', 'fall', 'winter'])
            
            noy = end_year - start_year
            base = int(_np.ceil(noy/10))
            a.xaxis.set_major_locator(_mpldates.YearLocator(base=base))
            a.xaxis.set_major_formatter(_mpldates.DateFormatter('%Y'))
            a.xaxis.set_minor_locator(_mpldates.MonthLocator(interval = base))
            a.xaxis.set_tick_params(reset = True)
            a.xaxis.tick_bottom()
            a.set_xlabel('')
        return f,a

    def plot_test_overview(self):
        f = _plt.figure()
        # f.set_figheight(f.get_figheight() * 1.5)
        f.set_figheight(7.4)
        
        #### Settyo
        if 1:
            gs = mplgrid.GridSpec(3,2, figure = f, 
                                  height_ratios= [1,1,1.5],
                                  hspace=0.35,  wspace = 0.3)
            # gs.update(
            #     #left=0.05, right=0.48, 
            #             # wspace = 0.3,
            #             # hspace = 0.5,
            # )
            a_residual = _plt.subplot(gs[0,:])
            a_acf = _plt.subplot(gs[1,:])
            a_rdist = _plt.subplot(gs[2,0])
            a_rprob = _plt.subplot(gs[2,1])
        
        aa = [a_residual, a_acf, a_rdist, a_rprob]
        
        #### Residual
        if 1:
            a = a_residual
            out_r = self.plot_test_residual(ax = a, show_std=False, gridsize=(180, 30))
        
            pc = a.get_children()[0]
            cm = _plt.cm.Oranges_r
            pc.set_cmap(cm)
        
            g = a.get_children()[1]
            g.set_color('black')
            
            a.set_ylabel('$\\log(AOD)$', labelpad=-3)
            # a.set_ylabel('Model residual', labelpad=-3)
            lim = 1.1
            a.set_ylim(-lim, lim)
        
        
            colors = _plt.cm.Oranges_r(_np.linspace(0, 1, 100))
            patch = [mpatches.Patch(facecolor=c) for c in colors]
            dummy = _plt.Line2D([0], [0], color=[1,0,0,1],
                               # marker='o', markersize=0, lw=0, ls = '',
                              )
            handles = [#dummy,
                        g,
                       patch,
                      ]
            labels =  [#'Residual represenations:',
                        'mean$_{rolling}$(residuals)',
                       r'$\rho$(residuals)',
                      ]
            
            legr = a.legend(handles, labels, ncol=len(handles), fontsize='small',
                     # title = 'Residual represenations',
                    handler_map = {list: HandlerTuple(None)},
                    # handletextpad = -5,
                    )
            
        
        ############
        #### Autocorrelation function
        if 1:
            a = a_acf
            
            self.plot_test_autocorr(ax = a, 
                                                           # acf_kwargs={'lags': 400, 'use_vlines': False},
                                                          )
            g = a.get_lines()[0]
            g.set_markersize(4)
            # g.set_marker('')
            # g.set_linestyle('-')
        
            a.set_ylim(-0.2, 1)
            a.set_title('')
            a.set_ylabel('r')
            a.set_xlabel('Time lag (days)', labelpad=0)
            a.legend(fontsize = 'small')
        
        ########
        ### distribution of residuals
        if 1:
            a= a_rdist
            outd = self.plot_test_dist_of_res(ax = a)
            
            p = a.get_children()[0]
            p.set_label('residuals')
            
            g = a.get_lines()[-1]
            g.set_linewidth(1.5)
            g.set_linestyle('--')
            g.set_label('fit')
            
            a.set_xlim(-lim, lim)
            a.set_yticklabels([])
            a.set_ylabel('Probability density')
            a.set_title('')
            a.set_xlabel('$\\log(AOD)$')
        
            a.legend(fontsize = 'small', handlelength = 1)
            
        ##############################
        #### Q-Q plot
        ##############
        if 1:
            a = a_rprob
            out = self.plot_test_qq(ax = a, gridsize=40)
            a.set_ylabel('Ordered $\\log(AOD)$ res.', labelpad=-3)
            # a.grid(False)
            a.set_title('')
        
            # leg = out[3]
        
        out = [f,aa,]
        return out
        

    
    def plot_test_residual_vs_obs(self, 
                                  ul=0.999,
                                  ll=0.001, 
                                  ax=None,
                                  gridsize = 100):
        if not isinstance(ax, type(None)):
            a = ax
            f = a.get_figure()
        else:
            f,a = _plt.subplots()
        
        pred = self.prediction
        # residual = self.data.iloc[:,0] - pred.prediction
        residual = pred.residual
    # a.scatter(gamcl.data.y_aod, residual)
        cmap = _plt.cm.gnuplot2
        cmap.set_under(color='none')
        #(xmin, xmax, ymin, ymax)
        maxi = self.data.iloc[:,0].quantile(ul)
        mini = self.data.iloc[:,0].quantile(ll)
        maxiy = residual.quantile(ul)
        miniy = residual.quantile(ll)
        extent = _np.array((mini, maxi, mini, maxi))
        extent = _np.array((mini, maxi, miniy, maxiy))
        pc = a.hexbin(self.data.iloc[:,0], residual, 
                      cmap = cmap, 
                      vmin = 0.000000001, 
                      linewidths=0.2, 
                      extent=extent,
                      gridsize = gridsize,
                     )
        f.colorbar(pc)
        
        a.set_xlabel('Observation')
        a.set_ylabel('Residual')
        a.set_title('Observation vs Residual')
        return f,a
        
    def plot_test_dist_of_res(self, bins = 50, ax = None):
        if not isinstance(ax, type(None)):
            a = ax
            f = a.get_figure()
        else:
            f,a = _plt.subplots()
        pred = self.prediction
        # residual = self.data.iloc[:,0] - pred.prediction
        residual = pred.residual
        a.hist(residual, bins = bins, density=True, alpha = 0.6)

        mu, sigma = scistats.norm.fit(residual)
        self.tp_residual = residual
        self.tp_bins = bins
        x = _np.linspace(float(residual.min()), float(residual.max()), bins*3)
        p = scistats.norm.pdf(x, mu, sigma)
        a.plot(x, p, 'k', linewidth=2)
        
        a.set_title('Distribution of residual')
        a.set_xlabel('Residual')
        a.axvline(0, color = [0,0,0,0.5], ls = '--')
        return f,a 
    
    def plot_test_qq(self, ax = None, textpos = (0.1, 0.9), sig_fig = 2, gridsize = 50):
        def format_significant(x, sig_digits=2):
            if x == 0:
                return f"{0:.{sig_digits-1}f}"  # handle zero specially
        
            if x > 1e-2 and x < 1e2:
                magnitude = int(_np.floor(_np.log10(abs(x))))  # find magnitude of the number
                factor = 10**(sig_digits - 1 - magnitude)
                pres = sig_digits - 1 - magnitude
                if pres < 0:
                    pres = 0
                return f"{round(x * factor) / factor:.{pres}f}"
            else:    
                # Get the order of magnitude of the number
                exponent = int(_np.floor(_np.log10(abs(x))))
                # Normalize the number to its significant digits
                scaled = round(x / 10**exponent, sig_digits - 1)
                # Format using LaTeX scientific notation style
                # if exponent == 0:
                
                    # return f"{scaled:.{sig_digits-1}f}"
                return r"${:.{prec}f} \cdot 10^{{{exp}}}$".format(scaled, prec=sig_digits-1, exp=exponent)
            
        if not isinstance(ax, type(None)):
            a = ax
            f = a.get_figure()
        else:
            f,a = _plt.subplots()
        pred = self.prediction
        # residual = self.data.iloc[:,0] - pred.prediction
        residual = pred.residual
        out = scistats.probplot(residual, dist="norm", plot=a)
        gs = a.get_lines()
        gs[0].remove()
        ((osm, osr),(slope, intercept, r)) = out
        self.prediction.residual.attrs['qq_intercept'] = intercept
        self.prediction.residual.attrs['qq_correlation'] = r
        # cmap = _plt.cm.gnuplot2
        cmap = _plt.cm.plasma
        cmap.set_under(color='none')
        pc = a.hexbin(osm, osr, cmap = cmap, vmin = 0.000000001, linewidths=0.2, gridsize = gridsize,
                 # zorder = 10
                )
        # a.grid()
        # if 0:
        #     a.text(*textpos, 
        #            (#f'm = {format_significant(slope, sig_fig)}\n'
        #             f'intercept $c$ = {format_significant(intercept, sig_fig)}\n'
        #             f'correlation $r$={format_significant(r, sig_fig)}'),
        #             transform = a.transAxes, 
        #             ha = 'left', va= 'top')
        
        g = a.get_lines()[-1]
        # g.set_label(('linear reg., '
        #              f'$c$ = {format_significant(intercept, sig_fig)}, '
        #              f'$r$={format_significant(r, sig_fig)}'),
        #            )
        pcolors = pc.get_cmap()(_np.linspace(0, 1, 100))
        patch = [mpatches.Patch(facecolor=c) for c in pcolors]
        
        handles=[patch,g]
        glabel = ('linear reg.\n'
                 f'$c$ = {format_significant(intercept, sig_fig)}\n'
                 f'$r$ = {format_significant(r, sig_fig)}')
        labels=[r'$\rho$(data)', glabel]
        
        leg = a.legend(handles, labels, 
                       # ncol=len(handles), 
                         fontsize='small',
                       handler_map = {list: HandlerTuple(None)},
                       loc = 4)
        col = [0,0,0,.5]
        a.axvline(0, ls = '--', color = col)
        a.axhline(0, ls = '--', color = col)
        return f,a,pc, leg
    
    def plot_test_residual(self, window = 2*365, ax = None, show_std = True, gridsize = 80):
        if not isinstance(ax, type(None)):
        #     a = ax
        #     f = a.get_figure()
        # else:
            if show_std:
                aa = ax
            else:
                aa = [ax,]
            # assert(False), "ax has to be none right now, sorry, programming required"
            f = aa[0].get_figure()
            
            # f,a = _plt.subplots()
        else:
            if show_std:
                f,aa = _plt.subplots(2, height_ratios=[2,1], gridspec_kw={'hspace':0})
    
            else:
                f, a = _plt.subplots()
                aa = [a,]
                
        # residual = self.data.iloc[:,0] - self.prediction.prediction
        residual = self.prediction.residual.to_dataframe()
        # residualdf = _pd.DataFrame(residual)
        # residualdf = residual.to_dataframe()
        
        roll = residual.rolling(window, center = True, win_type='gaussian')

       
        
        a = aa[0]
        out = a.hexbin(mdates.date2num(residual.index), residual.residual, linewidths=0.2,
                       gridsize = gridsize
                      )
        qm = out
        cmap = _plt.cm.gnuplot2
        cmap.set_under(color='none')
        # cm = _plt.cm.inferno_r
        # cm.set_under([1,1,1,0])
        qm.set_cmap(cmap)
        qm.set_clim(vmin = 0.0001)
        # a.xaxis_date()
        format_str = '%Y'
        format_ = mdates.DateFormatter(format_str)
        a.xaxis.set_major_formatter(format_)
        
        std = window/3
        roll.mean(std=std).plot(ax = a, color = 'white')

        if show_std:
            a  = aa[1]
            roll.std(std=std).plot(ax = a)
            a.set_xlim(aa[0].get_xlim())
            a.set_ylabel('Std')
            
            a = aa[0]
            a.set_xticklabels([]) 
            a.set_title('Residual time series')
            a.set_ylabel('Residual')
            
        for a in aa:
            a.legend().remove()
            a.set_xlabel('')
        
        return f,a 
    
    def plot_test_autocorr(self, 
                           lag = None,
                           # acf_kwargs = {'lags':60}, 
                           ax = None):
        """
        Parameters
        -----------
        lag: int
            Maximum lag. Note this is limited what has been created with autocorrelation_create_custom. Rerun
            autocorrelation_create_custom with adjusted parameters if neede.
        """
        if not isinstance(ax, type(None)):
            a = ax
            f = a.get_figure()
        else:
            f,a = _plt.subplots()
            
        acf = self.autocorrelation

        sel = acf.sel(acf_params = 'acf')
        sel.plot(ax = a, ls = '', marker = '.')
        g = a.get_lines()[-1]
        g.set_label('Autocorrelation function')

        a.vlines(acf.lag, 0, sel, alpha = 0.7)
        
        sel = acf.sel(acf_params = 'ci_bartlett')
        out = a.fill_between(acf.lag, sel, -sel, alpha = 0.5)
        out.set_label('Confidence interval')       
        
        a.axhline(0, ls = '--', color = '0.5')
        a.set_title('Autocorrelation function')
        a.set_ylabel('Autocorrelation $r$')
        a.set_xlabel('Lag')
        a.legend()

        
        # acf_kwargs: see https://www.statsmodels.org/stable/generated/statsmodels.graphics.tsaplots.plot_acf.html
        # residual = self.prediction.residual
        # statsmod_plot_acf(residual, ax = a, **acf_kwargs,
        #         )  # Adjust 'lags' as necessary
        # a.set_title('Autocorrelation Function')
        # a.set_ylim(auto = True)
        return f,a 
    
    def plot_test_obs_vs_pred(self, ll = 0.05, ul = 0.95, ax = None):
        if not isinstance(ax, type(None)):
            a = ax
            f = a.get_figure()
        else:
            f,a = _plt.subplots()
    
        maxi = self.data.iloc[:,0].quantile(ul)
        mini = self.data.iloc[:,0].quantile(ll)
        extent = _np.array((mini, maxi, mini, maxi))
        
        # a.scatter(gamcl.data.y_aod, residual)
        cmap = _plt.cm.gnuplot2
        cmap.set_under(color='none')
        
        pc = a.hexbin(self.data.iloc[:,0], self.prediction.prediction, 
                      cmap = cmap, 
                      vmin = 0.000000001, 
                      linewidths=0.2, 
                      gridsize=50,
                      extent = extent)
        f.colorbar(pc)
        a.set_xlabel('observation')
        a.set_ylabel('prediction')
        a.set_title('Observation vs Prediction')
        # a.set_
        return f,a 

class Climatology(object):
    def __init__(self, parent_stats = None, reference = 'index', frequency = 'M', bins = None):
        """

        Parameters
        ----------
        parent_stats
        reference: str ['index']
            if one wants to plot agains a particular column in the timeseries rather then the index.
        frequency: str (['M'], 'H')
            originally this was written for things like seasonality, diurnality, with otherwords with anyting related
            to the datetime.
        bins: array like
            This is only used if reference column (e.g. index) is not a datetime object, e.g. if other climatologic
            indeces that are not datetime related are used.
        """
        self._bins = bins
        self._reference = reference
        if not isinstance(parent_stats, type(None)):
            self._parent_stats = parent_stats
            self._parent_ts = parent_stats._parent_ts
        self._frequency = frequency
        self._reset()
        self._timezone = None

    def _reset(self):
        self._percentiles = None


    @property
    def frequency(self):
        return self._frequency

    @property
    def percentiles(self):
        if type(self._percentiles) == type(None):
            # rs = self._parent_ts.data.resample(self.frequency, label='left',
            #                                    # convention = 'start',
            #                                    # closed = 'left'
            #                                    )
            data = self._parent_ts.data.copy()

            # if another column instead of the index is to be used as a reference
            if not self._reference == 'index':
                data.index = data[self._reference]
                data.drop(labels = [self._reference], axis = 1, inplace = True)

            if type(data.index).__name__ == 'DatetimeIndex':
                if self._timezone:
                    data.index += _pd.Timedelta(self._timezone, 'h')
                if self.frequency == 'H':
                    data.index = data.index.hour
                elif self.frequency == 'M':
                    data.index = data.index.month
            else:
                data.index = _pd.Series(data.index).apply(lambda x: self._bins[abs(self._bins - x).argmin()])

            data.sort_index(inplace=True)
            rs = data.groupby(data.index)
            def percentile(x, q):
                # print(x.columns)

                # print('-----')
                x = x.dropna()
                if x.shape[0] == 0:
                    pct = _np.nan
                elif x.shape[1] == 1:
                    pct = _np.percentile(x, q)
                elif x.shape[1] == 2:
                    weights = x.weights
                    values = x.drop(['weights'], axis = 1).iloc[:,0]
                    # print(values.)
                    # print(x.columns)
                    # print([q/100.])
                    # print('--------')
                    pct = _array_tools.weighted_quantile(values, [q/100.], sample_weight = weights)[0]
                else:
                    txt = 'this currently works only for single column time series or 2 colmn if on of the columns give weights'
                    raise ValueError(txt)
                return pct

            def average(x):
                x = x.dropna()
                weights = x.weights
                values = x.drop(['weights'], axis=1).iloc[:, 0]
                avg = _np.average(values,weights = weights)
                return avg

            def median(x):
                x = x.dropna()
                cols = list(x.columns)
                # print(cols)
                cols.pop(cols.index('weights'))
                x.sort_values(cols[0], inplace=True)
                cumsum = x.weights.cumsum()
                cutoff = x.weights.sum() / 2.0
                median = x[cols[0]][cumsum >= cutoff].iloc[0]
                return median

            def number_of_valid_data_points(x):
                x = x.dropna()
                return x.shape[0]

            out = _pd.DataFrame()
            if data.shape[1] == 1:
                out['mean'] = rs.mean().iloc[:, 0]
                out['median'] = rs.median().iloc[:,0]
            elif data.shape[1] == 2:
                if 'weights' not in data.columns:
                    raise KeyError('If two columns are given one of them must have the label "weights"')
                out['mean'] = rs.apply(average)
                out['median'] = rs.apply(median)

            perc_list = [5, 25, 50, 75, 95]
            for perc in perc_list:
                outt = rs.apply(lambda x: percentile(x, perc))
                # print(outt.shape)
                out[perc] = outt#.iloc[:, 0]
            out['n_valid'] = rs.apply(number_of_valid_data_points)
            # out.index += _np.timedelta64(1, 'D')
            self._percentiles = out
        return self._percentiles
    
    @percentiles.setter 
    def percentiles(self, df):
        assert(df.index[0] == 1), f'first index should be 1 (is {df.index[0]}) for January. Make sure to use "index_col=0" when reading the csv file'
        df.columns = [int(c) if c.isnumeric() else c for c in df.columns]
        self._percentiles = df
    
    
    def plot_percentiles(self, ax=None, box_width=0.2, wisker_size=20, mean_size=10, median_size = 10 , line_width=1.5,
                         xoffset=0,
                         color=0, tickbase = 1):
        """
        This will plot the percentiles of the data ([5, 25, 50, 75, 95]. In addition, the mean (o) and median (_) are shown.

        Parameters
        ----------
        ax
        box_width
        wisker_size
        mean_size
        median_size
        line_width
        xoffset
        color
        tickbase

        Returns
        -------
        f, a, boxes, vlines, wisker_tips, mean
        """
        try:
            import plt_tools as _plt_tools
            if type(color) == int:
                color = _plt.rcParams['axes.prop_cycle'].by_key()['color'][color]
                col = _plt_tools.colors.Color(color, model='hex')
            elif type(color) == str:
                col = _plt_tools.colors.Color(color, model='hex')
            else:
                col = _plt_tools.colors.Color(color, model='rgb')
    
            col.saturation = 0.3
            color_bright = col.rgb
        except:
            color = _mcolors.hex2color(_plt.rcParams['axes.prop_cycle'].by_key()['color'][color])
            color_bright = color + (0.6,)

        if ax:
            a = ax
            f = a.get_figure()
        else:
            f, a = _plt.subplots()

        boxes = []
        vlines = []
        xordinal = []
        for row in self.percentiles.iterrows():
            #             width = 10
            x = row[0] + xoffset
            xordinal.append(x)

            # box
            # y = (row[1][75] + row[1][25]) / 2
            y = row[1][25]
            height = row[1][75] - row[1][25]
            box = _plt.Rectangle((x - box_width / 2, y), box_width, height,
                                 #                         ha = 'center'
                                 )
            box.set_facecolor([1, 1, 1, 1])
            a.add_patch(box)
            boxes.append(box)
            # wiskers
            y = (row[1][95] + row[1][5]) / 2
            vl = a.vlines(x, row[1][5], row[1][95])
            vlines.append(vl)

        for b in boxes:
            b.set_linewidth(line_width)
            b.set_facecolor(color_bright)
            b.set_edgecolor(color)
            b.set_zorder(2)

        for vl in vlines:
            vl.set_color(color)
            vl.set_linewidth(line_width)
            vl.set_zorder(1)

        # wm_lw = 2
        wisker_tips = []
        if wisker_size:
            g, = a.plot(xordinal, self.percentiles[5], ls='')
            wisker_tips.append(g)

            g, = a.plot(xordinal, self.percentiles[95], ls='')
            wisker_tips.append(g)

        for wt in wisker_tips:
            wt.set_markeredgewidth(line_width)
            wt.set_color(color)
            wt.set_markersize(wisker_size)
            wt.set_marker('_')

        mean = None
        if mean_size:
            g, = a.plot(xordinal, self.percentiles['mean'], ls='')
            g.set_marker('o')
            g.set_markersize(mean_size)
            g.set_zorder(20)
            g.set_markerfacecolor('None')
            g.set_markeredgewidth(line_width)
            g.set_markeredgecolor(color)
            mean = g

        median = None
        if median_size:
            g, = a.plot(xordinal, self.percentiles['median'], ls='')
            g.set_marker('_')
            g.set_markersize(median_size)
            g.set_zorder(20)
            g.set_markeredgewidth(line_width)
            g.set_markeredgecolor(color)
            median = g

        # a.xaxis.set_major_locator(_MonthLocator())
        # a.xaxis.set_major_formatter(_DateFormatter('%b'))
        try:
            a.set_ylim(_np.nanmin(self.percentiles.drop(['n_valid'], axis=1)),
                       _np.nanmax(self.percentiles.drop(['n_valid'], axis=1)))
        except:
            pass

        try:
            a.set_xlim(self.percentiles.index.min(), self.percentiles.index.max())
        except:
            pass

        mjl = _MultipleLocator(tickbase)
        a.xaxis.set_major_locator(mjl)




        # a.relim()
        # a.autoscale_view(tight=True)
        # f.autofmt_xdate()
        return f, a, boxes, vlines, wisker_tips, mean, median

            # def plot_percentiles(self, ax=None, box_width=10, wisker_size=20, mean_size = 10, line_width = 1.5, xoffset = 0, color=0):
    #     """
    #
    #     Parameters
    #     ----------
    #     ax
    #     box_width
    #     wisker_size
    #     mean_size
    #     color
    #
    #     Returns
    #     -------
    #     f, a, boxes, vlines, wisker_tips, mean
    #     """
    #     if type(color) == int:
    #         color = _plt.rcParams['axes.prop_cycle'].by_key()['color'][color]
    #         col = _plt_tools.colors.Color(color, model='hex')
    #     else:
    #         col = _plt_tools.colors.Color(color, model='rgb')
    #
    #     col.saturation = 0.3
    #     color_bright = col.rgb
    #
    #     if ax:
    #         a = ax
    #         f = a.get_figure()
    #     else:
    #         f, a = _plt.subplots()
    #
    #     boxes = []
    #     vlines = []
    #     xordinal = []
    #     for row in self.percentiles.iterrows():
    #         #             width = 10
    #         x = row[0].toordinal() + xoffset
    #         xordinal.append(x)
    #
    #         # box
    #         y = (row[1][75] + row[1][25]) / 2
    #         height = row[1][75] - row[1][25]
    #         box = _plt.Rectangle((x - box_width / 2, y), box_width, height,
    #                              #                         ha = 'center'
    #                              )
    #         box.set_facecolor([1, 1, 1, 1])
    #         a.add_patch(box)
    #         boxes.append(box)
    #         # wiskers
    #         y = (row[1][95] + row[1][5]) / 2
    #         vl = a.vlines(x, row[1][5], row[1][95])
    #         vlines.append(vl)
    #
    #     for b in boxes:
    #         b.set_linewidth(line_width)
    #         b.set_facecolor(color_bright)
    #         b.set_edgecolor(color)
    #         b.set_zorder(2)
    #
    #     for vl in vlines:
    #         vl.set_color(color)
    #         vl.set_linewidth(line_width)
    #         vl.set_zorder(1)
    #
    #     # wm_lw = 2
    #     wisker_tips = []
    #     if wisker_size:
    #         g, = a.plot(xordinal, self.percentiles[5], ls='')
    #         wisker_tips.append(g)
    #
    #         g, = a.plot(xordinal, self.percentiles[95], ls='')
    #         wisker_tips.append(g)
    #
    #     for wt in wisker_tips:
    #         wt.set_markeredgewidth(line_width)
    #         wt.set_color(color)
    #         wt.set_markersize(wisker_size)
    #         wt.set_marker('_')
    #
    #     mean = None
    #     if mean_size:
    #         g, = a.plot(xordinal, self.percentiles['mean'], ls = '')
    #         g.set_marker('o')
    #         g.set_markersize(mean_size)
    #         g.set_zorder(20)
    #         g.set_markerfacecolor('None')
    #         g.set_markeredgewidth(line_width)
    #         g.set_markeredgecolor(color)
    #         mean = g
    #
    #     a.xaxis.set_major_locator(_MonthLocator())
    #     a.xaxis.set_major_formatter(_DateFormatter('%b'))
    #
    #     a.set_ylim(_np.nanmin(self.percentiles.drop(['n_valid'], axis=1)), _np.nanmax(self.percentiles.drop(['n_valid'], axis=1)))
    #     a.set_xlim(self.percentiles.index.min().toordinal(), self.percentiles.index.max().toordinal())
    #     # a.relim()
    #     # a.autoscale_view(tight=True)
    #     # f.autofmt_xdate()
    #     return f, a, boxes, vlines, wisker_tips, mean

class Seasonality(Climatology):
    def plot_percentiles(self, *args, **kwargs):
        out = super().plot_percentiles(*args, **kwargs)
        a = out[1]
        def num2month(pos, num):
            month = ['', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', '', '']
            return month[num]

        if self.frequency == 'M':
            a.xaxis.set_major_formatter(_FuncFormatter(num2month))
            a.set_xlim(0.5, 12.5)
        a.set_xlabel('Month of year')
        return out

class Diurnality(Climatology):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, frequency='H')

    def plot_percentiles(self, *args, **kwargs):
        out = super().plot_percentiles(*args, **kwargs)
        a = out[1]

        if self.frequency == 'H':
            a.set_xlim(-0.5,23.5)

        a.set_xlabel('Hour of day')
        return out

    @property
    def timezone(self):
        return self._timezone

    @timezone.setter
    def timezone(self, value):
        self._reset()
        self._timezone = value

def create_lambda_set(no_col = 2, elim = (-3,3), size = 1000):
    # elim = (-0,2)
    lams = _np.random.rand(size, no_col) # random points on [0, 1], with shape (100, 3)
    lams = lams * (elim[1] - elim[0]) + elim[0] # shift values to -3, 3
    lams = 10 ** lams # transforms values to 1e-3, 1e3
    # f,a = plt.subplots()
    # out = a.hist(lams, bins = np.logspace(*elim, 20))
    # a.set_xscale('log')
    return lams