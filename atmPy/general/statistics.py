import pandas as _pd
import matplotlib.pylab as _plt
import numpy as _np
# from matplotlib.dates import MonthLocator  as _MonthLocator
# from matplotlib.dates import DateFormatter as _DateFormatter
from matplotlib.ticker import FuncFormatter as _FuncFormatter
from matplotlib.ticker import MultipleLocator as _MultipleLocator
import matplotlib.colors as _mcolors
#import plt_tools as _plt_tools
from atmPy.tools import array_tools as _array_tools
import atmPy.tools.plt_tool_kit.colors as _atmcols
# from pygam import LinearGAM, pygam_s, pygam_l

import datetime
import xarray as _xr 
import matplotlib.dates as mdates
import  matplotlib.lines as _mpllines
import matplotlib.dates as _mpldates

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
        self._splines_per_year = 15
        self._parent_stats = parent_stats
        self._data = None
        self._seasonality_lam = None
        self._seasonality_nsplines = None
        self._trend_lam = None
        self._trend_nsplines = None
        self._prediction_grid_size = None
        self._prediction_confidence = None
        
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
            
            row = data.iloc[0,:]
            start_date = row.name
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
                
            y.name = f'y_{y.name}'
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
    def gam_inst(self, linear = False, ):
        if isinstance(self._gam, type(None)):
            import pygam
            
            if linear:
                year = pygam.l(0)#, lam = resolution, n_splines=int(n_splines))
            else:
                year = pygam.s(0, lam = self.trend_lambda, n_splines=self.trend_nsplines)#int(n_splines))
            
            doy = pygam.s(1, basis = 'cp', lam = self.seasonality_lambda, n_splines= self.seasonality_nsplines)
            if self.distribution == 'normal':
                self._gam = pygam.LinearGAM(year + doy)
            else:
                self._gam = pygam.GAM(year + doy , distribution = self.distribution, link = self.link)
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
                
            XX = self.data.iloc[:,1:].values
            pred = gam.predict(XX)  
            dindex = self._start_date + _pd.to_timedelta(XX[:,0], 'd')
            ds['prediction'] = _xr.DataArray(pred, coords = {'datetime_full': dindex})
            
            if hasattr(gam, 'prediction_intervals'):
                pred = gam.prediction_intervals(XX)
                ds['prediction_confidence'] = _xr.DataArray(pred, coords = {'datetime_full': dindex,
                                                                            'boundary': ['top', 'bottom']})
            self._prediction = ds
        return self._prediction
    
    @prediction.setter 
    def prediction(self, value):
        self._prediction = value
    
    def plot_seasonality(self, ax = None,
                         offset = 0,
                         xticklablesmonth = True,
                         show_confidence = True,
                         **plot_kwargs):
        
        if isinstance(ax, type(None)):
            f, a = _plt.subplots()
        else:
            a = ax
            f = a.get_figure()
        
        # if 'label' not in plot_kwargs:
        #     plot_kwargs['label'] = 'seasonal'
        
        (self.prediction.partial_doy + offset).plot(ax = a, **plot_kwargs)   
        
        if show_confidence:
            bot,up = self.prediction.partial_doy_confidence.values.transpose()
            fill = a.fill_between(self.prediction.doy,bot,up, color = [0,0,0,0.3])
        else:
            fill = None
        
        if xticklablesmonth:
            a.xaxis.set_major_locator(mdates.MonthLocator())  # Set the major ticks to be at the beginning of each month
            a.xaxis.set_minor_locator(mdates.WeekdayLocator())  # Set the minor ticks to be at the beginning of each week
            a.xaxis.set_major_formatter(mdates.DateFormatter('%b'))  # Format the major ticks with abbreviated month names
            
            a.xaxis.set_minor_locator(_plt.NullLocator())
        a.set_xlabel('')
        # a.set_xlabel('Day of year')
        return f,a,fill
    
    def plot_prediction(self, ax = None, show_confidence=True, show_original_data=True, **plot_kwargs):
        if isinstance(ax, type(None)):
            f, a = _plt.subplots()
        else:
            a = ax
            f = a.get_figure()

    
        if 'label' not in plot_kwargs:
            plot_kwargs['label'] = 'prediction'
        

        self.prediction.prediction.plot(ax = a, zorder = 3, label = 'prediction')
        if show_confidence:
            a.fill_between(self.prediction.datetime_full, 
                           self.prediction.prediction_confidence.sel(boundary='top'), 
                           self.prediction.prediction_confidence.sel(boundary='bottom'), 
                           alpha=0.5, zorder = 2, color = '0.5', label = 'confidence')
        if show_original_data:
            data = self._parent_stats._parent_ts.data
            a.plot(data.index, data.iloc[:,0], ls = '', marker = '.', markersize = 1, zorder = 1, label = 'observation', color = '0.7')
        return f,a
        
    def plot_overview(self, axis = None, show_confidence = True, show_original_data = True):
        if isinstance(axis, type(None)):
            aa = []
            f = _plt.figure()
            f.set_figheight(f.get_figheight() * 1.5)
            aa.append(f.add_subplot(3,1,1))
            aa.append(f.add_subplot(3,1,2, sharex = aa[0]))
            aa.append(f.add_subplot(3,1,3))
            shade_seasons = {'color': '#02401A', 'alpha': 0.5}
        else:
            aa = axis
            f = aa[0].get_figure()
            shade_seasons = False
            
        xshift = 0
        
        self.plot_prediction(ax = aa[0], show_confidence=show_confidence, show_original_data=show_original_data)
        self.plot_trend(ax=aa[1], shade_seasons=shade_seasons, show_confidence=show_confidence)
        self.plot_seasonality(ax=aa[2], show_confidence = show_confidence)
        
        
        # a = aa[0]
        # fillb = a.get_children()[2]
        # fillb.set_color('0.3')
        # fillb.set_zorder(10)
        for e,a in enumerate(aa):
            if e == 0:
                continue
            poslast = aa[e-1].get_position().bounds
            pos = list(a.get_position().bounds)
            if e == 1:
                xshift += poslast[1] - (pos[1] + pos[3])
            pos[1] = pos[1] + xshift
            a.set_position(pos)
        aa[0].legend(fontsize = 'small')
        
        return f,aa
    
    def plot_trend(self, ax = None, 
                   shade_seasons = False, 
                   show_confidence = True,
                   offset= 0, **plot_kwargs):
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
            plot_kwargs['zorder'] = 3
        
        if isinstance(ax, type(None)):
            f, a = _plt.subplots()
            
        else:
            a = ax
            f = a.get_figure()
            
 
        (self.prediction.partial_datetime + offset).plot(ax = a, **plot_kwargs)   
         

            
        if show_confidence:
            a.fill_between(self.prediction.datetime, 
                   self.prediction.partial_datetime_confidence.sel(boundary='top'), 
                   self.prediction.partial_datetime_confidence.sel(boundary='bottom'), 
                   alpha=0.5, zorder = 2, color = '0.5', label = 'confidence')
            
            
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
            #                 _mpllines.Line2D([0], [0], color=f'{colorlam(seasons[2])}', lw=4),
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