from atmPy.general import measurement_site as _measurement_site
import pandas as _pd
import numpy as _np
from atmPy.radiation import solar as _solar
from atmPy.general import timeseries as _timeseries
import multiprocessing as _mp
import xarray as _xr
import matplotlib.pyplot as _plt
import scipy as _sp

_colors = _plt.rcParams['axes.prop_cycle'].by_key()['color']

class AOD_AOT(object):
    def __init__(self,
                 wavelengths = None,
                 site = None,
                 lat = None,
                 lon = None,
                 elevation = 0,
                 name = None,
                 name_short = None,
                 timezone = 0,
                 site_info = None):
        """This class is for column AOD or AOT meausrements at a fixed site. This class is most usfull for aerosol
        optical properties from a CIMEL (AERONET) or a MFRSR (SURFRAD)

        Parameters
        ----------
        wavelengths: dict
            Column names are often not reflecting the precise wavelength in the channel, but the typical wavelength.
            The dictionary translates column names to exact wavelength. If AOD is calculated and wavelengths is set
            wavelengths from this will be used instead of column names.
        site: atmPy.general.station instance
        lat, lon: location of site
            lat: deg north, lon: deg east
            e.g. lat = 40.05192, lon = -88.37309 for Bondville, IL, USA
        elevation: float
            elevation of site in meter (not very importent parameter, keeping it at zero is usually fine ...
            even vor Boulder
        name, name_short: str
            name and abbriviation of name for site
        timezone: int
            Timezon in houres, e.g. Bondville has the timezone -6.
            Note, this can lead to some confusion if the data is in UTC not in local time.... this needs to be improved
            ad some point"""

        self._aot = None
        self._aod = None
        self._sunposition = None
        self._timezone = timezone
        self.wavelengths = wavelengths

        if not isinstance(site, type(None)):
            self.site = site
        elif not isinstance(lat, type(None)):
            self.site = _measurement_site.Station(lat, lon, elevation, name=name, abbreviation=name_short, info = site_info)


        self.cloudmask = CloudDetection(self)

    @property
    def sun_position(self):
        if not self._sunposition:
            if self._timezone != 0:
                date = self._timestamp_index +  _pd.to_timedelta(-1 * self._timezone, 'h')
            else:
                date = self._timestamp_index
            self._sunposition = _solar.get_sun_position(self.site.lat, self.site.lon, date)
            self._sunposition.index = self._timestamp_index
            self._sunposition = _timeseries.TimeSeries(self._sunposition)
        return self._sunposition

    @property
    def AOT(self):
        if not self._aot:
            if not self._aod:
                raise AttributeError('Make sure either AOD or AOT is set.')
            aot = self.AOD.data.mul(self.sun_position.data.airmass, axis='rows')
            aot.columns.name = 'AOT@wavelength(nm)'
            aot = _timeseries.TimeSeries(aot)
            self._aot = aot
        return self._aot

    @ AOT.setter
    def AOT(self,value):
        self._aot = value
        self._aot.data.columns.name = 'AOT@wavelength(nm)'
        self._timestamp_index = self._aot.data.index

    @property
    def AOD(self):
        if not self._aod:
            if not self._aot:
                raise AttributeError('Make sure either AOD or AOT is set.')
            aod = self.AOT.data.div(self.sun_position.data.airmass, axis='rows')
            aod.columns.name = 'AOD@wavelength(nm)'
            aod = _timeseries.TimeSeries(aod)
            self._aod = aod
        return self._aod

    @ AOD.setter
    def AOD(self,value):
        self._aod = value
        self._aod.data.columns.name = 'AOD@wavelength(nm)'
        self._timestamp_index = self._aod.data.index

    @property
    def ang_exp(self):
        return self._ang_exp

    @ang_exp.setter
    def ang_exp(self, value):
        self._ang_exp = value

    def aod2angstrom_exponent(self, column_1=500, column_2=870,
                              use_wavelength_from_column_names = None,
                              # wavelength_1=None, wavelength_2=None
                              ):
        """
        Calculates the angstrom exponents based on the AOD data.

        Parameters
        ----------
        column_1: type of column name
            column name of one of the two points used for the AOD calculation
        column_2: type of column name
            column name of the other of the two points used for the AOD calculation
        use_wavelength_from_column_names: bool [None]
            When the wavelength dictionary is set. Wavelengths from the dictionary are used instead of column names.
            Set this kwarg to True to ignore the wavelengths dictionary and use column names instead.

        Parameters (deprecated)
        -----------------------
        wavelength_1: float
            if the column name of column_1 is not accurate enough set the wavelenth used to calculate AOD here.
        wavelength_2: float
            as above for column_2

        Returns
        -------

        """
        if isinstance(self.wavelengths, type(None)) or use_wavelength_from_column_names:
            # if wavelength_1 == None:
            wavelength_1 = column_1
            # if wavelength_2 == None:
            wavelength_2 = column_2
        else:
            wavelength_1 = self.wavelengths[column_1]
            wavelength_2 = self.wavelengths[column_2]
        c1 = column_1
        c2 = column_2
        c1ex = wavelength_1
        c2ex = wavelength_2
        out = - _np.log10(self.AOD.data.loc[:, c1] / self.AOD.data.loc[:, c2]) / _np.log10(c1ex / c2ex)
        out = _timeseries.TimeSeries(_pd.DataFrame(out))
        setattr(self, 'ang_exp_{}_{}'.format(column_1, column_2), out)
        return out

    def cal_deviations(arglist):
        #     print(arglist)
        [idx, aod], [idxt, wls] = arglist

        #         wls = sfr.wavelength_matching_table.data.loc[idx,:]
        assert (len(wls) == len(aod))

        sfr.ang_exp.data.loc[idx, :]

        fitres = sp.stats.linregress(np.log10(wls.loc[[500, 870]]), np.log10(aod.loc[[500, 870]]))

        ff = lambda x: 10 ** (np.log10(x) * fitres.slope + fitres.intercept)

        aod_fit = ff(wls)
        diff = (aod - aod_fit)
        return diff

    def calculate_deviation_from_angstrom(sfr, wl1=500, wl2=870):
        """
        calcu
        Parameters
        ----------
        wl1
        wl2

        Returns
        -------

        """
        pool = _mp.Pool(6)

        # out = pool.map(cal_deviations, sfr.AOD.data.iloc[135:140].iterrows())
        out = pool.map(cal_deviations, zip(sfr.AOD.data.iterrows(), sfr.wavelength_matching_table.data.iterrows()))

        # make dataframe from output
        df = pd.concat(out, axis=1)
        df = df.transpose()

        # xarray does not like timezone ... remove it
        df.index = pd.to_datetime(df.index.astype(str).str[:-6])

        # xarray does not like integer column ames
        df.columns = df.columns.astype(str)

        # other stuff
        df.sort_index(inplace=True)
        df.index.name = 'datetime'
        return df

class CloudDetection(object):
    def __init__(self, parent):
        self.parent = parent
        self._cloudmask_AOD = None
        self._masked_data_AOD = None
        self._cloud_classifyiers_AOD = None
        self._cloudmask_angstrom = None
        self._masked_data_angstrom = None
        self._cloud_classifyiers_angstrom = None
        self._cloudmask_nativ = None
    
    @property
    def cloudmask_AOD(self):
        if isinstance(self._cloudmask_AOD, type(None)):
            self.get_custom_cloudmask_AOD()
        return self._cloudmask_AOD        
    
    @property
    def masked_data_AOD(self):
        if isinstance(self._masked_data_AOD, type(None)):
            self.get_custom_cloudmask_AOD()
        return self._masked_data_AOD   
    
    @property
    def cloud_classifyiers_AOD(self):
        if isinstance(self._cloud_classifyiers_AOD, type(None)):
            self.get_custom_cloudmask_AOD()
        return self._cloud_classifyiers_AOD    
        
        
    def get_custom_cloudmask_AOD(self, 
                                      mean_discriminator = 0.035,
                                      mean_window = 3,
                                      linreg_discriminator = 0.00006,
                                      linreg_window = 4,
                                      deriv_discriminator = 5e-4,
                                      min_consec_valid = 5,
                                     ):
        

        settings = _pd.DataFrame({'mean':{'discriminator': mean_discriminator, 'window': mean_window, 'min_consec_valid': min_consec_valid},
                                  'linreg':{'discriminator': linreg_discriminator, 'window': linreg_window, 'min_consec_valid': min_consec_valid},
                                  'deriv':{'discriminator': deriv_discriminator, 
                                           # 'window': linreg_window,
                                           'min_consec_valid': min_consec_valid}})
       
        
        self._setting_AOD = settings

        data = self.parent.AOD.data
        channels = data.columns.values.astype(int)

        masked_data = _xr.DataArray(coords= {'channel':channels,
                                            'datetime': data.index.values,
                                            'criteria': ['mean', 'linreg', 'deriv']}, dims = ['channel', 'datetime', 'criteria'])
        cloudmask = _xr.DataArray(coords= {'channel':channels,
                                            'datetime': data.index.values,
                                            'criteria': ['mean', 'linreg', 'deriv']}, dims = ['channel', 'datetime', 'criteria'])

        cloud_classifyiers = _xr.DataArray(coords= {'channel':channels,
                                                   'datetime': data.index.values,
                                                   'criteria': ['mean', 'linreg', 'deriv'],
                                                   'direction': ['forward', 'backward']}, dims = ['channel', 'datetime', 'criteria', 'direction'])

        for ch, data_ncs in data.iteritems():

            for cn, ser in settings.iteritems():
                # get classifyers ... both directions
                rstd = self.get_classifyer(data_ncs, ctype = cn, window = ser.window)['classifyer']
                cloud_classifyiers.loc[dict(channel = ch, criteria = cn, direction = 'forward')] = rstd

                rstd_r = self.get_classifyer(data_ncs[::-1], ctype = cn, window = ser.window)['classifyer']
                cloud_classifyiers.loc[dict(channel = ch, criteria = cn, direction = 'backward')] = rstd_r

                # cloudmasked by discriminator
                cloud_mask_t = _np.logical_and(rstd > ser.discriminator, rstd_r > ser.discriminator)
                
                # add to cloudmasked by minimum number of points
                cloud_mask_t = self.check_min_consec_valid(cloud_mask_t, min_consec_valid=ser.min_consec_valid)
                cloudmask.loc[dict(channel = ch, criteria = cn)] = cloud_mask_t
                
                # apply coudmask to data, is can probably go?!?
                data_cs_new_j = data_ncs.copy()
                data_cs_new_j[cloud_mask_t == 1] = _np.nan
                masked_data.loc[dict(channel = ch, criteria = cn)] = data_cs_new_j

        self._cloudmask_AOD = cloudmask
        self._masked_data_AOD = masked_data
        self._cloud_classifyiers_AOD = cloud_classifyiers

    def plot_classifyer_AOD(self, criteria = 'linreg', direction = 'forward', channel = None, ax = None):
        if isinstance(ax, type(None)):
            a = _plt.subplot()
        else:
            a = ax
        # df = self.cloud_classifyiers_AOD.sel(criteria = criteria, direction = direction).to_pandas().transpose()
        df = self.cloud_classifyiers_AOD.sel(criteria = criteria, direction = direction)
        if not isinstance(channel, type(None)):
             df = df.sel(channel = channel)
        df = df.to_pandas().transpose()

        
        df.plot(ax = a)
        for g in a.get_lines():
            g.set_linestyle('')
            g.set_marker('.')
        a.set_title(f'classifyer - criteria: {criteria}, direction: {direction}')
        a.set_xlabel('')
        a.axhline(self._setting_AOD.loc['discriminator', criteria], color = 'black', ls = '--')
        a.set_yscale('log')
        a.legend(fontsize = 'small', title = 'channel (nm)')
        return a
        
    @property
    def cloudmask_angstrom(self):
        if isinstance(self._cloudmask_angstrom, type(None)):
            self.get_custom_cloudmask_angstrom()
        return self._cloudmask_angstrom        
    
    @property
    def masked_data_angstrom(self):
        if isinstance(self._masked_data_angstrom, type(None)):
            self.get_custom_cloudmask_angstrom()
        return self._masked_data_angstrom   
    
    @property
    def cloud_classifyiers_angstrom(self):
        if isinstance(self._cloud_classifyiers_angstrom, type(None)):
            self.get_custom_cloudmask_angstrom()
        return self._cloud_classifyiers_angstrom    
        
        
    def get_custom_cloudmask_angstrom(self, 
                                      mean_discriminator = 0.1,
                                      mean_window = 3,
                                      linreg_discriminator = 20e-05,
                                      deriv_discriminator = 5e-4,
                                      linreg_window = 4,
                                      min_consec_valid = 5,
                                     ):
        

        settings = _pd.DataFrame({'mean':{'discriminator': mean_discriminator, 'window': mean_window, 'min_consec_valid': min_consec_valid},
                                  'linreg':{'discriminator': linreg_discriminator, 'window': linreg_window, 'min_consec_valid': min_consec_valid},
                                  'deriv':{'discriminator': deriv_discriminator, 
                                           # 'window': linreg_window,
                                           'min_consec_valid': min_consec_valid}
                                  })
        self._setting_angstrom = settings
        # settings

        # cloudmask = pd.DataFrame()
        # cloud_classifyiers = pd.DataFrame()
        # masked_data = pd.DataFrame()

        data_ncs = self.parent.ang_exp.data.iloc[:,0]
#         channels = data.columns.values.astype(int)

        masked_data = _xr.DataArray(coords= {#'channel':channels,
                                            'datetime': data_ncs.index.values,
                                            'criteria': ['mean', 'linreg', 'deriv']}, dims = ['datetime', 'criteria'])
        cloudmask = _xr.DataArray(coords= {#'channel':channels,
                                            'datetime': data_ncs.index.values,
                                            'criteria': ['mean', 'linreg', 'deriv']}, dims = ['datetime', 'criteria'])

        cloud_classifyiers = _xr.DataArray(coords= {#'channel':channels,
                                                   'datetime': data_ncs.index.values,
                                                   'criteria': ['mean', 'linreg', 'deriv'],
                                                   'direction': ['forward', 'backward']}, dims = ['datetime', 'criteria', 'direction'])

#         for ch, data_ncs in aod_ncs.data.iteritems():
        #     break

        for cn, ser in settings.iteritems():
            # get the values to be juged
            rstd = abs(self.get_classifyer(data_ncs, ctype = cn, window = ser.window)['classifyer'])
            cloud_classifyiers.loc[dict(criteria = cn, direction = 'forward')] = rstd

            rstd_r = abs(self.get_classifyer(data_ncs[::-1], ctype = cn, window = ser.window)['classifyer'])
            rstd_r = rstd_r[::-1]
            cloud_classifyiers.loc[dict(criteria = cn, direction = 'backward')] = rstd_r

            cloud_mask_t = _np.logical_and(rstd > ser.discriminator, rstd_r > ser.discriminator) 
#             self.tp_mask = cloud_mask_t.copy()
#             self.tp_disc = ser.discriminator
#             self.tp_rstd = rstd.copy()
#             self.tp_rstd_r = rstd_r.copy()
#             self.tp_cn = cn
#             self.tp_data_ncs = data_ncs
            
            cloud_mask_t = self.check_min_consec_valid(cloud_mask_t,  min_consec_valid=ser.min_consec_valid)
            cloudmask.loc[dict(criteria = cn)] = cloud_mask_t

            data_cs_new_j = data_ncs.copy()
            data_cs_new_j[cloud_mask_t == 1] = _np.nan
            masked_data.loc[dict(criteria = cn)] = data_cs_new_j
            # minimum number of points between cloud detection?
#             data_cs_new_j_mcv = check_min_consec_valid(data_cs_new_j, min_consec_valid=ser.min_consec_valid)
        #     masked_data[cn] = data_cs_new_j_mcv

        #     cloudmask[cn] = cloud_mask_t
                
        self._cloudmask_angstrom = cloudmask
        self._masked_data_angstrom = masked_data
        self._cloud_classifyiers_angstrom = cloud_classifyiers

    def plot_classifyer_angstrom(self, criteria = 'linreg', direction = 'forward', ax = None):
        if isinstance(ax, type(None)):
            a = _plt.subplot()
        else:
            a = ax
        df = self.cloud_classifyiers_angstrom.sel(criteria = criteria, direction = direction).to_pandas()#.transpose()
        
        df.plot(ax = a)
        for g in a.get_lines():
            g.set_linestyle('')
            g.set_marker('.')
        a.set_title(f'classifyer - criteria: {criteria}, direction: {direction}')
        a.set_xlabel('')
        a.axhline(self._setting_angstrom.loc['discriminator', criteria], color = 'black', ls = '--')
        a.set_yscale('log')
        return a
    
    def plot_masked_data(self, data2plot = 'AOD', classifyer = False, data2discriminate = 'AOD', channel_plot = 500, channel_discriminate = None, criteria = 'mean', show_unmasked = True, ax = None, marker_size_scale = 2):
        """
        

        Parameters
        ----------
        data2plot : TYPE, optional
            DESCRIPTION. The default is 'AOD'.
        classifyer : TYPE, optional
            DESCRIPTION. The default is False.
        data2discriminate : TYPE, optional
            Which data to the classivication is based on.The default is 'AOD'.
            AOD:
            angstrom:
            nativ: this will apply the cloudmask that comes with the data product
        channel_plot : TYPE, optional
            DESCRIPTION. The default is 500.
        channel_discriminate : TYPE, optional
            DESCRIPTION. The default is None.
        criteria : TYPE, optional
            DESCRIPTION. The default is 'mean'.
        show_unmasked : TYPE, optional
            DESCRIPTION. The default is True.
        ax : TYPE, optional
            DESCRIPTION. The default is None.
        marker_size_scale : TYPE, optional
            DESCRIPTION. The default is 2.

        Returns
        -------
        None.

        """
        if isinstance(ax, type(channel_discriminate)):
            channel_discriminate = channel_plot
            
        if isinstance(ax, type(None)):
            a = _plt.subplot()
        else:
            a = ax

        if data2discriminate == 'AOD':
            cloudmask = self.cloudmask_AOD.sel(criteria = criteria, channel = channel_discriminate).values == 1
    
        elif data2discriminate == 'angstrom':
            cloudmask = self.cloudmask_angstrom.sel(criteria = criteria).values == 1

        elif data2discriminate == 'nativ':
            cloudmask = self.cloudmask_nativ.data.values[:,0] == 1
        
        else:
            assert(False), f'{data2discriminate} not an option for kwarg data2discriminate'

        # print(f'cloudmask.shape: {cloudmask.shape}')

        if data2plot == 'AOD':
            if classifyer:
                data_ncs = self.cloud_classifyiers_AOD.sel(criteria = criteria,  channel = channel_discriminate, direction = classifyer).to_pandas()
            else:
                data_ncs = self.parent.AOD.data.loc[:,channel_plot]
            ylabel = f'AOD@{channel_plot}'
    #         data_cs = self.masked_data_AOD.sel(criteria = criteria, channel = channel).to_pandas().transpose()

        elif data2plot == 'angstrom':
            if classifyer:
                data_ncs = self.cloud_classifyiers_angstrom.sel(criteria = criteria, direction = classifyer).to_pandas()
            else:
                data_ncs = self.parent.ang_exp.data
            ylabel = '$\\AA$'
    #         data_cs = self.masked_data_angstrom.sel(criteria = criteria).to_pandas().transpose()
        else:
            assert(False), f'{data2plot} not an option for kwarg data2plot'

    #     return data_ncs, cloudmask
        #apply the cloudmask
        data_cs = data_ncs.copy()
        data_cs[cloudmask] = _np.nan

        if show_unmasked:
            data_ncs.plot(ax = a)

    #     data = self.parent.AOD.data.loc[:,500].copy()
    #     data

        data_cs.plot(ax = a)
        # data_cs.plot(ax = a)
        # data_cs_new_j.plot(ax = a)
        # data_cs_new_j_mcv.plot(ax = a)
        for e,g in enumerate(a.get_lines()):
            g.set_ls('')
            g.set_marker('.')
            ms = 8
            g.set_markersize(ms - (marker_size_scale * e))
            g.set_markerfacecolor(_colors[e])
            g.set_markeredgewidth(0)

        title = f'Cloudmask based on {data2discriminate} and {criteria}'
        a.set_title(title)
        a.set_ylabel(ylabel)
        a.set_xlabel('')
        return a
    
    @staticmethod
    def get_classifyer(data, window = 3, ctype = 'mean'):
        """
        

        Parameters
        ----------
        data : TYPE
            DESCRIPTION.
        window : TYPE, optional
            DESCRIPTION. The default is 3.
        ctype : TYPE, optional
            What classifying technique to use. 
            linreg: stderror of the linear regression for a rolling window.
            mean: std of the mean for a rolling window
            deriv: takes the derivative of the data (window is ignored)

        Returns
        -------
        out : TYPE
            DESCRIPTION.

        """
        out = {}
        if ctype in ['mean', 'linreg']:
            roll = data.rolling(int(window), 
                                # min_periods=int(window-1),
                                   ) 
        
            rm = roll.mean()
            if ctype == 'linreg':
                stderr = roll.apply(lambda grp: _sp.stats.linregress(abs(grp.index - grp.index[0]).seconds, grp.values).stderr)
                rstd = stderr/rm
            elif ctype == 'mean':
                rstd = roll.std()/rm
            
        elif ctype == 'deriv':
            # data = mfrsr.AOD.data#.loc[:, 500]
            rstd = data.diff().div(data.index.to_series().diff().dt.total_seconds(), axis = 'index').abs()
            # rstd /= data
            roll = None
        else:
            assert(False), f'ctype = {ctype} is not an option'
        
        out['classifyer'] = rstd.copy()
        out['roll'] = roll
        return out
    
    @staticmethod
    def check_min_consec_valid(df, min_consec_valid = 5):
        df = df.copy()
        df_no_consec = df.copy()
        df_no_consec[:] = 0
        
    #     group = df.notna().groupby(df.isna().cumsum())
        group = df.groupby((df != 0).cumsum())
        
        # df.groupby()
    
        for grpidx, grp in group:
            df_no_consec.loc[grp.index] = grp.shape[0]
        #     break
        df[df_no_consec <= min_consec_valid] = 1.
        return df
    #         leg = a.legend(fontsize = 'small', title = 'channel (nm)')