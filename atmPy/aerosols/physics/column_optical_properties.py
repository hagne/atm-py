from atmPy.general import measurement_site as _measurement_site
import pandas as _pd
import numpy as _np
from atmPy.radiation import solar as _solar
from atmPy.general import timeseries as _timeseries
import atmPy.aerosols.physics.optical_properties as _atmop
import atmPy.aerosols.size_distribution.sizedistribution as _atmsd
# import multiprocessing as _mp
import xarray as _xr
import matplotlib.pyplot as _plt
import scipy as _sp
import pathlib as _pl
import xarray as xr
import statsmodels.nonparametric.smoothers_lowess as smlowess

_colors = _plt.rcParams['axes.prop_cycle'].by_key()['color']

class AOD_AOT(object):
    def __init__(self,
                 data = None,
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
            If site is given the remaining site relevant kwargs (lat, lon, elevation, name, name_short) are ignored.
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
        # self._aod = None
        data = data.copy()
        
        if isinstance(data, xr.core.dataarray.DataArray):
            ds = xr.Dataset({'aod': data})
        elif isinstance(data, xr.core.dataset.Dataset):
            ds = data
        elif isinstance(data, _pd.core.frame.DataFrame):
            data.columns.name = 'channel'
            data.index.name = 'datetime'
            ds = xr.Dataset({'aod': data})
        else:
            raise TypeError(f'"{type(data).__name__}" is an invalid type for argument data. It should be: xr.Dataset, xr.DataArray, or pd.DataFrame')
        
        self.dataset = ds
        self._sunposition = None
        self._timezone = timezone
        self.wavelengths = wavelengths

        if not isinstance(site, type(None)):
            self.site = site
        elif not isinstance(lat, type(None)):
            self.site = _measurement_site.Station(lat, lon, elevation, name=name, abbreviation=name_short, info = site_info)


        self.cloudmask = CloudDetection(self)
        self.aerosolmask = AerosolAndCloudDetection(self)

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

    # @property
    # def AOT(self):
    #     if not self._aot:
    #         if not self._aod:
    #             raise AttributeError('Make sure either AOD or AOT is set.')
    #         aot = self.AOD.data.mul(self.sun_position.data.airmass, axis='rows')
    #         aot.columns.name = 'AOT@wavelength(nm)'
    #         aot = _timeseries.TimeSeries(aot)
    #         self._aot = aot
    #     return self._aot

    # @ AOT.setter
    # def AOT(self,value):
    #     self._aot = value
    #     self._aot.data.columns.name = 'AOT@wavelength(nm)'
    #     self._timestamp_index = self._aot.data.index

    # @property
    # def AOD(self):
    #     if not self._aod:
    #         if not self._aot:
    #             raise AttributeError('Make sure either AOD or AOT is set.')
    #         aod = self.AOT.data.div(self.sun_position.data.airmass, axis='rows')
    #         aod.columns.name = 'AOD@wavelength(nm)'
    #         aod = _timeseries.TimeSeries(aod)
    #         self._aod = aod
    #     return self._aod

    # @ AOD.setter
    # def AOD(self,value):
    #     if isinstance(value, type(None)):            
    #         self._aod = value
    #         return
    #     elif isinstance(value, _pd.DataFrame):
    #         value = _timeseries.TimeSeries(value)
            
    #     self._aod = value
    #     self._aod.data.columns.name = 'AOD@wavelength(nm)'
    #     self._timestamp_index = self._aod.data.index

    @property
    def ang_exp(self):
        return self._ang_exp

    @ang_exp.setter
    def ang_exp(self, value):
        self._ang_exp = value
        
        
    @property
    def aod(self):
        return self.dataset.aod
    
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
        
        # out = - _np.log10(self.AOD.data.loc[:, c1] / self.AOD.data.loc[:, c2]) / _np.log10(c1ex / c2ex)
        aod1 = self.dataset.aod.sel(channel = c1)
        aod2 = self.dataset.aod.sel(channel = c2)
        out = - _np.log(aod1/aod2) / _np.log(c1ex/c2ex)

        # out = _timeseries.TimeSeries(_pd.DataFrame(out))
        setattr(self, 'ang_exp_{}_{}'.format(column_1, column_2), out)
        return out

    def cal_deviations(arglist):
        """
        currently not working, fix it if needed

        Parameters
        ----------
        arglist : TYPE
            DESCRIPTION.

        Returns
        -------
        diff : TYPE
            DESCRIPTION.

        """
        # [idx, aod], [idxt, wls] = arglist

        # assert (len(wls) == len(aod))
        # sfr.ang_exp.data.loc[idx, :]
        # fitres = sp.stats.linregress(np.log10(wls.loc[[500, 870]]), np.log10(aod.loc[[500, 870]]))
        # ff = lambda x: 10 ** (np.log10(x) * fitres.slope + fitres.intercept)
        # aod_fit = ff(wls)
        # diff = (aod - aod_fit)
        # return diff

    def calculate_deviation_from_angstrom(sfr, wl1=500, wl2=870):
        """
        currently not working, fix it if needed
        ----------
        wl1
        wl2

        Returns
        -------

        """
        # pool = _mp.Pool(6)

        # out = pool.map(cal_deviations, zip(sfr.AOD.data.iterrows(), sfr.wavelength_matching_table.data.iterrows()))

        # # make dataframe from output
        # df = pd.concat(out, axis=1)
        # df = df.transpose()

        # # xarray does not like timezone ... remove it
        # df.index = pd.to_datetime(df.index.astype(str).str[:-6])

        # # xarray does not like integer column ames
        # df.columns = df.columns.astype(str)

        # # other stuff
        # df.sort_index(inplace=True)
        # df.index.name = 'datetime'
        # return df
    
    @property
    def AODinversion(self):
        if isinstance(self._aodinv, type(None)):
            self.derive_size_distribution()
        return self._aodinv
        
    
    def invertAOD(self, width_of_aerosol_mode = (0.2, 0.25), 
                                 channels = [500,670, 870, 1625],
                                 all_valid = True, 
                                 verbose = False, 
                                 test = False,
                                 cut_off = 1e-10, 
                                 returntype = 'InversionAOD'):
        """
        get the size distribution

        Parameters
        ----------
        width_of_aerosol_mode : TYPE, optional
            DESCRIPTION. The default is 0.15.
        all_valid : TYPE, optional
            All AOD values need to be valid for a retrieval to be derived. Otherwise fit results are invalid. The default is True.
        verbose : TYPE, optional
            DESCRIPTION. The default is True.
        test : TYPE, optional
            DESCRIPTION. The default is False.
        cut_off: float,
            cut_off for fit termination conditions. 1e-10 is recommended. At the default of 1e-8 convergence is not always good.
        returntype: str, ['InversionAOD', 'InversionSizeDistribution']
            What is returned, 
            InversionAOD: new format with focus on AOD partitioning
            InversionSizeDistribution: older output when I was trying to retriev a sizedistribution.
                          
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        assert(all_valid), 'all_valid needs to be true, programming required if you insist on False.'
        # dist_df = _pd.DataFrame()
        # fit_res_df = _pd.DataFrame()
        aoddata = self.aod.to_pandas()
        count_down = aoddata.shape[0] #self.AOD.data.shape[0]
        result_ds = None
        self.tp_test1 = False
        if verbose:
            print(f'Starting to make inversion for {count_down} cases.')
        for e,(ts, test_dp) in enumerate(aoddata.iterrows()): #self.AOD.data.iterrows():
            # if e != 1500:
            #     continue
            
            test_dp = test_dp.loc[channels] #exclude the 940!!! And maybe the 415?
            self.tp_test_dp = test_dp
            self.tp_ts = ts
            
            # generate dataset
            ds = xr.Dataset()
            
            df = _pd.DataFrame(columns= channels, index = ['fine', 'coarse'])
            df.index.name = 'mode'            
            df.columns.name = 'channel'
            ds['aod_fine_coarse'] = df
            
            ds['cost'] = _np.nan
            ds['cpu_usage'] = _np.nan
            
            df = _pd.DataFrame(columns= ['diameter', 'number', 'width'], index = ['fine', 'coarse'])
            df.index.name = 'mode'
            df.columns.name = 'sd_param'
            ds['size_distribution_parameters'] = df
            
            if verbose:
                print(test_dp)
            
            # do not fit if not all aod values are given -> Values remain nan
            if not test_dp.isna().any():
                inv = _atmop.Inversion2SizeDistribution(test_dp, width_of_aerosol_mode = width_of_aerosol_mode, verbose=verbose)
                inv.fit_cutoff = cut_off
                
                # if inv.fit_result.full_result.cost > 1e-7:
                #     print('other starting params')
                #     inv = _atmop.Inversion2SizeDistribution(test_dp, width_of_aerosol_mode = width_of_aerosol_mode, 
                #                                             start_args = [120, 2000.0, 500, 1],
                #                                             verbose=verbose)
                #     inv.fit_cutoff = cut_off
                
                self.tp_inv = inv
                
                
                #### set values
                # aod_fine_coarse
                df = inv.fit_result.size_distribution_aods
                df.index.name = 'mode'
                ds['aod_fine_coarse'] = df
                ds.aod_fine_coarse.attrs['info'] = f'AOD as attributed to the fine and coarse mode by the {_atmop.__name__} module.'
                # stderr
                ds['cost'] = _np.sqrt(inv.fit_result.full_result.cost)
                ds.cost.attrs['info'] = f'Residual cost after optimization. For more information checkout the cost function in the Inversion2SizeDiatribution class of the {_atmop.__name__} module.'
                # ds['optimization_time'] = inv.fit_result.optimization_time / _pd.to_timedelta(1,'s')
                ds['cpu_usage'] = inv.fit_result.cpu_usage
                ds.cpu_usage.attrs['info'] = "Measure of how long the optimization took; cpu_usage in percent times the duration of the computation in seconds."
                # size_distribution_parameters
                df = inv.fit_result.size_distribution_parameters
                df.columns.name = 'sd_param'
                df.index = ['fine', 'coarse'] # this should be passed by the inv class?
                df.index.name = 'mode'
                ds['size_distribution_parameters'] = df
                
            # add the timestamp for concatination later
            ds = ds.assign_coords(datetime = test_dp.name)
            ds = ds.expand_dims(['datetime',])
            # self.tp_test1 = True
            # self.tp_ds = ds.copy()
            if isinstance(result_ds, type(None)):
                result_ds = ds
            else:
                result_ds = xr.concat([result_ds, ds], 'datetime')
            if test:
                if self.tp_test1:
                    break
        self._aodinv = result_ds
        return result_ds
        

class CloudDetection(object):
    def __init__(self, parent):
        self.parent = parent
        self._cloudmask_AOD = None
        self._cloudmask_angstrom = None
        self._cloudmask_combined = None
        self._cloudmask_michalsky = None
        self._cloudmask_augustine = None
        self._masked_data_AOD = None
        self._cloud_classifyiers_AOD = None
        self._masked_data_angstrom = None
        self._cloud_classifyiers_angstrom = None
        self._cloudmask_nativ = None
        self._cloud_classifyiers_combined = None
    
    @property
    def cloudmask_michalsky(self):
        if isinstance(self._cloudmask_michalsky, type(None)):
            self._cloudmask_michalsky = cloud_screening_michalsky(self.parent.dataset)
        return self._cloudmask_michalsky
    
    @property
    def cloudmask_augustine(self):
        if isinstance(self._cloudmask_augustine, type(None)):
            self._cloudmask_augustine = cloud_screening_michalsky_smoker(self.parent.dataset)
        return self._cloudmask_augustine
                
    @property
    def cloudmask_nativ(self):
        return self._cloudmask_nativ
    
    @cloudmask_nativ.setter
    def cloudmask_nativ(self, value):
        self._cloudmask_nativ = value
        return
    
    @property
    def cloudmask_AOD(self):
        if isinstance(self._cloudmask_AOD, type(None)):
            self.get_custom_cloudmask_AOD()
        return self._cloudmask_AOD        
    
    @property
    def cloudmask_angstrom(self):
        if isinstance(self._cloudmask_angstrom, type(None)):
            self.get_custom_cloudmask_angstrom()
        return self._cloudmask_angstrom    
    
    @property
    def cloudmask_combined(self):
        if isinstance(self._cloudmask_combined, type(None)):
            self.get_custom_cloudmask_combined()
        return self._cloudmask_combined    
    
    @property
    def masked_data_AOD(self):
        if isinstance(self._masked_data_AOD, type(None)):
            self.get_custom_cloudmask_AOD()
        return self._masked_data_AOD   
    
    @property
    def masked_data_angstrom(self):
        if isinstance(self._masked_data_angstrom, type(None)):
            self.get_custom_cloudmask_angstrom()
        return self._masked_data_angstrom  
    
    @property
    def cloud_classifyiers_AOD(self):
        if isinstance(self._cloud_classifyiers_AOD, type(None)):
            self.get_custom_cloudmask_AOD()
        return self._cloud_classifyiers_AOD    
        
        
 
    
    @property
    def cloud_classifyiers_angstrom(self):
        if isinstance(self._cloud_classifyiers_angstrom, type(None)):
            self.get_custom_cloudmask_angstrom()
        return self._cloud_classifyiers_angstrom    

    @property
    def cloud_classifyiers_combined(self):
        if isinstance(self._cloud_classifyiers_combined, type(None)):
            self.get_custom_cloudmask_combined()
        return self._cloud_classifyiers_combined



    def get_custom_cloudmask_combined(self, 
                                      deriv_corr_discriminator = -0.3,
                                      deriv_corr_window = 15,
                                      # mean_discriminator = 0.1,
                                      # mean_window = 3,
                                      # linreg_discriminator = 20e-05,
                                      # deriv_discriminator = 5e-4,
                                      # linreg_window = 4,
                                      direction = 'forward',
                                       min_consec_valid = 15,
                                     ):
        """
        

        Parameters
        ----------
        deriv_corr_discriminator : TYPE, optional
            DESCRIPTION. The default is -2.
        deriv_corr_window : TYPE, optional
            DESCRIPTION. The default is 15.
        # mean_discriminator : TYPE, optional
            DESCRIPTION. The default is 0.1.
        # mean_window : TYPE, optional
            DESCRIPTION. The default is 3.
        # linreg_discriminator : TYPE, optional
            DESCRIPTION. The default is 20e-05.
        # deriv_discriminator : TYPE, optional
            DESCRIPTION. The default is 5e-4.
        # linreg_window : TYPE, optional
            DESCRIPTION. The default is 4.
        direction : str, optional
            Which direction to consider (forward, backward, both). The default is 'forward'.
        min_consec_valid : TYPE, optional
            DESCRIPTION. The default is 15.
         : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        

        settings = _pd.DataFrame({'deriv_corr':{'discriminator': deriv_corr_discriminator, 'window': deriv_corr_window, 
                                                'min_consec_valid': min_consec_valid,
                                                },
                                  # 'linreg':{'discriminator': linreg_discriminator, 'window': linreg_window, 'min_consec_valid': min_consec_valid},
                                  # 'deriv':{'discriminator': deriv_discriminator, 
                                  #          # 'window': linreg_window,
                                  #          'min_consec_valid': min_consec_valid}
                                  })
        self._setting_combined = settings
        # settings

        # cloudmask = pd.DataFrame()
        # cloud_classifyiers = pd.DataFrame()
        # masked_data = pd.DataFrame()

        data_ncs = self.parent.ang_exp.data.iloc[:,0]
#         channels = data.columns.values.astype(int)

        # masked_data = _xr.DataArray(coords= {#'channel':channels,
        #                                     'datetime': data_ncs.index.values,
        #                                     'criteria': ['deriv_corr']}, dims = ['datetime', 'criteria'])
        cloudmask = _xr.DataArray(coords= {#'channel':channels,
                                            'datetime': data_ncs.index.values,
                                            'criteria': ['deriv_corr']}, dims = ['datetime', 'criteria'])

        cloud_classifyiers = _xr.DataArray(coords= {#'channel':channels,
                                                   'datetime': data_ncs.index.values,
                                                   'criteria': ['deriv_corr'],
                                                   'direction': ['forward', 'backward']}, dims = ['datetime', 'criteria', 'direction'])

#         for ch, data_ncs in aod_ncs.data.iteritems():
        #     break

        for cn, ser in settings.iteritems():
            # get the values to be juged
            df_ang = self.cloud_classifyiers_angstrom.sel(criteria = 'deriv', direction = 'forward').to_pandas()
            df_aod = self.cloud_classifyiers_AOD.sel(criteria = 'deriv', direction = 'forward', channel = 500).to_pandas()
            roll = df_ang.rolling(int(ser.window))
            rstd = roll.corr(df_aod)
            
            df_ang = self.cloud_classifyiers_angstrom.sel(criteria = 'deriv', direction = 'backward').to_pandas()
            df_aod = self.cloud_classifyiers_AOD.sel(criteria = 'deriv', direction = 'backward', channel = 500).to_pandas()
            roll = df_ang.rolling(int(ser.window))
            rstd_r = roll.corr(df_aod)
            
            # rstd = self.get_classifyer(data_ncs, ctype = cn, window = ser.window)['classifyer']
            cloud_classifyiers.loc[dict(criteria = cn, direction = 'forward')] = rstd

            # rstd_r = self.get_classifyer(data_ncs[::-1], ctype = cn, window = ser.window)['classifyer']
            # rstd_r = rstd_r[::-1]
            cloud_classifyiers.loc[dict(criteria = cn, direction = 'backward')] = rstd_r

            if direction in ['both', 'backward']:
                assert(False), 'direction==both and backward is not working right now'
                cloud_mask_t = _np.logical_and(rstd < ser.discriminator, rstd_r < ser.discriminator) 
            elif direction == 'forward':
                cloud_mask_t = rstd < ser.discriminator
            else:
                assert(False), f'{direction} not an option for direction kwarg'
            
            cloud_mask_t = self.check_min_consec_valid(cloud_mask_t,  min_consec_valid=ser.min_consec_valid)
            cloudmask.loc[dict(criteria = cn)] = cloud_mask_t

            # data_cs_new_j = data_ncs.copy()
            # data_cs_new_j[cloud_mask_t == 1] = _np.nan
            # masked_data.loc[dict(criteria = cn)] = data_cs_new_j
                
        self._cloudmask_combined = cloudmask
        # self._masked_data_angstrom = masked_data
        self._cloud_classifyiers_combined = cloud_classifyiers
        
    def get_custom_cloudmask_angstrom(self, 
                                      mean_discriminator = 0.09,
                                      mean_window = 15,
                                      linreg_discriminator = 5e-5,
                                      linreg_window = 15,
                                      deriv_discriminator = [2e-3, 0.01, 35],
                                      min_consec_valid = 15,
                                     ):
        
        
        deriv_discriminator = deriv_discriminator[0] + (deriv_discriminator[1] * _np.exp(-( deriv_discriminator[2] * self.parent.AOD.data[500])))
        self.tp_deriv_discriminator = deriv_discriminator
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
            rstd = self.get_classifyer(data_ncs, ctype = cn, window = ser.window)['classifyer']
            cloud_classifyiers.loc[dict(criteria = cn, direction = 'forward')] = rstd

            rstd_r = self.get_classifyer(data_ncs[::-1], ctype = cn, window = ser.window)['classifyer']
            rstd_r = rstd_r[::-1]
            cloud_classifyiers.loc[dict(criteria = cn, direction = 'backward')] = rstd_r

            cloud_mask_t = _np.logical_or(rstd.abs() > ser.discriminator, rstd_r.abs() > ser.discriminator) 
            cloud_mask_t[_np.logical_or(_np.isnan(rstd), _np.isnan(rstd_r))] = True
#             self.tp_mask = cloud_mask_t.copy()
#             self.tp_disc = ser.discriminator
#             self.tp_rstd = rstd.copy()
#             self.tp_rstd_r = rstd_r.copy()
#             self.tp_cn = cn
#             self.tp_data_ncs = data_ncs

            if 1: # lets switch this of for now
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
        return

    def get_custom_cloudmask_AOD(self, 
                                      mean_discriminator = 6e-2,
                                      mean_window = 15,
                                      linreg_discriminator = 2.5e-5,
                                      linreg_window = 15,
                                      deriv_discriminator = 6e-5,
                                      min_consec_valid = 15,
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
                cloud_mask_t = _np.logical_or(rstd.abs() > ser.discriminator, rstd_r.abs() > ser.discriminator)
                cloud_mask_t[_np.logical_or(_np.isnan(rstd), _np.isnan(rstd_r))] = True
                
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
        return 

    # def plot_classifyer_angstrom(self, criteria = 'linreg', direction = 'forward', ax = None):
    #     if criteria == 'all':
    #         criteria = self.cloud_classifyiers_AOD.criteria.values
    #     else:
    #         criteria = [criteria]
        
    #     if isinstance(ax, type(None)):
    #         f,aa = _plt.subplots(len(criteria), sharex = True, gridspec_kw={'hspace':0})
    #     else:
    #         aa = ax
    #         f = aa[0].get_figure()
    #     # df = self.cloud_classifyiers_AOD.sel(criteria = criteria, direction = direction).to_pandas().transpose()

            
    #     for e,crit in enumerate(criteria):
        
        
        
        
    #     if isinstance(ax, type(None)):
    #         a = _plt.subplot()
    #     else:
    #         a = ax
    #     df = self.cloud_classifyiers_angstrom.sel(criteria = criteria, direction = direction).to_pandas()#.transpose()
        
    #     df.plot(ax = a)
    #     for g in a.get_lines():
    #         g.set_linestyle('')
    #         g.set_marker('.')
    #     a.set_title(f'classifyer - criteria: {criteria}, direction: {direction}')
    #     a.set_xlabel('')
    #     a.axhline(self._setting_angstrom.loc['discriminator', criteria], color = 'black', ls = '--')
    #     a.set_yscale('log')
    #     return a
    
    def plot_classifyer(self, classifyer = 'AOD', criteria = 'linreg', direction = 'forward', channel = 500, ax = None):
        """
        

        Parameters
        ----------
        classifyer: str, optional
            AOD, angstrom, combined
        criteria : str, optional
            'linreg', 'mean', 'deriv', or 'all'. The default is 'linreg'.
        direction : TYPE, optional
            DESCRIPTION. The default is 'forward'.
        channel : TYPE, optional
            DESCRIPTION. The default is None.
        ax : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        if classifyer == 'AOD':
            parameter = self.cloud_classifyiers_AOD
            settings = self._setting_AOD
        elif classifyer == 'angstrom':
            parameter = self.cloud_classifyiers_angstrom
            settings = self._setting_angstrom
        elif classifyer == 'combined':
            parameter = self.cloud_classifyiers_combined
            settings = self._setting_combined
        else:
            assert(False), 'nonononon, what did you do????'
        
        if criteria == 'all':
            criteria = parameter.criteria.values
        else:
            criteria = [criteria]
        
        if isinstance(ax, type(None)):
            f,aa = _plt.subplots(len(criteria), sharex = True, gridspec_kw={'hspace':0})
            if len(criteria) == 1:
                aa = [aa]
        else:
            aa = ax
            f = aa[0].get_figure()
        # df = self.cloud_classifyiers_AOD.sel(criteria = criteria, direction = direction).to_pandas().transpose()

            
        for e,crit in enumerate(criteria):
            
            df = parameter.sel(criteria = crit, direction = direction)
            if not isinstance(channel, type(None)) and classifyer == 'AOD':
                 df = df.sel(channel = channel)
            df = df.to_pandas().transpose()
    
            a = aa[e]
            df.plot(ax = a, color = _colors[e])
            for g in a.get_lines():
                g.set_linestyle('')
                g.set_marker('.')
            # a.set_title(f'classifyer - criteria: {criteria}, direction: {direction}')
            a.set_ylabel(crit)
            a.set_xlabel('')
            a.axhline(settings.loc['discriminator', crit], color = 'black', ls = '--')
            if classifyer in ['AOD', 'angstrom']:
                a.set_yscale('log')
            a.legend(fontsize = 'small', title = 'channel (nm)').remove()
        return f,aa
    
    # def plot_classifyer_combined(self, direction = 'forward', ax = None):
    #     criteria = 'deriv_corr'
    #     if isinstance(ax, type(None)):
    #         a = _plt.subplot()
    #     else:
    #         a = ax
    #     df = self.cloud_classifyiers_combined.sel(criteria = criteria, direction = direction).to_pandas()#.transpose()
        
    #     df.plot(ax = a)
    #     for g in a.get_lines():
    #         g.set_linestyle('')
    #         g.set_marker('.')
    #     a.set_title(f'classifyer - criteria: {criteria}, direction: {direction}')
    #     a.set_xlabel('')
    #     a.axhline(self._setting_combined.loc['discriminator', criteria], color = 'black', ls = '--')
    #     # a.set_yscale('log')
    #     return a
    
    def plot_cloudmask(self, classifyer = 'combined', ax = None):
        if classifyer == 'combined':
            df = self.cloudmask_combined.to_pandas()
        elif classifyer == 'AOD':
            df = self.cloudmask_AOD.sel(channel = 500).to_pandas()
        elif classifyer == 'angstrom':
            df = self.cloudmask_angstrom.to_pandas()
        # elif classifyer == 'nativ':
        #     df = self.cloudmask_nativ
        
        if isinstance(ax, type(None)):
            a = _plt.subplot()
        else:
            a = ax
    
        offset = -0.1
        # f,a = plt.subplots()
    
        for e, col in enumerate(df):
            (df[col] + (e*offset)).plot(ax = a)
            g = a.get_lines()[-1]
            g.set_label(col)
        
        if classifyer =='AOD':
            e+=1
            (self.cloudmask_nativ.data + (e*offset)).plot(ax = a)
            g = a.get_lines()[-1]
            g.set_label('nativ')
    
        for g in a.get_lines():
            g.set_linestyle('')
            g.set_marker('.')
    
        yrange = offset * e
    
        # a.set_ylim(0 - yrange * 0.1, yrange * (1 + 0.1))
        a.set_ylim(yrange * (1 + 0.1), 0 - yrange * 0.1)
        a.set_yticklabels([])
        a.set_ylabel(classifyer)
        leg = a.legend(loc = (1.05, 0.05)).remove()
        a.set_xlabel('')
        return a
    
    def plot_masked_data(self, data2plot = 'AOD', classifyer = False, 
                         data2discriminate = 'AOD', channel_plot = 500, 
                         channel_discriminate = None, criteria = 'mean', 
                         show_unmasked = True, 
                         invert_mask = True,
                         ax = None, marker_size_scale = 2):
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
            
        if invert_mask:
            cloudmask = ~cloudmask
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
        a.legend().remove()
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
                rstd = stderr#/rm
            elif ctype == 'mean':
                rstd = roll.std()#/rm
            
        elif ctype == 'deriv':
            # data = mfrsr.AOD.data#.loc[:, 500]
            rstd = data.diff().div(data.index.to_series().diff().dt.total_seconds(), axis = 'index')#.abs()
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
    
def cloud_screening_michalsky(dataset, testdev = 0.05,rollwindow = 15, channel = 500):
    """
    This is a translation of Joe Michalskies R code. The translation is not 
    excact, so some small discrapencies are expected.

    Parameters
    ----------
    dataset : xarray.Dataset
        A dataset that has a variable named aod with a koordinate named channel 
        that has a value of the keyword "channel".
    testdev : float, optional
        Threshold for initial test. The default is 0.05.
    rollwindow : int, optional
        This the window size for the rolling functions. Note this value gives 
        to the number of minutes that need to be inbetween two cloud passes in 
        order to see the clouds separately. The default is 15.
    channel : int, optional
        This is the name of the wavelength channel at which the cloud screening 
        is perfomed. The default is 500.

    Returns
    -------
    cloudmask : TYPE
        DESCRIPTION.

    """
    # select the channel to perform this on
    daaod5 = dataset.aod.sel(channel = channel)
    saod5 = daaod5.to_pandas()
    
    # This first test is merely a rough one to create the lowess smoothening of aod only 
    # In a rolling window test if any difference between one value and its neighbor exceeds testdev. If not, !!all!! timestamps are excepted (not just the single timestamp that corresponds with the rolling window!! Effectively this will prevent breaks between clouds that are longer than the window to be excepted as clear. This also measn that you can get very close to a cloud!!
    roll = saod5.rolling(f'{rollwindow} min',                         )
    clear_t1 = _pd.DatetimeIndex([])
    for r in roll:    
        clear = r.diff().abs().dropna() < testdev  
        if not _np.all(clear):
            continue
        else:
            clear_t1 = clear_t1.append(clear.index).drop_duplicates()


    # make a timeseries of 0 with the length of the original data
    # Set all timestamps that are excepted as clear to 1
    test1 = saod5.copy()
    test1[:] = 0
    test1.loc[clear_t1] = 1

    # generate the smooth aod only function
    daaod5_t1 = daaod5.where(test1.values).dropna('datetime')
    out = smlowess.lowess(daaod5_t1.values, daaod5_t1.datetime.values,
                 frac = 2/3)
    lowess = _xr.DataArray(out[:,1], coords={'datetime':daaod5_t1.datetime})
    
    # make a DataFrame with all timestamps over which we can roll again using the actual cloud testing
    # this requires the interpolation of the lowess data
    df = _pd.DataFrame()
    df['lowess'] =  lowess
    df.index = daaod5_t1.datetime
    df['aod'] = _pd.DataFrame({'aod':daaod5.to_pandas()})
    df = _pd.DataFrame({'aod':daaod5.to_pandas()})
    dflow = _pd.DataFrame({'lowess':lowess.to_pandas()}).reindex(df.index)
    dflow = dflow.interpolate()
    df['lowess'] = dflow

    # Second round of thesting. This is more suffisticated and applies a threshold that depends on the actual AOD value
    rollt2 = df.rolling(f'{rollwindow} min', 
                        # center = True,
                       )

    clear_t2 = _pd.DatetimeIndex([])

    for r in rollt2:  
        testdev = r.lowess.mean()

        if testdev > 0.10:
            testdev *= 0.1
        elif (testdev < 0.03) and (testdev > 0.025):
            testdev *= 0.3
        elif (testdev <= 0.025) and (testdev > 0.014):
            testdev *= 0.6
        elif testdev <= 0.014:
            testdev *= 0.8
        else:
            testdev *= 0.2

        test = r.aod.diff().abs().max()
        if test < testdev:
            clear_t2 = clear_t2.append(r.aod.dropna().index).drop_duplicates()
        else:
            continue


    test2 = saod5.copy()
    test2[:] = 0
    test2.loc[clear_t2] = 1

    testall = _pd.DataFrame(index = dataset.datetime.to_pandas())
    testall['cloudy'] = 0
    testall.loc[test2[test2 == 0].index] = 1

    # turn into xarray Dataarray
    testall.index.name = 'datetime'
    cloudmask = testall.cloudy.to_xarray()
    cloudmask = cloudmask.where(~daaod5.isnull())
    return cloudmask

def cloud_screening_michalsky_smoker(dataset, testdev = 0.05,rollwindow = 15, channel = '500_870'):
    """
    This is an adaptation of the Michalsky cloudmask to be applied on the angstrom exponent to detect smoke. This was proposed by John Augustine.
    Parameters
    ----------
    dataset : xarray.Dataset
        A dataset that has a variable named aod with a koordinate named channel 
        that has a value of the keyword "channel".
    testdev : float, optional
        Threshold for initial test. The default is 0.05.
    rollwindow : int, optional
        This the window size for the rolling functions. Note this value gives 
        to the number of minutes that need to be inbetween two cloud passes in 
        order to see the clouds separately. The default is 15.
    channel : int, optional
        This is the name of the wavelength channel at which the cloud screening 
        is perfomed. The default is 500.

    Returns
    -------
    cloudmask : TYPE
        DESCRIPTION.

    """
    # select the channel to perform this on
    daaod5 = dataset.angstrom_exp.sel(ang_channels = channel)
    saod5 = daaod5.to_pandas()
    
    # This first test is merely a rough one to create the lowess smoothening of aod only 
    # In a rolling window test if any difference between one value and its neighbor exceeds testdev. If not, !!all!! timestamps are excepted (not just the single timestamp that corresponds with the rolling window!! Effectively this will prevent breaks between clouds that are longer than the window to be excepted as clear. This also measn that you can get very close to a cloud!!
    roll = saod5.rolling(f'{rollwindow} min',                         )
    clear_t1 = _pd.DatetimeIndex([])
    for r in roll:    
        clear = r.diff().abs().dropna() < testdev  
        if not _np.all(clear):
            continue
        else:
            clear_t1 = clear_t1.append(clear.index).drop_duplicates()


    # make a timeseries of 0 with the length of the original data
    # Set all timestamps that are excepted as clear to 1
    test1 = saod5.copy()
    test1[:] = 0
    test1.loc[clear_t1] = 1

    # generate the smooth aod only function
    daaod5_t1 = daaod5.where(test1.values).dropna('datetime')
    out = smlowess.lowess(daaod5_t1.values, daaod5_t1.datetime.values,
                 frac = 2/3)
    lowess = _xr.DataArray(out[:,1], coords={'datetime':daaod5_t1.datetime})
    
    # make a DataFrame with all timestamps over which we can roll again using the actual cloud testing
    # this requires the interpolation of the lowess data
    df = _pd.DataFrame()
    df['lowess'] =  lowess
    df.index = daaod5_t1.datetime
    df['aod'] = _pd.DataFrame({'aod':daaod5.to_pandas()})
    df = _pd.DataFrame({'aod':daaod5.to_pandas()})
    dflow = _pd.DataFrame({'lowess':lowess.to_pandas()}).reindex(df.index)
    dflow = dflow.interpolate()
    df['lowess'] = dflow

    # Second round of thesting. This is more suffisticated and applies a threshold that depends on the actual AOD value
    rollt2 = df.rolling(f'{rollwindow} min', 
                        # center = True,
                       )

    clear_t2 = _pd.DatetimeIndex([])

    for r in rollt2:  
        testdev = r.lowess.mean() * 0.015
        test = r.aod.diff().abs().max()
        if test < testdev:
            clear_t2 = clear_t2.append(r.aod.dropna().index).drop_duplicates()
        else:
            continue


    test2 = saod5.copy()
    test2[:] = 0
    test2.loc[clear_t2] = 1

    testall = _pd.DataFrame(index = dataset.datetime.to_pandas())
    testall['cloudy'] = 0
    testall.loc[test2[test2 == 0].index] = 1

    # turn into xarray Dataarray
    testall.index.name = 'datetime'
    cloudmask = testall.cloudy.to_xarray()
    cloudmask = cloudmask.where(~daaod5.isnull())
    return cloudmask

class AerosolAndCloudDetection(object):
    def __init__(self, parent):
        self.parent = parent
        # self._cloudmask_AOD = None
        self._badmask_angstrom = None
        # self._cloudmask_combined = None
        # self._masked_data_AOD = None
        # self._cloud_classifyiers_AOD = None
        self._masked_data_angstrom = None
        self._cloud_classifyiers_angstrom = None
        self._cloudmask_nativ = None
        self._minutes2bad_angstrom = None
        self._badmask_angstrom_padded = None
        # self._cloud_classifyiers_combined = None
    
    @property
    def cloudmask_nativ(self):
        return self._cloudmask_nativ
    
    @cloudmask_nativ.setter
    def cloudmask_nativ(self, value):
        self._cloudmask_nativ = value
        return
    
    # @property
    # def cloudmask_AOD(self):
    #     if isinstance(self._cloudmask_AOD, type(None)):
    #         self.get_custom_cloudmask_AOD()
    #     return self._cloudmask_AOD        
    
    @property
    def badmask_angstrom(self):
        if isinstance(self._badmask_angstrom, type(None)):
            self.get_custom_badmask_angstrom()
        return self._badmask_angstrom  
    
    @property
    def badmask_angstrom_padded(self):
        if isinstance(self._badmask_angstrom_padded, type(None)):
            self.get_custom_badmask_angstrom()
        return self._badmask_angstrom_padded 
    
    @property
    def minutes2bad_angstrom(self):
        if isinstance(self._minutes2bad_angstrom, type(None)):
            self.get_custom_badmask_angstrom()
        return self._minutes2bad_angstrom    
    
    # @property
    # def cloudmask_combined(self):
    #     if isinstance(self._cloudmask_combined, type(None)):
    #         self.get_custom_cloudmask_combined()
    #     return self._cloudmask_combined    
    
    # @property
    # def masked_data_AOD(self):
    #     if isinstance(self._masked_data_AOD, type(None)):
    #         self.get_custom_cloudmask_AOD()
    #     return self._masked_data_AOD   
    
    @property
    def masked_data_angstrom(self):
        if isinstance(self._masked_data_angstrom, type(None)):
            self.get_custom_badmask_angstrom()
        return self._masked_data_angstrom  
    
    # @property
    # def cloud_classifyiers_AOD(self):
    #     if isinstance(self._cloud_classifyiers_AOD, type(None)):
    #         self.get_custom_cloudmask_AOD()
    #     return self._cloud_classifyiers_AOD    
        
        
 
    
    @property
    def bad_classifyiers_angstrom(self):
        if isinstance(self._cloud_classifyiers_angstrom, type(None)):
            self.get_custom_badmask_angstrom()
        return self._cloud_classifyiers_angstrom    

    # @property
    # def cloud_classifyiers_combined(self):
    #     if isinstance(self._cloud_classifyiers_combined, type(None)):
    #         self.get_custom_cloudmask_combined()
    #     return self._cloud_classifyiers_combined



#     def get_custom_cloudmask_combined(self, 
#                                       deriv_corr_discriminator = -0.3,
#                                       deriv_corr_window = 15,
#                                       # mean_discriminator = 0.1,
#                                       # mean_window = 3,
#                                       # linreg_discriminator = 20e-05,
#                                       # deriv_discriminator = 5e-4,
#                                       # linreg_window = 4,
#                                       direction = 'forward',
#                                        min_consec_valid = 15,
#                                      ):
#         """
        

#         Parameters
#         ----------
#         deriv_corr_discriminator : TYPE, optional
#             DESCRIPTION. The default is -2.
#         deriv_corr_window : TYPE, optional
#             DESCRIPTION. The default is 15.
#         # mean_discriminator : TYPE, optional
#             DESCRIPTION. The default is 0.1.
#         # mean_window : TYPE, optional
#             DESCRIPTION. The default is 3.
#         # linreg_discriminator : TYPE, optional
#             DESCRIPTION. The default is 20e-05.
#         # deriv_discriminator : TYPE, optional
#             DESCRIPTION. The default is 5e-4.
#         # linreg_window : TYPE, optional
#             DESCRIPTION. The default is 4.
#         direction : str, optional
#             Which direction to consider (forward, backward, both). The default is 'forward'.
#         min_consec_valid : TYPE, optional
#             DESCRIPTION. The default is 15.
#          : TYPE
#             DESCRIPTION.

#         Returns
#         -------
#         None.

#         """
        

#         settings = _pd.DataFrame({'deriv_corr':{'discriminator': deriv_corr_discriminator, 'window': deriv_corr_window, 
#                                                 'min_consec_valid': min_consec_valid,
#                                                 },
#                                   # 'linreg':{'discriminator': linreg_discriminator, 'window': linreg_window, 'min_consec_valid': min_consec_valid},
#                                   # 'deriv':{'discriminator': deriv_discriminator, 
#                                   #          # 'window': linreg_window,
#                                   #          'min_consec_valid': min_consec_valid}
#                                   })
#         self._setting_combined = settings
#         # settings

#         # cloudmask = pd.DataFrame()
#         # cloud_classifyiers = pd.DataFrame()
#         # masked_data = pd.DataFrame()

#         data_ncs = self.parent.ang_exp.data.iloc[:,0]
# #         channels = data.columns.values.astype(int)

#         # masked_data = _xr.DataArray(coords= {#'channel':channels,
#         #                                     'datetime': data_ncs.index.values,
#         #                                     'criteria': ['deriv_corr']}, dims = ['datetime', 'criteria'])
#         cloudmask = _xr.DataArray(coords= {#'channel':channels,
#                                             'datetime': data_ncs.index.values,
#                                             'criteria': ['deriv_corr']}, dims = ['datetime', 'criteria'])

#         cloud_classifyiers = _xr.DataArray(coords= {#'channel':channels,
#                                                    'datetime': data_ncs.index.values,
#                                                    'criteria': ['deriv_corr'],
#                                                    'direction': ['forward', 'backward']}, dims = ['datetime', 'criteria', 'direction'])

# #         for ch, data_ncs in aod_ncs.data.iteritems():
#         #     break

#         for cn, ser in settings.iteritems():
#             # get the values to be juged
#             df_ang = self.cloud_classifyiers_angstrom.sel(criteria = 'deriv', direction = 'forward').to_pandas()
#             df_aod = self.cloud_classifyiers_AOD.sel(criteria = 'deriv', direction = 'forward', channel = 500).to_pandas()
#             roll = df_ang.rolling(int(ser.window))
#             rstd = roll.corr(df_aod)
            
#             df_ang = self.cloud_classifyiers_angstrom.sel(criteria = 'deriv', direction = 'backward').to_pandas()
#             df_aod = self.cloud_classifyiers_AOD.sel(criteria = 'deriv', direction = 'backward', channel = 500).to_pandas()
#             roll = df_ang.rolling(int(ser.window))
#             rstd_r = roll.corr(df_aod)
            
#             # rstd = self.get_classifyer(data_ncs, ctype = cn, window = ser.window)['classifyer']
#             cloud_classifyiers.loc[dict(criteria = cn, direction = 'forward')] = rstd

#             # rstd_r = self.get_classifyer(data_ncs[::-1], ctype = cn, window = ser.window)['classifyer']
#             # rstd_r = rstd_r[::-1]
#             cloud_classifyiers.loc[dict(criteria = cn, direction = 'backward')] = rstd_r

#             if direction in ['both', 'backward']:
#                 assert(False), 'direction==both and backward is not working right now'
#                 cloud_mask_t = _np.logical_and(rstd < ser.discriminator, rstd_r < ser.discriminator) 
#             elif direction == 'forward':
#                 cloud_mask_t = rstd < ser.discriminator
#             else:
#                 assert(False), f'{direction} not an option for direction kwarg'
            
#             cloud_mask_t = self.check_min_consec_valid(cloud_mask_t,  min_consec_valid=ser.min_consec_valid)
#             cloudmask.loc[dict(criteria = cn)] = cloud_mask_t

#             # data_cs_new_j = data_ncs.copy()
#             # data_cs_new_j[cloud_mask_t == 1] = _np.nan
#             # masked_data.loc[dict(criteria = cn)] = data_cs_new_j
                
#         self._cloudmask_combined = cloudmask
#         # self._masked_data_angstrom = masked_data
#         self._cloud_classifyiers_combined = cloud_classifyiers
        
    def get_custom_badmask_angstrom(self, 
                                      mean_discriminator = None,#0.09,
                                      mean_window = 15,
                                      linreg_discriminator = None, #5e-5,
                                      linreg_window = 15,
                                      deriv_discriminator = [2e-3, 0.01, 35],
                                      distance2bad = 5,
                                        min_consec_valid = None,
                                     ):
        
        
        deriv_discriminator = deriv_discriminator[0] + (deriv_discriminator[1] * _np.exp(-( deriv_discriminator[2] * self.parent.AOD.data[500])))
        self.tp_deriv_discriminator = deriv_discriminator
        
        settings = {}
        if not isinstance(mean_discriminator, type(None)):
            settings['mean'] = {'discriminator': mean_discriminator, 'window': mean_window, 'min_consec_valid': min_consec_valid}
        if not isinstance(linreg_discriminator, type(None)):
            settings['linreg'] = {'discriminator': linreg_discriminator, 'window': linreg_window, 'min_consec_valid': min_consec_valid}
        if not isinstance(deriv_discriminator, type(None)):
            settings['deriv'] ={'discriminator': deriv_discriminator, 
                                 # 'window': linreg_window,
                                 # 'min_consec_valid': min_consec_valid
                                 }
        # settings = _pd.DataFrame({'mean':{'discriminator': mean_discriminator, 'window': mean_window, 'min_consec_valid': min_consec_valid},
        #                           'linreg':{'discriminator': linreg_discriminator, 'window': linreg_window, 'min_consec_valid': min_consec_valid},
        #                           'deriv':{'discriminator': deriv_discriminator, 
        #                                    # 'window': linreg_window,
        #                                    'min_consec_valid': min_consec_valid
        #                                    }
        #                           })
        settings = _pd.DataFrame(settings)
        
        self._setting_angstrom = settings
        # settings

        # cloudmask = pd.DataFrame()
        # cloud_classifyiers = pd.DataFrame()
        # masked_data = pd.DataFrame()

        data_ncs = self.parent.ang_exp.data.iloc[:,0]
#         channels = data.columns.values.astype(int)
        # criteria =  ['mean', 'linreg', 'deriv']
        criteria = settings.columns
        masked_data = _xr.DataArray(coords= {#'channel':channels,
                                            'datetime': data_ncs.index.values,
                                            'criteria':criteria}, dims = ['datetime', 'criteria'])
        cloudmask = _xr.DataArray(coords= {#'channel':channels,
                                            'datetime': data_ncs.index.values,
                                            'criteria': criteria}, dims = ['datetime', 'criteria'])
        badmask_padded = _xr.DataArray(coords= {#'channel':channels,
                                            'datetime': data_ncs.index.values,
                                            'criteria': criteria}, dims = ['datetime', 'criteria'])
        
        min2bad = _xr.DataArray(coords= {#'channel':channels,
                                            'datetime': data_ncs.index.values,
                                            'criteria': criteria}, dims = ['datetime', 'criteria'])

        cloud_classifyiers = _xr.DataArray(coords= {#'channel':channels,
                                                   'datetime': data_ncs.index.values,
                                                   'criteria': criteria,
                                                   'direction': ['forward', 'backward']}, dims = ['datetime', 'criteria', 'direction'])

#         for ch, data_ncs in aod_ncs.data.iteritems():
        #     break

        for cn, ser in settings.iteritems():
            # get the values to be juged
            if hasattr(ser, 'window'):
                window = ser.window
            else:
                window = None
            rstd = self.get_classifyer(data_ncs, ctype = cn, window = window)['classifyer']
            cloud_classifyiers.loc[dict(criteria = cn, direction = 'forward')] = rstd

            rstd_r = self.get_classifyer(data_ncs[::-1], ctype = cn, window = window)['classifyer']
            rstd_r = rstd_r[::-1]
            cloud_classifyiers.loc[dict(criteria = cn, direction = 'backward')] = rstd_r

            cloud_mask_t = _np.logical_or(rstd.abs() > ser.discriminator, rstd_r.abs() > ser.discriminator) 

            if 0: # lets switch this of for now
                cloud_mask_t = self.check_min_consec_valid(cloud_mask_t,  min_consec_valid=ser.min_consec_valid)
            cloudmask.loc[dict(criteria = cn)] = cloud_mask_t
            
            #### minutes to next bad point
            min2bad_df = self.get_minutes_to_bad(cloud_mask_t)
            min2bad.loc[dict(criteria = cn)] = min2bad_df
            
            #### padded bad mask
            badmask_ext = min2bad_df.copy()
            where = badmask_ext < distance2bad
            badmask_ext[where] = 1
            badmask_ext[~where] = 0
            badmask_padded.loc[dict(criteria = cn)] = badmask_ext
            
            
            data_cs_new_j = data_ncs.copy()
            data_cs_new_j[cloud_mask_t == 1] = _np.nan
            masked_data.loc[dict(criteria = cn)] = data_cs_new_j
            # minimum number of points between cloud detection?
#             data_cs_new_j_mcv = check_min_consec_valid(data_cs_new_j, min_consec_valid=ser.min_consec_valid)
        #     masked_data[cn] = data_cs_new_j_mcv

        #     cloudmask[cn] = cloud_mask_t
                
        self._badmask_angstrom = cloudmask
        self._badmask_angstrom_padded = badmask_padded
        self._masked_data_angstrom = masked_data
        self._cloud_classifyiers_angstrom = cloud_classifyiers
        self._minutes2bad_angstrom = min2bad
        return

  
    
    def plot_classifyer(self, classifyer = 'AOD', criteria = 'linreg', direction = 'forward', channel = 500, ax = None):
        """
        

        Parameters
        ----------
        classifyer: str, optional
            AOD, angstrom, combined
        criteria : str, optional
            'linreg', 'mean', 'deriv', or 'all'. The default is 'linreg'.
        direction : TYPE, optional
            DESCRIPTION. The default is 'forward'.
        channel : TYPE, optional
            DESCRIPTION. The default is None.
        ax : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        if classifyer == 'AOD':
            parameter = self.cloud_classifyiers_AOD
            settings = self._setting_AOD
        elif classifyer == 'angstrom':
            parameter = self.cloud_classifyiers_angstrom
            settings = self._setting_angstrom
        elif classifyer == 'combined':
            parameter = self.cloud_classifyiers_combined
            settings = self._setting_combined
        else:
            assert(False), 'nonononon, what did you do????'
        
        if criteria == 'all':
            criteria = parameter.criteria.values
        else:
            criteria = [criteria]
        
        if isinstance(ax, type(None)):
            f,aa = _plt.subplots(len(criteria), sharex = True, gridspec_kw={'hspace':0})
            if len(criteria) == 1:
                aa = [aa]
        else:
            aa = ax
            f = aa[0].get_figure()
        # df = self.cloud_classifyiers_AOD.sel(criteria = criteria, direction = direction).to_pandas().transpose()

            
        for e,crit in enumerate(criteria):
            
            df = parameter.sel(criteria = crit, direction = direction)
            if not isinstance(channel, type(None)) and classifyer == 'AOD':
                 df = df.sel(channel = channel)
            df = df.to_pandas().transpose()
    
            a = aa[e]
            df.plot(ax = a, color = _colors[e])
            for g in a.get_lines():
                g.set_linestyle('')
                g.set_marker('.')
            # a.set_title(f'classifyer - criteria: {criteria}, direction: {direction}')
            a.set_ylabel(crit)
            a.set_xlabel('')
            a.axhline(settings.loc['discriminator', crit], color = 'black', ls = '--')
            if classifyer in ['AOD', 'angstrom']:
                a.set_yscale('log')
            a.legend(fontsize = 'small', title = 'channel (nm)').remove()
        return f,aa
    
    # def plot_classifyer_combined(self, direction = 'forward', ax = None):
    #     criteria = 'deriv_corr'
    #     if isinstance(ax, type(None)):
    #         a = _plt.subplot()
    #     else:
    #         a = ax
    #     df = self.cloud_classifyiers_combined.sel(criteria = criteria, direction = direction).to_pandas()#.transpose()
        
    #     df.plot(ax = a)
    #     for g in a.get_lines():
    #         g.set_linestyle('')
    #         g.set_marker('.')
    #     a.set_title(f'classifyer - criteria: {criteria}, direction: {direction}')
    #     a.set_xlabel('')
    #     a.axhline(self._setting_combined.loc['discriminator', criteria], color = 'black', ls = '--')
    #     # a.set_yscale('log')
    #     return a
    
    def plot_cloudmask(self, classifyer = 'combined', ax = None):
        if classifyer == 'combined':
            df = self.cloudmask_combined.to_pandas()
        elif classifyer == 'AOD':
            df = self.cloudmask_AOD.sel(channel = 500).to_pandas()
        elif classifyer == 'angstrom':
            df = self.badmask_angstrom.to_pandas()
        # elif classifyer == 'nativ':
        #     df = self.cloudmask_nativ
        
        if isinstance(ax, type(None)):
            a = _plt.subplot()
        else:
            a = ax
    
        offset = -0.1
        # f,a = plt.subplots()
    
        for e, col in enumerate(df):
            (df[col] + (e*offset)).plot(ax = a)
            g = a.get_lines()[-1]
            g.set_label(col)
        
        if classifyer =='AOD':
            e+=1
            (self.cloudmask_nativ.data + (e*offset)).plot(ax = a)
            g = a.get_lines()[-1]
            g.set_label('nativ')
    
        for g in a.get_lines():
            g.set_linestyle('')
            g.set_marker('.')
    
        yrange = offset * e
    
        # a.set_ylim(0 - yrange * 0.1, yrange * (1 + 0.1))
        a.set_ylim(yrange * (1 + 0.1), 0 - yrange * 0.1)
        a.set_yticklabels([])
        a.set_ylabel(classifyer)
        leg = a.legend(loc = (1.05, 0.05)).remove()
        a.set_xlabel('')
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
            cloudmask = self.badmask_angstrom.sel(criteria = criteria).values == 1

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
        a.legend().remove()
        return a
    
    @staticmethod
    def get_minutes_to_bad(cmt):
        cmtones = cmt[cmt == 1]
        if cmtones.shape[0] == 0:
            dist_df = _pd.Series(index = cmt.index)
            dist_df[:] = 9999
        else:
            dist_df = _pd.DataFrame(index = cmt.index)
            for tstp in cmt[cmt == 1].index:
                dist_df[tstp] =  abs(dist_df.index - tstp)
        
            dist_df = dist_df.min(axis = 1)
        
            dist_df = dist_df / _pd.to_timedelta(1, 'min')
        return dist_df
    
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
            rstd = data.diff().div(data.index.to_series().diff().dt.total_seconds(), axis = 'index')#.abs()
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
    
    
class AOD_AOT_20221216(object):
    """This is an older version that might still be used (subclassed) by the surfrad module"""
    def __init__(self,
                 AOD = None,
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
            If site is given the remaining site relevant kwargs (lat, lon, elevation, name, name_short) are ignored.
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
        # self._aod = None
        self.AOD = AOD
        self._sunposition = None
        self._timezone = timezone
        self.wavelengths = wavelengths

        if not isinstance(site, type(None)):
            self.site = site
        elif not isinstance(lat, type(None)):
            self.site = _measurement_site.Station(lat, lon, elevation, name=name, abbreviation=name_short, info = site_info)


        self.cloudmask = CloudDetection(self)
        self.aerosolmask = AerosolAndCloudDetection(self)

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
        if isinstance(value, type(None)):            
            self._aod = value
            return
        elif isinstance(value, _pd.DataFrame):
            value = _timeseries.TimeSeries(value)
            
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
        """
        currently not working, fix it if needed

        Parameters
        ----------
        arglist : TYPE
            DESCRIPTION.

        Returns
        -------
        diff : TYPE
            DESCRIPTION.

        """
        # [idx, aod], [idxt, wls] = arglist

        # assert (len(wls) == len(aod))
        # sfr.ang_exp.data.loc[idx, :]
        # fitres = sp.stats.linregress(np.log10(wls.loc[[500, 870]]), np.log10(aod.loc[[500, 870]]))
        # ff = lambda x: 10 ** (np.log10(x) * fitres.slope + fitres.intercept)
        # aod_fit = ff(wls)
        # diff = (aod - aod_fit)
        # return diff

    def calculate_deviation_from_angstrom(sfr, wl1=500, wl2=870):
        """
        currently not working, fix it if needed
        ----------
        wl1
        wl2

        Returns
        -------

        """
        # pool = _mp.Pool(6)

        # out = pool.map(cal_deviations, zip(sfr.AOD.data.iterrows(), sfr.wavelength_matching_table.data.iterrows()))

        # # make dataframe from output
        # df = pd.concat(out, axis=1)
        # df = df.transpose()

        # # xarray does not like timezone ... remove it
        # df.index = pd.to_datetime(df.index.astype(str).str[:-6])

        # # xarray does not like integer column ames
        # df.columns = df.columns.astype(str)

        # # other stuff
        # df.sort_index(inplace=True)
        # df.index.name = 'datetime'
        # return df
        
    def derive_size_distribution(self, width_of_aerosol_mode = 0.15, all_valid = True, verbose = True, test = False):
        """
        get the size distribution

        Parameters
        ----------
        width_of_aerosol_mode : TYPE, optional
            DESCRIPTION. The default is 0.15.
        all_valid : TYPE, optional
            All AOD values need to be valid for a retrieval to be derived. Otherwise fit results are invalid. The default is True.
        verbose : TYPE, optional
            DESCRIPTION. The default is True.
        test : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        assert(all_valid), 'all_valid needs to be true, programming required if you insist on False.'
        dist_df = _pd.DataFrame()
        fit_res_df = _pd.DataFrame()
        count_down = self.AOD.data.shape[0]
        if verbose:
            print(f'Starting to make inversion for {count_down} cases.')
        for ts, test_dp in self.AOD.data.iterrows():
            if verbose:
                print(test_dp)
            
            inv = _atmop.Inversion2SizeDistribution(test_dp, width_of_aerosol_mode = width_of_aerosol_mode, verbose=verbose)
            inv.fit_cutoff = 1e-8
        # inv.start_conditions.size_distribution.convert2dVdlogDp().plot()
        # inv.fit_result.size_distribution.convert2dVdlogDp().plot()
        
        # [0.0015510533, 0.0014759477, 0.0012328415, 0.00094773073, 0.00032288354]
        # inv.start_conditions.extinction_coeff
        
            inv.fit_result.size_distribution.data.index = [ts]
            sdt = inv.fit_result.size_distribution.data.copy()
            # dist_df = dist_df.append(sdt)
            dist_df =_pd.concat([dist_df, sdt])
        
            data = _np.array([_np.append(inv.fit_result.args, inv.fit_result.sigma)])
            # fit_res_df = fit_res_df.append(_pd.DataFrame(data, columns=['pos1', 'amp1', 'pos2', 'amp2', 'sigma'], index = [ts]))
            dft = _pd.DataFrame(data, columns=['pos1', 'amp1', 'pos2', 'amp2', 'sigma'], index = [ts])
            fit_res_df = _pd.concat([fit_res_df, dft])
            count_down -=1
            print(count_down, end = ', ')
            
            if test:
                break
        
        
        
        distt = inv.fit_result.size_distribution.distributionType
        bins = inv.fit_result.size_distribution.bins
        # dist_ts = sd.SizeDist_TS(dist_df, bins, distt)
        dist_ts = _atmsd.SizeDist_TS(dist_df, bins, distt, ignore_data_gap_error=True)
        out= {}
        out['dist_ts'] = dist_ts
        out['last_inv_inst'] = inv
        fit_res_df.index.name = 'Time'
        fit_res_df.columns.name = 'fit_params'
        out['fit_res'] = fit_res_df
        return InversionSizeDistribution(self, dist_ts, fit_res_df, width_of_aerosol_mode)