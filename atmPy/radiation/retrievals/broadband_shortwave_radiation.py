"""This is a collection of classes the are uses as the basis of many retrievals."""

import numpy as np
import xarray as xr
import pandas as pd
import atmPy.radiation.retrievals.clearsky as atmcsk
import matplotlib.pyplot as _plt
import warnings



default_config = dict(# General filters
                        mu0_min = 0.05,
                        # NSW magnitude test parameters
                        nsw_exp = 1.2,
                        nsw_min = 800.0,
                        nsw_max = 1400.0,
                        # Diffuse magnitude test parameters
                        diffuse_max_coeff = 600,
                        diffuse_max_exp = 0.5,
                        # Change-with-time test parameters
                        max_dsw_dt = 75.0,  # W m-2 per minute
                        # NDR variability test parameters
                        ndr_exp = -0.8,
                        ndr_std_max = 0.005,
                        ndr_window = 11,
                        # # Output options
                        # return_tests = False,
                     )

class _DatasetRef(object):
    def __init__(self, dataset):
        self.dataset = dataset


class SolarIrradiation(object):
    def __init__(self, dataset, site = None, verbose = False, _dataset_ref = None):
        self.verbose = verbose  
        self._dataset_ref = _dataset_ref if _dataset_ref is not None else _DatasetRef(dataset)
        if _dataset_ref is None:
            self.dataset = dataset
        self.site = site
        self._sun_position_variables = ['zenith', 'zenith_geometric', 'apparent_elevation', 'elevation', 'azimuth', 'equation_of_time', 'airmass', 'airmass_absolute', 'sun_earth_distance']
        assert('datetime' in dataset), 'Time coordinate has to be called datetime .... sorry, i know that is an unconventional choise for the time cooridinate name.'

    @property
    def dataset(self):
        return self._dataset_ref.dataset

    @dataset.setter
    def dataset(self, value):
        self._dataset_ref.dataset = value

    def drop_vars(self, names, *, errors = 'raise'):
        """Drop variables and update the shared dataset reference."""
        self.dataset = self.dataset.drop_vars(names, errors = errors)
        return self.dataset

    def get_attr(self, attr):
        if attr not in self.dataset.attrs:
            if self.verbose:
                print(f'{attr} attribute is not set, using default')
            self.dataset.attrs[attr] = default_config[attr]
        return self.dataset.attrs[attr]

    @property
    def mu0(self):
        if 'mu0' not in self.dataset:
            self.dataset['mu0'] = np.cos(self.sun_position.zenith) 
        return self.dataset['mu0']
    
    @property
    def sun_position(self):
        if not np.all([v in self.dataset for v in self._sun_position_variables]):
        # if isinstance(self._sun_position, type(None)):
            sp = self.site.get_sun_position(self.dataset.datetime)
            for v in self._sun_position_variables:
                self.dataset[v] = sp[v]
        return self.dataset[self._sun_position_variables]
    
    @sun_position.setter
    def sun_position(self, value):
        self.dataset.datetime.identical(value.datetime)
        for v in self._sun_position_variables:
            self.dataset[v] = value[v]



class GlobalHorizontalIrradiation(SolarIrradiation):
    def __init__(self, dataset, **kwargs):
        super().__init__(dataset, **kwargs)
        assert('global_horizontal' in dataset), f'global_horizontal variable is missing.'

    @property
    def mask_variability(self,
                                            # global_irradiance: xr.DataArray,
                                            # time_dim: str,
                                            # max_dsw_dt: float,
                                        ):
        if 'mask_global_irradiance_variability' not in self.dataset:
            if self.verbose:
                print('Running global irradiance variability test')
            mask = atmcsk.global_irradiance_variability_test(self, 
                                                             self.dataset.global_horizontal, 
                                                             max_dsw_dt = self.get_attr('max_dsw_dt'))
            self.dataset['mask_global_irradiance_variability'] = mask
        return self.dataset.mask_global_irradiance_variability
    
    @property
    def mask_normalized_global_magnitude(self):
        if 'mask_normalized_global_magnitude' not in self.dataset:
            if self.verbose:
                print('Running normalized global magnitude test')
            self.dataset['mask_normalized_global_magnitude'] = atmcsk.normalized_global_magnitude_test(self.dataset.global_horizontal,
                                                                                                        self.mu0,
                                                                                                        mu0_min = self.get_attr('mu0_min'), 
                                                                                                        nsw_exp = self.get_attr('nsw_exp'),
                                                                                                        nsw_min = self.get_attr('nsw_min'),
                                                                                                        nsw_max = self.get_attr('nsw_max'),)
        return self.dataset.mask_normalized_global_magnitude


class DiffuseHorizontalIrradiation(SolarIrradiation):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert('diffuse_horizontal' in self.dataset), f'diffuse_horizontal variable is missing.'

    @property
    def mask_magnitude(self):
        if 'mask_diffuse_magnitude' not in self.dataset:
            if self.verbose:
                print('Running diffuse magnitude test')
            self.dataset['mask_diffuse_magnitude'] = atmcsk.diffuse_magnitude_test(self.dataset.diffuse_horizontal,
                                                                                    self.mu0,
                                                                                    mu0_min = self.get_attr('mu0_min'),
                                                                                    diffuse_max_coeff = self.get_attr('diffuse_max_coeff'),
                                                                                    diffuse_max_exp = self.get_attr('diffuse_max_exp'),)
        return self.dataset.mask_diffuse_magnitude
    



class DirectNormalIrradiation(SolarIrradiation):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert(('direct_normal' in self.dataset) or ('direct_horizontal' in self.dataset)), f'direct_normal or alternatively direct_horizontal variable is missing.'

        if 'direct_horizontal' in self.dataset and 'direct_normal' not in self.dataset:
            self.direct_normal_from_direct_horizontal()


    def direct_normal_from_direct_horizontal(self):
        """Converts direct horizontal irradiance to direct normal irradiance. 'direct_normal' will be added to the dataset if it is not already present. """
        if 'direct_normal' not in self.dataset:
            ds = self.dataset
            ds['direct_normal'] = ds.direct_horizontal / xr.DataArray(np.cos(self.sun_position.zenith))
        return ds['direct_normal']
    

class CombinedGlobalDiffuseDirect(SolarIrradiation):
    def __init__(self, dataset, **kwargs):
        """
        A class that combines the three irradiation components into one dataset and provides some specific functions
        Parameters
        ----------
        dataset : xr.Dataset
            A dataset that contains the three components as variables named 'global_horizontal', 'diffuse_horizontal', 'direct_normal'
            and the coordinates 'datetime' and 'channel'
        -----------
        2024-06-10, HTelg   
        """

        super().__init__(dataset, **kwargs)
        self._global_horizontal_irradiation = None
        self._diffuse_horizontal_irradiation = None
        self._direct_normal_irradiation = None
        self._kwargs = kwargs

    @property
    def clearsky_global_irradiation_powerlow_fit_params(self):
        """Performes a empirical power law fit to the clear sky global irradiance
        and returns the fit parameters."""

        if 'clearsky_global_irradiation_powerlow_fit_params' not in self.dataset:
            if self.verbose:
                print('performing powerlow fit for clear sky global irradiance')
            params = atmcsk.fit_global_powerlaw_mu0(self.mu0, self.dataset.global_horizontal, self.mask_clear_sky_radflux, 
                                                        mu0_min = self.get_attr('mu0_min'),
                                                        min_points = 100)
            if isinstance(params, type(None)) and self.verbose:
                print('fit failed, probably not enough clear sky points')
                
            self.dataset['clearsky_global_irradiation_powerlow_fit_params'] = params
        return self.dataset['clearsky_global_irradiation_powerlow_fit_params']
        # return params
        #     def fit_func(mu0_values):
        # return A * np.power(mu0_values, b)

    def plot_clearsky_global_irradiation_powerlow_fit(self, ax = None):
        params = self.clearsky_global_irradiation_powerlow_fit_params
        if params.isnull():
            print('Not enough clearsky points for valid fit.')
            return None, None
        if isinstance(ax, type(None)):
            f, a= _plt.subplots()    
        else:
            a = ax
            f = a.get_figure()
        
        A = params.to_pandas().a
        b = params.to_pandas().b
        fit = A * np.power(self.mu0, b)
        fit.plot(ax = a, label = 'fit')
        return f,a

    @property
    def clearsky_parameters(self):
        return self._get_clearsky_parameters()
    
    @clearsky_parameters.setter
    def clearsky_parameters(self, cs_params: dict):
        self.dataset = self.dataset.drop_vars(['mask_normalized_global_magnitude',
                                               'mask_normalized_diffuse_ratio_variability',
                                               'mask_clear_sky_shortwave_radflux', 
                                               'clearsky_global_irradiation_powerlow_fit_params',
                                               'mask_diffuse_magnitude',
                                               'mask_global_irradiance_variability'], errors= 'ignore')
        for k,v in cs_params.items():
            self.dataset.attrs[k] = v

    def _get_clearsky_parameters(self, include_estimates = True):
        cs_params = ['mu0_min', 'nsw_exp', 'nsw_min', 'nsw_max', 'diffuse_max_coeff', 'diffuse_max_exp', 'max_dsw_dt', 'ndr_exp', 'ndr_std_max', 'ndr_window']
        if include_estimates:
            cs_params += ['ndr_std_max_estimated',
                            'diffuse_max_coeff_estimated',
                            'diffuse_max_exp_estimated',
                            'max_dsw_dt_estimated'
                            ]
        attr = self.dataset.attrs
        cs_params_dict = {k: attr[k] for k in attr if k in cs_params}
        return cs_params_dict




    
    @property
    def global_horizontal_irradiation(self):
        if isinstance(self._global_horizontal_irradiation, type(None)):
            self._global_horizontal_irradiation = GlobalHorizontalIrradiation(
                self.dataset,
                _dataset_ref = self._dataset_ref,
                **self._kwargs,
            )
        return self._global_horizontal_irradiation
    
    @property
    def diffuse_horizontal_irradiation(self):
        if isinstance(self._diffuse_horizontal_irradiation, type(None)):
            self._diffuse_horizontal_irradiation = DiffuseHorizontalIrradiation(
                self.dataset,
                _dataset_ref = self._dataset_ref,
                **self._kwargs,
            )
        return self._diffuse_horizontal_irradiation

    @property
    def direct_normal_irradiation(self):
        if isinstance(self._direct_normal_irradiation, type(None)):
            self._direct_normal_irradiation = DirectNormalIrradiation(
                self.dataset,
                _dataset_ref = self._dataset_ref,
                **self._kwargs,
            )
        return self._direct_normal_irradiation
    

    #TODO: come up with a better name
    @property
    def mask_normalized_diffuse_ratio_variability(self):
        if 'mask_normalized_diffuse_ratio_variability' not in self.dataset:
            if self.verbose:
                print('Running normalized diffuse ratio variability test')
            test_mask = atmcsk.normalized_diffuse_ratio_variability_test(global_irradiance=self.dataset.global_horizontal,
                                                                         diffuse_irradiance=self.dataset.diffuse_horizontal,
                                                                         mu0 = self.mu0,        
                                                                         mu0_min = self.get_attr('mu0_min'),
                                                                         ndr_exp = self.get_attr('ndr_exp'),
                                                                         ndr_std_max = self.get_attr('ndr_std_max'),
                                                                         window = self.get_attr('ndr_window'),)
            self.dataset['mask_normalized_diffuse_ratio_variability'] = test_mask
        return self.dataset.mask_normalized_diffuse_ratio_variability
    
    @property
    def mask_clear_sky_radflux(self) -> xr.DataArray:
        """
        Detect clear-sky periods in shortwave radiation using four Long & Ackerman–
        style tests.
    
        This implements **one iteration** of the clear-sky detection logic used in
        Radflux / SWFLUXANAL:
    
          1. Normalized shortwave magnitude test (NSW test)
          2. Diffuse magnitude test
          3. Change-with-time test on global SW
          4. Normalized diffuse ratio (NDR) variability test
    
        The thresholds and exponents are configurable via keyword arguments so that
        you can:
          - Use "generic" values for a first pass, or
          - Plug in iteration- and site-specific configuration later when you build
            the full Radflux-like system.
    

    
        Returns
        -------
        clear_mask : xr.DataArray
            Boolean mask along `time_dim` where True indicates clear-sky candidates
            according to all four tests.
        tests_dict : dict of xr.DataArray, optional
            Only if `return_tests` is True. Contains masks for each individual test:
              - "nsw"
              - "diffuse_mag"
              - "change_with_time"
              - "ndr_var"
    
        Notes
        -----
        - This is a **single iteration** detector; Radflux repeats detect–fit–
          interpolate cycles and updates thresholds based on fitted clear-sky
          functions. For now, you can treat this as the "iteration 0" or as a
          configurable building block.
        - Threshold values provided here are reasonable starting points but should
          be tuned/overridden to match the exact ARM/Radflux configuration.
        """
       
        # Combine all tests: clear if all tests pass
        if 'mask_clear_sky_shortwave_radflux' not in self.dataset:
            if self.verbose:
                print('Running clear sky tests (RADFLUX equivalent)')
            self.dataset['mask_clear_sky_shortwave_radflux'] = (self.global_horizontal_irradiation.mask_normalized_global_magnitude 
                                                                & self.diffuse_horizontal_irradiation.mask_magnitude 
                                                                & self.global_horizontal_irradiation.mask_variability 
                                                                & self.mask_normalized_diffuse_ratio_variability)
            self.dataset.mask_clear_sky_shortwave_radflux.attrs = {}
            
            self.dataset.mask_clear_sky_shortwave_radflux.attrs["info"] = "Radflux clear sky mask according to Long & Ackerman (2000) and subsequent publication iterations."
            self.dataset.mask_clear_sky_shortwave_radflux.attrs["unit"] = "1", 
            self.dataset.mask_clear_sky_shortwave_radflux.attrs["long_name"] = "clear sky classification mask",
            self.dataset.mask_clear_sky_shortwave_radflux.attrs["flag_values"] = '0, 1',
            self.dataset.mask_clear_sky_shortwave_radflux.attrs["flag_meanings"] = "0: fails radflux clear-sky test (cloudy), 1: passes radflux clear-sky test (possible clear-sky)"

        if ('clear_sky_params_optimized' not in self.dataset.attrs) or (self.dataset.attrs['clear_sky_params_optimized'] == 'False'):
            warnings.warn('Clear_sky parameters have not been optimized! It is recommended to run optimize_clearsky_parameters().')
        return self.dataset['mask_clear_sky_shortwave_radflux']



    def plot_overview(self, 
                    #   channel = 500, 
                      ax = None, 
                      apply_mask_clear_sky = False,
                      show_sunelevation = False,
                      plot_kwargs = {},):
        
        if isinstance(ax, type(None)):
            f, a= _plt.subplots()    
        else:
            a = ax
            f = a.get_figure()
        
        self.global_horizontal_irradiation
        self.diffuse_horizontal_irradiation
        self.direct_normal_irradiation
        
        dssel = self.dataset

        dssel.global_horizontal.plot(ax = a, label = 'global_horizontal')
        if apply_mask_clear_sky:
            g = a.get_lines()[-1]
            col = g.get_color()
            g.set_alpha(0.5)
            g.set_label('_nolegend_')
            dssel.global_horizontal.where(self.mask_clear_sky_radflux).plot(ax = a, label = 'global_horizontal', color = col)
        dssel.diffuse_horizontal.plot(ax = a, label = 'diffuse_horizontal')
        if apply_mask_clear_sky:
            g = a.get_lines()[-1]
            col = g.get_color()
            g.set_alpha(0.5)
            g.set_label('_nolegend_')
            dssel.diffuse_horizontal.where(self.mask_clear_sky_radflux).plot(ax = a, label = 'diffuse_horizontal', color = col)
        dssel.direct_normal.plot(ax = a, label = 'direct')
        if apply_mask_clear_sky:
            g = a.get_lines()[-1]
            col = g.get_color()
            g.set_alpha(0.5)
            g.set_label('_nolegend_')
            dssel.direct_normal.where(self.mask_clear_sky_radflux).plot(ax = a, label = 'direct_normal', color = col)
        
        if show_sunelevation:
            at = a.twinx()
            np.rad2deg(self.sun_position.elevation).plot(ax = at, color = 'black', ls = '--')
            # at.set_ylim(top = 0.9, bottom = 0)
        
        # a.set_xlim(left = pd.to_datetime('20220103 14:00:00'))
        a.grid()
        a.legend()
        return f,a

    def optimize_clearsky_parameters(self,
                                     n_iterations = 3,
                                     min_clear_for_update = 200,):
        
        self.dataset.attrs['clear_sky_params_optimized'] = 'Failed'
        if self.verbose:
            print((f'Original values -- nsw_exp:{self.get_attr('nsw_exp'):0.3f},'
                f'nsw_min: {self.get_attr('nsw_min'):0.1f},' 
                f'nsw_max: {self.get_attr('nsw_max'):0.1f},'
                f'ndr_exp: {self.get_attr('ndr_exp'):0.3f}'))

        for it in range(n_iterations):
            # 1 clear sky
            n_clear = int(self.mask_clear_sky_radflux.sum())
            if self.verbose:
                print('Number of clearsky (valid) points: ', n_clear)
            if n_clear < min_clear_for_update:
                if self.verbose:
                    print('Not enough clear sky points -> skip optimazation, keep old params')
                self.dataset.attrs['clear_sky_params_optimized'] = 'Not enough clear sky points'
                return None
            
            # 2. Update NSW thresholds/exponent
            
             #TODO is this even a thing? Do we need to catch bad fits?
            # if fit is None:
            #     if self.verbose:
            #         print('Fit failed: keep old params')
            #     assert(False), (params, None, None, None)
            
            
            mask_ngm = atmcsk.normalized_global_magnitude_test(self.dataset.global_horizontal.where(self.mask_clear_sky_radflux),
                                                            self.mu0,
                                                            mu0_min = self.get_attr('mu0_min'), 
                                                            nsw_exp = self.clearsky_global_irradiation_powerlow_fit_params.to_pandas().b,
                                                            nsw_min = self.get_attr('nsw_min'), #in the original those this and the following are not directly considered, they are indirectly considered throught the cloudmask though
                                                            nsw_max = self.get_attr('nsw_max'),
                                                        )
            

            
            # 3. update NDR exponent
            
            dgr_fit = atmcsk.fit_diffuse_global_ratio_mu0_powerlaw(self.mu0,
                                                                self.dataset.diffuse_horizontal,
                                                                self.dataset.global_horizontal,
                                                                self.mask_clear_sky_radflux, 
                                                                mu0_min = self.get_attr('mu0_min'),)

            self.dataset.attrs['nsw_exp'] = float(self.clearsky_global_irradiation_powerlow_fit_params.to_pandas().b)
            self.dataset.attrs['nsw_min'] = mask_ngm.nsw_min
            self.dataset.attrs['nsw_max'] = mask_ngm.nsw_max
            self.dataset.attrs['ndr_exp'] = float(dgr_fit.to_pandas().b)
            self.dataset.attrs['ndr_std_max_estimated'] = float(dgr_fit.ndr_std_max_estimated)
            if self.verbose:
                print((f'New values -- nsw_exp:{self.dataset.attrs['nsw_exp']:0.3f},'
                    f'nsw_min: {self.dataset.attrs['nsw_min']:0.1f},' 
                    f'nsw_max: {self.dataset.attrs['nsw_max']:0.1f},'
                    f'ndr_exp: {self.dataset.attrs['ndr_exp']:0.3f}'))
            self.dataset = self.dataset.drop_vars(['mask_normalized_global_magnitude',
                                                   'mask_normalized_diffuse_ratio_variability',
                                                    'mask_clear_sky_shortwave_radflux',
                                                    'mask_diffuse_magnitude',
                                                     'clearsky_global_irradiation_powerlow_fit_params',
                                                     'mask_global_irradiance_variability'])

        dcswd = atmcsk.fit_diffuse_mu0_powerlaw(self.mu0,self.dataset.diffuse_horizontal, self.mask_clear_sky_radflux, mu0_min = self.get_attr('mu0_min'))

        self.dataset.attrs['diffuse_max_coeff_estimated'] = float(dcswd.to_pandas().a * 1.2)
        self.dataset.attrs['diffuse_max_exp_estimated'] = float(dcswd.to_pandas().b)

        dsw_dt = self.dataset.global_horizontal.where(self.mask_clear_sky_radflux).differentiate('datetime') 
        ns_per_minute = np.float64(60 * 1e9)
        dsw_dt_per_min = dsw_dt * ns_per_minute
        dsw_dt_abs = np.abs(dsw_dt_per_min)
        self.dataset.attrs['max_dsw_dt_estimated'] = float(dsw_dt_abs.mean()*3)

        if self.verbose:
            print('suggested adjustements to user controlled configurations')
            print(f'\tdiffuse_max_coeff: {self.dataset.diffuse_max_coeff_estimated:0.2f} (is: {self.get_attr('diffuse_max_coeff')})')
            print(f'\tdiffuse_max_exp: {self.dataset.diffuse_max_exp_estimated:0.2f} (is: {self.get_attr('diffuse_max_exp')})')
            print(f'\tmax_dsw_dt: {self.dataset.max_dsw_dt_estimated:0.2f} (is: {self.get_attr('max_dsw_dt')})')
            print(f'\tndr_std_max: {self.dataset.ndr_std_max_estimated:0.4f} (is: {self.get_attr('ndr_std_max')})')

        self.dataset.attrs['clear_sky_params_optimized'] = 'True'