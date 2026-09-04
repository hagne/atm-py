"""This is a collection of classes the are uses as the basis of many retrievals."""

import numpy as np
import sklearn
import xarray as xr
import pandas as pd
import atmPy.radiation.retrievals.clearsky as atmcsk
# import matplotlib.pyplot as _plt
from atmPy.opt_imports import matplotlib as mpl
import atmPy.general.measurement_site as atmgms
import warnings
from atmPy.radiation.retrievals import tiltcorrection
import pathlib as pl
import atmPy.radiation.radflux.lab as atmradflux

SOLAR_CONSTANT = 1361.0

# default_config = dict(# General filters
#                         mu0_min = 0.05,
#                         # normalized shortwave (nsw) magnitude test parameters
#                         nsw_exp = 1.2,
#                         nsw_min = 800.0,
#                         nsw_max = 1400.0,
#                         # Diffuse magnitude test parameters
#                         diffuse_max_coeff = 600,
#                         diffuse_max_exp = 0.5,
#                         # Change-with-time test parameters
#                         max_dsw_dt = 75.0,  # W m-2 per minute
#                         # NDR variability test parameters
#                         ndr_exp = -0.8,
#                         ndr_std_max = 0.005,
#                         ndr_window = 11,
#                         # # Output options
#                         # return_tests = False,
#                      )

def fit_powerlaw_mu0(
    mu0: xr.DataArray,
    values: xr.DataArray, 
    mask_clearsky: xr.DataArray,
    *,
    mu0_min: float,
    weight_by_mu0: bool = False,
    weight_by_1_over_mu0: bool = False,
    min_points: int = 100,
    verbose: bool = True) -> xr.DataArray | None:
    """
    Fit a simple power law to `values` using the robust HuberRegressor from sklearn. 
    The model is of the form:
    values = A * mu0^b
    using a linear regression in log space:
    log(values) = log(A) + b * log(mu0)
    over points where `mask_clearsky` is True, mu0 >= mu0_min, and
    values > 0.

    Parameters
    ----------
    weigh_by_mu0 : bool, optional
        If True, weight the regression by mu0. This is used for the global in the final iterations.
    weigh_by_1_over_mu0 : bool, optional
        If True, weight the regression by 1/mu0. This is used for the diffuse in the final iterations.
    Returns
    -------
    xr.DataArray or None
        Labeled output with:
        - tcswd_a: coefficient
        - tcswd_b: exponent
        - r_squared: coefficient of determination in log space
        None if not enough valid points.
    """
    cond = (
        mask_clearsky
        & (mu0 >= mu0_min)
        & mu0.notnull()
        & values.notnull()
        & (values > 0)
    )
    cond_vals = cond.values
    if int(cond_vals.sum()) < min_points:
        return None

    mu0_sel = mu0.values[cond_vals]
    y_sel = values.values[cond_vals]

    # Flatten and drop NaNs
    valid = np.isfinite(mu0_sel) & np.isfinite(y_sel) & (mu0_sel > 0) & (y_sel > 0)
    if valid.sum() < min_points:
        return None

    if 1:
        # robust linear regression
        assert(not (weight_by_mu0 and weight_by_1_over_mu0)), "Cannot weight by both mu0 and 1/mu0."
        if weight_by_mu0:
            if verbose:
                print("Weighting by mu0")
            # y_sel = y_sel * mu0_sel
            weights = mu0_sel[valid]
        elif weight_by_1_over_mu0:
            if verbose:
                print("Weighting by 1/mu0")
            # y_sel = y_sel / mu0_sel
            weights = 1 / mu0_sel[valid]
        else:
            weights = None
            pass

        x = np.log(mu0_sel[valid])
        z = np.log(y_sel[valid])

        fit = sklearn.linear_model.HuberRegressor().fit(x[:, None], z, 
                                                        sample_weight=weights
                                                        )

        b = fit.coef_[0]
        # if weight_by_mu0:
        #     b -= 1.0
        # elif weight_by_1_over_mu0:
        #     b += 1.0

        logA = fit.intercept_
        A = np.exp(logA)

        z_fit = fit.predict(x[:, None])
        residual = z - z_fit

        r2 = fit.score(x[:, None], z)
        rmse = np.sqrt(np.mean(residual**2))
        mae = np.mean(np.abs(residual))
        median_abs_residual = np.median(np.abs(residual))
        outlier_fraction = fit.outliers_.mean()


        da =  xr.DataArray(
            np.array((A, b, r2, rmse, mae, median_abs_residual, outlier_fraction), dtype=np.float64),
            dims=("fit_params_tcswd",),
            coords={"fit_params_tcswd": np.array(("a", "b", "r2","rmse","mae","median_abs_residual","outlier_fraction"), dtype=object)},
            name="global_powerlaw_mu0_fit",
        )
        da.attrs['info'] = 'Fit result for values = a * mu0^b under clearsky conditions.'
        ds = xr.Dataset()
        ds['fit_result'] = da
        ds['valid_points'] = xr.DataArray(z, coords = {'x': x})
        ds['fit'] = xr.DataArray(z_fit, coords = {'x': x})
        ds['residual'] = xr.DataArray(residual, coords = {'x': x})
    else:
        # Simple least-squares fit: z = log(A) + b * x
        n = x.size
        x_sum = x.sum()
        z_sum = z.sum()
        x_mean = x_sum / n
        z_mean = z_sum / n
        sxx = np.dot(x, x) - x_sum * x_mean
        sxz = np.dot(x, z) - x_sum * z_mean
        szz = np.dot(z, z) - z_sum * z_mean
        if sxx <= 0:
            return None

        b = sxz / sxx
        logA = z_mean - b * x_mean
        A = np.exp(logA)
        r2 = np.nan
        if szz > 0:
            r2 = float(np.clip((sxz * sxz) / (sxx * szz), 0.0, 1.0))

        da =  xr.DataArray(
            np.array((A, b, r2), dtype=np.float64),
            dims=("fit_params_tcswd",),
            coords={"fit_params_tcswd": np.array(("a", "b", "r2"), dtype=object)},
            name="global_powerlaw_mu0_fit",
        )
        da.attrs['info'] = 'Fit result for global_irradiance = a * mu0^b under clearsky conditions.'

    return ds

class _DatasetRef(object):
    def __init__(self, dataset):
        """A simple class to hold a reference to a dataset. This is used to allow multiple classes to share the same dataset without having to pass it around explicitly."""
        self.dataset = dataset


class SolarIrradiation(object):
    def __init__(self, dataset, site = None, verbose = False, _dataset_ref = None):
        """Base class for solar irradiation datasets. Provides shared functionality and properties for global horizontal, diffuse horizontal, and direct normal irradiation.
        Parameters
        ----------
        dataset : xr.Dataset
            A dataset that contains the relevant irradiation variable(s) and the coordinates 'datetime'.
        site : atmPy.general.measurement_sites.Station, atmPy.general.measurement_sites.MovingPlatform, optional
            A measurement site object that can be used to calculate sun position and other site-specific parameters. 
            """
        self.verbose = verbose  
        self._dataset_ref = _dataset_ref if _dataset_ref is not None else _DatasetRef(dataset)
        if _dataset_ref is None:
            self.dataset = dataset
        self._site = site
        self._sun_position_variables = ['zenith', 'zenith_geometric', 'elevation_geometric', 'elevation', 'azimuth', 'equation_of_time', 'airmass', 'airmass_absolute', 'sun_earth_distance']
        self._sun_position_variables_ds = [f'solar_{v}' for v in self._sun_position_variables] # for internal use, to avoid name clashes with potential variables in the dataset.
        assert('datetime' in dataset or 'datetime' in dataset.dims), 'Time coordinate has to be called datetime .... sorry, i know that is an unconventional choise for the time cooridinate name.'

    @property
    def temporal_resolution_minutes(self):
        """Returns the temporal resolution of the dataset in minutes. This is computed as the median difference between consecutive time points."""
        return float((self.dataset.datetime - self.dataset.datetime.shift(datetime = 1)).dropna('datetime').median()/pd.to_timedelta(1, 'm'))

    @property
    def site(self):
        if isinstance(self._site, type(None)):
            if self.verbose:
                print('No site provided, trying to infer site information from dataset.')
            if 'latitude' in self.dataset.attrs and 'longitude' in self.dataset.attrs:
                if self.verbose:
                    print('Found latitude and longitude in dataset attributes. Assuming fixed station.')
                if 'altitude' in self.dataset.attrs:
                    alt = self.dataset.attrs['altitude']
                else:
                    alt = 0
                self._site = atmgms.Station(lat = self.dataset.attrs['latitude'], lon = self.dataset.attrs['longitude'], alt = alt)
            elif 'latitude' in self.dataset and 'longitude' in self.dataset:
                if self.verbose:
                    print('Found latitude and longitude in dataset variables. ', end = '')
                if 'altitude' in self.dataset:
                    alt = self.dataset.altitude
                else:
                    alt = 0
                if 'datetime' in self.dataset.latitude.dims:
                    if self.verbose:
                        print('Latitude and longitude are time dependent. Assuming moving platform.')
                    self._site = atmgms.MovingPlatform(lat = self.dataset['latitude'], lon = self.dataset['longitude'], alt = alt)
                else:
                    if self.verbose:
                        print('Latitude and longitude are not time dependent. Assuming fixed station.')
                    self._site = atmgms.Station(lat = self.dataset['latitude'], lon = self.dataset['longitude'], alt = alt)

            else:
                raise ValueError('No site information found. Set site keyword or provide latitude and longitude in dataset attributes or variables. Note, time dependent latitude and longitude will be interpreted as a moving platform.')
        return self._site

            
            
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
            # self.dataset.attrs[attr] = default_config[attr]
            self.clearsky_parameters = 'default'
        return self.dataset.attrs[attr]

    @property
    def mu0(self):
        if 'mu0' not in self.dataset:
            self.dataset['mu0'] = np.cos(self.sun_position.solar_zenith) 
        return self.dataset['mu0']
    
    @property
    def sun_position(self):
        if not np.all([v in self.dataset for v in self._sun_position_variables_ds]):
        # if isinstance(self._sun_position, type(None)):
            sp = self.site.get_sun_position(self.dataset.datetime)
            self.tp_sp = sp
            for v in self._sun_position_variables:
                self.dataset[f'solar_{v}'] = sp[v]
        return self.dataset[self._sun_position_variables_ds]
    
    @sun_position.setter
    def sun_position(self, value):
        self.dataset.datetime.identical(value.datetime)
        for v in self._sun_position_variables:
            self.dataset[f'solar_{v}'] = value[v]



class GlobalHorizontalIrradiation(SolarIrradiation):
    def __init__(self, dataset, **kwargs):
        super().__init__(dataset, **kwargs)
        assert('global_horizontal' in dataset), f'global_horizontal variable is missing.'

    @property
    def mask_temporal_gradient(self,
                                            # global_irradiance: xr.DataArray,
                                            # time_dim: str,
                                            # max_dsw_dt: float,
                                        ):
        if 'mask_global_irradiance_temporal_gradient' not in self.dataset:
            if self.verbose:
                print('Running global irradiance temporal gradient test')
            self._global_irradiance_temporal_gradient_test( # self.dataset.global_horizontal, 
                                                            # max_dsw_dt = self.get_attr('max_dsw_dt')
                                                            )
            # self.dataset['mask_global_irradiance_temporal_gradient'] = out['mask']
            # self._mask_global_irradiance_temporal_gradient = out['auxiliary']
        return self.dataset.mask_global_irradiance_temporal_gradient

    def _global_irradiance_temporal_gradient_test(self,
        # global_irradiance: xr.DataArray,
        # *,
        # max_dsw_dt: float,
        # # time_dim: str = "datetime",
    ) -> xr.DataArray:
        """
        #todo: adjust docstring to reflect recent changes

        Originally called global_irradiance_change_with_time_test

        Temporal gradient test on global shortwave. Note, this is not a variability test as the 
        normalized_diffuse_ratio_variability_test. This is a test that looks for abrupt changes
        in terms of the rate, that is the derivative.

        Long & Ackerman compare the rate of change of surface SW to that of TOA
        SW; here we implement a generic magnitude limit on |d(sw_global)/dt|:

            |d sw_global / dt| <= max_dsw_dt

        where dt is computed from the time coordinate.

        Parameters
        ----------
        sw_global : xr.DataArray
            Downwelling global shortwave flux [W m-2].
        max_dsw_dt : float
            Maximum allowed |d(sw_global)/dt| in W m-2 per minute.

        Returns
        -------
        xr.DataArray
            Boolean mask where True indicates the point passes this test.
            Endpoints (first/last point) are treated as failing if derivative
            cannot be computed.
        """
        max_dsw_dt = self.get_attr('max_dsw_dt')
        global_irradiance = self.dataset.global_horizontal.copy(deep=True)
        # Differentiate wrt time; xarray returns W m-2 per nanosecond for datetime64,
        # so convert to per minute.
        dsw_dt_per_min = global_irradiance.differentiate('datetime',
                                                 datetime_unit = 'm') 

        # Take absolute value and pad endpoints with NaNs
        dsw_dt_abs = np.abs(dsw_dt_per_min)
        dsw_dt_abs = dsw_dt_abs.reindex_like(global_irradiance)  # align with original time axis

        # build limits based on top of the atmosphere radiation
        I_t = (SOLAR_CONSTANT / self.sun_position.solar_sun_earth_distance**2) / self.sun_position.solar_airmass
        dI_t = abs(I_t.differentiate('datetime', 
                         datetime_unit = 'm'
                        ))
           
        lim_max = dI_t + max_dsw_dt * 1/self.sun_position.solar_airmass
        mu_noon = 1/self.sun_position.solar_airmass.min()
        lim_min = dI_t - (self.temporal_resolution_minutes * (mu_noon + 0.01) * self.sun_position.solar_airmass)

        test_mask = (dsw_dt_abs <= lim_max) & (dsw_dt_abs >= lim_min) & dsw_dt_abs.notnull()
        test_mask.name = "test_global_irradiance_temporal_gradient"


        test_mask.attrs = {}
        test_mask.attrs["info"] = "Mask based on change-with-time test on global shortwave (global variability test). See Long & Ackerman (2000) and subsequent iterations for details."
        test_mask.attrs["unit"] = "1", 
        test_mask.attrs["long_name"] = "clear sky classification mask",
        test_mask.attrs["flag_values"] = '0, 1',
        test_mask.attrs["flag_meanings"] = "0: fails change-with-time test (cloudy), 1: passes change-with-time test (possible clear-sky)"
        axiliary = {'dsw_dt_abs': dsw_dt_abs}
        axiliary['lim_max'] = lim_max
        axiliary['lim_min'] = lim_min
        out = {'mask': test_mask, 'auxiliary': axiliary}
        self.dataset['mask_global_irradiance_temporal_gradient'] = test_mask
        self._mask_global_irradiance_temporal_gradient_auxiliary = axiliary
        return out

    def plot_global_irradiance_temporal_gradient_test(self, 
                                                    #   ax = None
                                                      ):
        """similar to Fig. 3 in Long & Ackerman 2000"""
        self.mask_temporal_gradient # just to initiate the calculation 
        # if isinstance(ax, type(None)):
        f, aa= mpl.pyplot.subplots(2, sharex = True, height_ratios = [3,1], gridspec_kw = {'hspace':0} )

        # the test
        a = aa[0]
        lim_max = self._mask_global_irradiance_temporal_gradient_auxiliary['lim_max']
        lim_min = self._mask_global_irradiance_temporal_gradient_auxiliary['lim_min']
        dsw_dt_abs = self._mask_global_irradiance_temporal_gradient_auxiliary['dsw_dt_abs']
        dsw_dt_abs.plot(ax = a, label = '|dSW/dt|')
        a.fill_between(lim_max.datetime, lim_max, lim_min, alpha = 0.4, label = 'limit envelope')
        a.set_ylabel('dSW/dt [W m-2 per minute]')
        a.legend()
        a.set_title('Global irradiance temporal gradient test')

        # the mask
        # at = a.twinx(zorder = -1)
        a = aa[1]
        self.mask_temporal_gradient.plot(ax = a, 
                                        #  color = 'black', zorder = -1
                                         )
        a.set_ylabel('mask')
        # a.set_ylim(0,10)
        return f,aa


            
    @property
    def mask_normalized_global_magnitude(self):
        if 'mask_normalized_global_magnitude' not in self.dataset:
            if self.verbose:
                print('Running normalized global magnitude test')
            out = self._normalized_global_magnitude_test() #self.dataset.global_horizontal,
                                                                                                        # self.mu0,
                                                                                                        # mu0_min = self.get_attr('mu0_min'), 
                                                                                                        # nsw_exp = self.get_attr('nsw_exp'),
                                                                                                        # nsw_min = self.get_attr('nsw_min'),
                                                                                                        # nsw_max = self.get_attr('nsw_max'),)


            self.dataset['mask_normalized_global_magnitude'] = out['mask']                        
            self._normalized_global_magnitude_test_auxiliary = out['auxiliary']
        return self.dataset.mask_normalized_global_magnitude


    def _normalized_global_magnitude_test(self,
        # sw_global: xr.DataArray,
        # mu0: xr.DataArray,
        # *,
        # mu0_min: float,
        # nsw_exp: float,
        # nsw_min: float,
        # nsw_max: float,
    ) -> xr.DataArray:
        """
        Normalized shortwave magnitude test (NSW test).

        Implements a simplified version of the Long & Ackerman 'normalized SW'
        constraint:

            NSW = sw_global / mu0**nsw_exp

        A sample passes this test if:
            nsw_min <= NSW <= nsw_max   and   mu0 >= mu0_min

        Parameters
        ----------
        sw_global : xr.DataArray
            Downwelling global shortwave flux on a horizontal surface [W m-2].
        mu0 : xr.DataArray
            Cosine of solar zenith angle (unitless).
        mu0_min : float
            Minimum mu0 required to consider points (exclude very low sun).
        nsw_exp : float
            Exponent used in the normalization (often ~1.2 in initial iterations).
        nsw_min, nsw_max : float
            Lower and upper allowed bounds for NSW (W m-2 * (unitless)^(-nsw_exp)).

        Returns
        -------
        xr.DataArray
            Boolean mask where True indicates the point passes this test.
        """
        auxiliary = {}
        mu0_min = self.get_attr('mu0_min')
        nsw_exp = self.get_attr('nsw_exp')
        # nsw_min = self.get_attr('nsw_min')
        nsw_min_low = self.get_attr('nsw_min')
        nsw_min_high = self.get_attr('normalized_total_shortwave_lower_limit_high_sun')


        nsw_max = self.get_attr('nsw_max')
        sw_global = self.dataset.global_horizontal.copy(deep = True)
        mu0 = self.mu0
        valid = (mu0 >= mu0_min) & sw_global.notnull()

        # Avoid division by zero
        mu0_safe = mu0.where(mu0 > 0)
        nsw = sw_global / (mu0_safe ** nsw_exp)
        nsw = nsw.where(valid)

        # self.tp_nsw = nsw
        # self.tp_nsw_min_low = nsw_min_low
        # self.tp_nsw_min_high = nsw_min_high
        # self.tp_nsw_max = nsw_max

        nsw_min = xr.zeros_like(nsw, dtype=float )
        nsw_min = nsw_min.where(self.mu0 < 0.2, other = nsw_min_high)
        nsw_min = nsw_min.where(self.mu0 > 0.2, other = nsw_min_low)

        test_mask = valid & (nsw >= nsw_min) & (nsw <= nsw_max)
        auxiliary['nsw_min'] = nsw_min
        auxiliary['nsw_max'] = nsw_max
        auxiliary['nsw'] = nsw
        test_mask.name = "test_nsw"

        # Robust center and spread
        median = float(np.nanmedian(nsw))
        q25, q75 = np.nanpercentile(nsw, [25, 75])
        iqr = float(q75 - q25)
        if not np.isfinite(iqr) or iqr <= 0:
            iqr = max(0.1 * median, 10.0)  # fall-back
        iqr_k = 3
        nsw_min = max(median - iqr_k * iqr, 0.0)
        nsw_max = median + iqr_k * iqr

        test_mask.attrs = {}
        test_mask.attrs["info"] = "Mask based on normalized shortwave magnitude (NSW) test. See Long & Ackerman (2000) and subsequent iterations for details."
        test_mask.attrs["unit"] = "1", 
        test_mask.attrs["long_name"] = "clear sky classification mask",
        test_mask.attrs["flag_values"] = '0, 1',
        test_mask.attrs["flag_meanings"] = "0: fails NSW test (cloudy), 1: passes NSW test (possible clear-sky)",
        test_mask.attrs["nsw_min"] = nsw_min
        test_mask.attrs["nsw_max"] = nsw_max
        out = dict(mask = test_mask,
                   auxiliary = auxiliary
                   )
        

        return out

    def plot_normalized_global_magnitude_test(self):
        self.mask_normalized_global_magnitude # to initiate calculation
        aux = self._normalized_global_magnitude_test_auxiliary
        nsw_min = aux['nsw_min']
        nsw_max = aux['nsw_max']
        nsw = aux['nsw']
        f, aa = mpl.pyplot.subplots(2, sharex=True, height_ratios=[3,1], gridspec_kw={'hspace':0})
        #################
        a = aa[0]
        nsw.plot(ax = a)
        a.fill_between(nsw_min.datetime,nsw_min, 
                nsw_min.where(False, other=nsw_max), 
                alpha = 0.4, zorder = 0)
        a.set_ylim(nsw_min.min()*0.9, nsw_max* 1.1)
        a.set_ylabel('Normalized SW [W m-2]')
        a.set_title('Normalized global magnitude test')

        ##################
        a = aa[1]
        self.mask_normalized_global_magnitude.plot(ax = a)
        a.set_ylabel('mask')
        #################
        a = aa[0]
        return f,aa




class DiffuseHorizontalIrradiation(SolarIrradiation):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert('diffuse_horizontal' in self.dataset), f'diffuse_horizontal variable is missing.'

    @property
    def mask_magnitude(self):
        if 'mask_diffuse_magnitude' not in self.dataset:
            if self.verbose:
                print('Running diffuse magnitude test')
            out = self._diffuse_magnitude_test(
                                                                                    # self.dataset.diffuse_horizontal,
                                                                                    # self.mu0,
                                                                                    # mu0_min = self.get_attr('mu0_min'),
                                                                                    # diffuse_max_coeff = self.get_attr('diffuse_max_coeff'),
                                                                                    # diffuse_max_exp = self.get_attr('diffuse_max_exp'),
                                                                                    )
            self.dataset['mask_diffuse_magnitude'] = out['mask']
            self._diffuse_magnitude_test_auxiliary = out['auxiliary']
        return self.dataset.mask_diffuse_magnitude

    def _diffuse_magnitude_test(self,
                            #     diffuse_irradiance: xr.DataArray,
                            # mu0: xr.DataArray,
                            # *,
                            # mu0_min: float,
                            # diffuse_max_coeff: float,
                            # diffuse_max_exp: float,
                            ) -> dict:
        """
        Diffuse shortwave magnitude test.

        Conceptually follows Long & Ackerman's requirement that diffuse SW
        not exceed a mu0-dependent envelope under clear-sky:

            sw_diffuse <= diffuse_max_coeff * mu0**diffuse_max_exp

        Parameters
        ----------
        sw_diffuse : xr.DataArray
            Downwelling diffuse shortwave flux [W m-2].
        mu0 : xr.DataArray
            Cosine of solar zenith angle (unitless).
        mu0_min : float
            Minimum mu0 required to consider points (exclude very low sun).
        diffuse_max_coeff : float
            Coefficient setting the magnitude of the diffuse envelope.
        diffuse_max_exp : float
            Exponent controlling how the envelope scales with mu0.

        Returns
        -------
        xr.DataArray
            Boolean mask where True indicates the point passes this test.
        """
        mu0 = self.mu0
        mu0_min = self.get_attr('mu0_min')
        diffuse_max_coeff = self.get_attr('diffuse_max_coeff')
        diffuse_max_exp = self.get_attr('diffuse_max_exp')
        diffuse_irradiance = self.dataset.diffuse_horizontal
        
        valid = (mu0 >= mu0_min) & diffuse_irradiance.notnull()

        mu0_safe = mu0.where(mu0 > 0)
        diffuse_limit = diffuse_max_coeff * (mu0_safe ** diffuse_max_exp)

        test_mask = valid & (diffuse_irradiance <= diffuse_limit)
        test_mask.name = "test_diffuse_mag"


        test_mask.attrs = {}
        test_mask.attrs["info"] = "Mask based on diffuse shortwave magnitude test. See Long & Ackerman (2000) and subsequent iterations for details."
        test_mask.attrs["unit"] = "1", 
        test_mask.attrs["long_name"] = "clear sky classification mask",
        test_mask.attrs["flag_values"] = '0, 1',
        test_mask.attrs["flag_meanings"] = "0: fails diffuse magnitude test (cloudy), 1: passes diffuse magnitude test (possible clear-sky)"
        out = {'mask': test_mask, 
               'auxiliary': {'diffuse_limit': diffuse_limit,
                             }}
        return out

    def plot_normalized_global_magnitude_test(self, x = 'datetime'):
        self.mask_magnitude # to initiate calculation
        if x == 'datetime':
            x = self.dataset.datetime
        elif x == 'mu0':
            x = self.mu0
        else:
            assert(False) 
        
        diffuse = self.dataset.diffuse_horizontal
        aux = self._diffuse_magnitude_test_auxiliary
        diffuse_lim = aux['diffuse_limit']

        f, aa = mpl.pyplot.subplots(2, sharex=True, height_ratios=[3,1], gridspec_kw={'hspace':0})
        #################
        a = aa[0]
        a.plot(x, diffuse, label = 'diffuse SW')
        a.plot(x,diffuse_lim, label = 'diffuse limit')
        a.set_ylim(0, diffuse_lim.max() * 1.1)
        a.set_ylabel('Diffuse SW [W m-2]')
        a.legend()
        a.set_title('Diffuse magnitude test')
        ##################
        a = aa[1]
        a.plot(x,self.mask_magnitude,)
        a.set_ylabel('mask')
        #################
        a = aa[0]
        return f,aa



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
            ds['direct_normal'] = ds.direct_horizontal / xr.DataArray(np.cos(self.sun_position.solar_zenith))
            ds.direct_normal.attrs.update(ds.direct_horizontal.attrs)
            ds.direct_normal.attrs["long_name"] += ' (calculated from direct_horizontal and solar zenith angle)'
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
    def clearsky_global_horizontal(self):
        if 'clearsky_global_horizontal' not in self.dataset:
            # params = self.clearsky_global_irradiation_powerlow_fit_params
            assert('nsw_coeff' in self.dataset.attrs and 'nsw_exp' in self.dataset.attrs), 'clearsky parameters nsw_coeff and nsw_exp are not set. Either set them manually (e.g. from database) or run optimize_clearsky_parameters() first.'
            A = self.dataset.attrs['nsw_coeff']
            b = self.dataset.attrs['nsw_exp']
            fit = A * np.power(self.mu0, b)
            self.dataset['clearsky_global_horizontal'] = fit
            self.dataset.clearsky_global_horizontal.attrs = {}
            self.dataset.clearsky_global_horizontal.attrs['long_name'] = 'Clear sky global horizontal irradiance (empirical power law fit)'
            self.dataset.clearsky_global_horizontal.attrs['unit'] = 'W m-2'
        return self.dataset['clearsky_global_horizontal']
    
    # @property
    # def clearsky_global_irradiation_powerlow_fit_params(self):
    #     """Performes a empirical power law fit to the clear sky global irradiance
    #     and returns the fit parameters."""

    #     if 'clearsky_global_irradiation_powerlow_fit_params' not in self.dataset:
    #         if self.verbose:
    #             print('performing powerlow fit for clear sky global irradiance')
    #         params = atmcsk.fit_global_powerlaw_mu0(self.mu0, self.dataset.global_horizontal, self.mask_clear_sky_radflux, 
    #                                                     mu0_min = self.get_attr('mu0_min'),
    #                                                     min_points = 100)
    #         if isinstance(params, type(None)) and self.verbose:
    #             print('fit failed, probably not enough clear sky points')
                
    #         self.dataset['clearsky_global_irradiation_powerlow_fit_params'] = params
    #     return self.dataset['clearsky_global_irradiation_powerlow_fit_params']

    # @property
    # def clearsky_diffuse_irradiation_powerlow_fit_params(self):
    #     """Performes a empirical power law fit to the clear sky diffuse irradiance
    #     and returns the fit parameters."""

    #     if 'clearsky_diffuse_irradiation_powerlow_fit_params' not in self.dataset:
    #         if self.verbose:
    #             print('performing powerlow fit for clear sky diffuse irradiance')
    #         dt_in_m = np.median(self.dataset.datetime.values[1:]-self.dataset.datetime.values[:-1])/pd.to_timedelta(1,'m')

    #         params = atmcsk.fit_global_powerlaw_mu0(self.mu0, self.dataset.diffuse_horizontal, self.mask_clear_sky_radflux, 
    #                                                     mu0_min = self.get_attr('mu0_min'),
    #                                                     min_points = dt_in_m * 100)
    #         if isinstance(params, type(None)) and self.verbose:
    #             print('fit failed, probably not enough clear sky points')
                
    #         self.dataset['clearsky_diffuse_irradiation_powerlow_fit_params'] = params
    #     return self.dataset['clearsky_diffuse_irradiation_powerlow_fit_params']
    
    @property
    def clearsky_diffuse_horizontal(self):
        if 'clearsky_diffuse_horizontal' not in self.dataset:
            # params = self.clearsky_diffuse_irradiation_powerlow_fit_params
            assert('normalized_diffuse_fit_coeff' in self.dataset.attrs and 'normalized_diffuse_fit_exp' in self.dataset.attrs), 'clearsky parameters ndiff_coeff and ndiff_exp are not set. Either set them manually (e.g. from database) or run optimize_clearsky_parameters() first.'
            A = self.dataset.attrs['normalized_diffuse_fit_coeff']
            b = self.dataset.attrs['normalized_diffuse_fit_exp']
            # A = params.to_pandas().a
            # b = params.to_pandas().b
            fit = A * np.power(self.mu0, b)
            self.dataset['clearsky_diffuse_horizontal'] = fit
            self.dataset.clearsky_diffuse_horizontal.attrs = {}
            self.dataset.clearsky_diffuse_horizontal.attrs['long_name'] = 'Clear sky diffuse horizontal irradiance (empirical power law fit)'
            self.dataset.clearsky_diffuse_horizontal.attrs['unit'] = 'W m-2'
        return self.dataset['clearsky_diffuse_horizontal']

    # def plot_clearsky_global_irradiation_powerlow_fit(self, ax = None):
    #     params = self.clearsky_global_irradiation_powerlow_fit_params
    #     if any(params.isnull()):
    #         print('Not enough clearsky points for valid fit.')
    #         return None, None
    #     if isinstance(ax, type(None)):
    #         f, a= mpl.pyplot.subplots()    
    #     else:
    #         a = ax
    #         f = a.get_figure()
        
    #     A = params.to_pandas().a
    #     b = params.to_pandas().b
    #     fit = A * np.power(self.mu0, b)
    #     fit.plot(ax = a, label = 'fit')
    #     return f,a

    @property
    def clearsky_parameters(self):
        return self._get_clearsky_parameters()
    
    @clearsky_parameters.setter
    def clearsky_parameters(self, cs_params: dict|str|pl.Path):
        """Sets the clearsky parameters in the dataset attributes. 
        Parameters
        ----------
        cs_params: dict, str, or Path
            If dict, should contain the clearsky parameters as keys and their values. 
            If str or Path, should be a path to a toml file containing the clearsky parameters.
            Note if cs_params ==  "default" the default clearsky parameters will be used.
        """
        self.reset_all_masks()
        if isinstance(cs_params, (str, pl.Path)):
            cs_params = atmradflux.read_clear_sky_shortwave_settings(cs_params)
            self.tp_sp = cs_params
            params = {}
            params['mu0_min'] = np.cos(np.deg2rad(cs_params.fixed.maximum_solar_zenith_angle))
            params['nsw_exp'] = cs_params.initial.normalized_total_shortwave_power_exponent
            # params['nsw_coeff'] = cs_params.fixed.maximum_diffuse_shortwave_irradiance
            # params['nsw_r2'] = cs_params.fixed.maximum_diffuse_shortwave_cosine_exponent
            params['nsw_min'] = cs_params.first_iteration_limits.normalized_total_shortwave_lower_limit_low_sun
            params['normalized_total_shortwave_lower_limit_high_sun'] = cs_params.first_iteration_limits.normalized_total_shortwave_lower_limit_high_sun
            params['nsw_max'] = cs_params.first_iteration_limits.normalized_total_shortwave_upper_limit
            params['diffuse_max_coeff'] = cs_params.fixed.maximum_diffuse_shortwave_irradiance
            params['diffuse_max_exp'] = cs_params.fixed.maximum_diffuse_shortwave_cosine_exponent
            # params['normalized_diffuse_fit_coeff'] = cs_params.fixed.normalized_total_shortwave_final_iteration_half_width
            # params['normalized_diffuse_fit_exp'] = cs_params.fixed.normalized_total_shortwave_low_sun_cosine_boundary
            params['max_dsw_dt'] = cs_params.fixed.shortwave_temporal_gradient_envelope_constant
            # minimum_clear_sky_duration_for_daily_fit
            params['ndr_exp'] = cs_params.initial.normalized_diffuse_ratio_power_exponent
            params['ndr_std_max'] = cs_params.fixed.normalized_diffuse_ratio_standard_deviation_limit
            params['ndr_window'] = cs_params.fixed.normalized_diffuse_ratio_standard_deviation_window_minutes
            cs_params = params
            
        # elif isinstance(cs_params, dict):
        for k,v in cs_params.items():
            if k== 'ndr_window':
                if self.verbose:
                    print('casting ndr_window to int')
                v = int(v)
            self.dataset.attrs[k] = v

    def _get_clearsky_parameters(self, include_estimates = True):
        cs_params = ['mu0_min', 'nsw_exp','nsw_coeff', 'nsw_r2',
                     'nsw_min', 
                     'normalized_total_shortwave_lower_limit_high_sun',
                     'nsw_max', 'diffuse_max_coeff', 'diffuse_max_exp', 'normalized_diffuse_fit_coeff', 'normalized_diffuse_fit_exp', 'max_dsw_dt', 'ndr_exp', 'ndr_std_max', 'ndr_window']
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
            self.run_normalized_diffuse_ratio_variability_test()
        return self.dataset.mask_normalized_diffuse_ratio_variability



    def run_normalized_diffuse_ratio_variability_test(self,                                                                  
                                                        mu0_min = None,
                                                        ndr_exp = None,
                                                        ndr_std_max = None,
                                                        window = None
                                                        ) -> dict:
        """
        Normalized diffuse ratio (NDR) variability test.

        This mirrors the core idea from Long & Ackerman:
            - Compute diffuse ratio DR = diffuse_irradiance / global_irradiance
            - Normalize by a power of mu0: NDR = DR * mu0**(-ndr_exp)
            - Compute rolling-window std(NDR); clear-sky requires that this std
            remains below a small threshold.

        Parameters
        ----------
        global_irradiance : xr.DataArray
            Downwelling global shortwave [W m-2].
        diffuse_irradiance : xr.DataArray
            Downwelling diffuse shortwave [W m-2].
        mu0 : xr.DataArray
            Cosine of solar zenith angle (unitless).
        mu0_min : float
            Minimum mu0 to consider (exclude very low sun).
        ndr_exp : float
            Exponent used in the normalization (often around -0.8).
        ndr_std_max : float
            Maximum allowed standard deviation of NDR in the rolling window.
        window : int
            Rolling window length in **number of samples**.
        time_dim : str
            Name of the time dimension.

        Returns
        -------
        xr.DataArray
            Boolean mask where True indicates the point passes the NDR variability test.
        """
        global_irradiance=self.dataset.global_horizontal
        diffuse_irradiance=self.dataset.diffuse_horizontal
        mu0 = self.mu0
        mu0_min = self.get_attr('mu0_min')
        ndr_exp = self.get_attr('ndr_exp')
        ndr_std_max = self.get_attr('ndr_std_max')
        window = self.get_attr('ndr_window')

        valid = (
            (mu0 >= mu0_min)
            & global_irradiance.notnull()
            & diffuse_irradiance.notnull()
            & (global_irradiance > 0)
        )
        
        mu0_safe = mu0.where(mu0 > 0)

        # Diffuse ratio
        dr = (diffuse_irradiance / global_irradiance).where(valid)

        # Normalized diffuse ratio
        ndr = dr * (mu0_safe ** (-ndr_exp))
        ndr = ndr.compute()  # ensure it's not lazy from now on
        ndr.name = "ndr"

        # Rolling std over the chosen window (centered to mimic a symmetric window)
        # window is in samples; ensure it's odd for symmetry if you care
        ndr_roll = ndr.rolling({'datetime': window}, center=True)# TODO variable name time will need to be adjusted
        ndr_std = ndr_roll.std()  
        ndr_mean = ndr_roll.mean()

        test_mask = (ndr_std <= ndr_std_max) & ndr_std.notnull() & valid
        test_mask.name = "test_ndr_var"

        test_mask.attrs = {}
        test_mask.attrs["info"] = "Mask based on normalized diffuse ratio (NDR) variability test. See Long & Ackerman (2000) and subsequent iterations for details."
        test_mask.attrs["unit"] = "1", 
        test_mask.attrs["long_name"] = "clear sky classification mask",
        test_mask.attrs["flag_values"] = '0, 1',
        test_mask.attrs["flag_meanings"] = "0: fails NDR variability test (cloudy), 1: passes NDR variability test (possible clear-sky)"

        out = {
            'mask': test_mask,
            'normalized_diffuse_ratio': ndr,
            'normalized_diffuse_ratio_std': ndr_std,
            'normalized_diffuse_ratio_mean': ndr_mean,
            'normalized_diffuse_ratio_std_max': ndr_std_max,
            'window': window,
        }
        self.dataset['mask_normalized_diffuse_ratio_variability'] = out['mask']
        self._normalized_diffuse_ratio_variability_test_results = out
        return out


    def plot_normalized_diffuse_ratio_variability_test(self):
        self.mask_clear_sky_radflux # to initiate calculation
        res = self._normalized_diffuse_ratio_variability_test_results
        f, aa = mpl.pyplot.subplots(3, sharex=True, gridspec_kw={'hspace':0})
        f.set_figheight(f.get_figheight() * 1.5)

        ndr = res['normalized_diffuse_ratio']
        ndr_std = res['normalized_diffuse_ratio_std']
        ndr_mean = res['normalized_diffuse_ratio_mean']
        ndr_std_max = res['normalized_diffuse_ratio_std_max']

        ###################
        a = aa[0]
        ndr.plot(ax = a, color = '0.5', lw = 0.8, label = 'NDR')
        ndr_mean.plot(ax = a, label = 'NDR mean')
        a.fill_between(ndr.datetime, ndr_mean - ndr_std, ndr_mean + ndr_std, alpha = 0.4, zorder = 0, label = 'NDR Mean +- std')
        # fb = a.fill_between(ndr.datetime, ndr_mean - ndr_std_max, ndr_mean + ndr_std_max, alpha = 0.4, zorder = 0)
        for i in (ndr_mean - ndr_std_max, ndr_mean + ndr_std_max):
            i.plot(color = 'red', ls = '--', ax = a, label = 'NDR mean +- max std')
        g = a.get_lines()[-1]
        g.set_label('no_legend')

        # dlim = ndr_std_max * 2
        dlim = ndr_mean.std() *2
        center = ndr_mean.median()
        ymin = center-dlim
        ymax = center+dlim
        a.set_ylim(ymin, ymax)
        a.set_ylabel('Normalized diffuse ratio')
        a.set_title('Normalized diffuse ratio variability test')
        #################
        a = aa[1]
        ndr_std.plot(ax= a)
        a.axhline(ndr_std_max, color = 'red', ls = '--', label = 'max std')
        a.set_ylim(0, ndr_std_max * 1.2)
        a.set_label('NDR std')
        a.legend()
        ###################
        a = aa[2]
        self.mask_normalized_diffuse_ratio_variability.plot(ax = a)
        a.set_ylabel('mask')
        return f, aa


    


    
    def reset_all_masks(self):
        """Resets all masks in the dataset. This is useful if you want to re-run the clear-sky detection with different parameters.
        Parameters
        ----------
        error_when_missing : str, optional
            If 'raise', an error will be raised if a mask is not found. If '"""

        for var in ['mask_normalized_global_magnitude',
                    'mask_normalized_diffuse_ratio_variability',
                    'mask_clear_sky_shortwave_radflux', 
                    'mask_diffuse_magnitude',
                    'mask_global_irradiance_temporal_gradient']:
            if var not in self.dataset:
                if self.verbose:
                    print(f'Warning: {var} not found in dataset. Skipping.')
            else:
                self.dataset = self.dataset.drop_vars(var)
                if self.verbose:
                    print(f'Reset {var} in dataset.')

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
                                                                & self.global_horizontal_irradiation.mask_temporal_gradient 
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
                      show_clearsky = False, 
                      apply_mask_clear_sky = False,
                      show_sunelevation = False,
                      plot_kwargs = {},):
        
        if isinstance(ax, type(None)):
            f, a= mpl.pyplot.subplots()    
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
        
        if show_clearsky:
            self.clearsky_global_horizontal.plot(ax = a, ls = '--', label = 'clearsky_global_horizontal')
            self.clearsky_diffuse_horizontal.plot(ax = a, ls = '--', label = 'clearsky_diffuse_horizontal')

        if show_sunelevation:
            at = a.twinx()
            np.rad2deg(self.sun_position.elevation).plot(ax = at, color = 'black', ls = '--')
            # at.set_ylim(top = 0.9, bottom = 0)
        
        # a.set_xlim(left = pd.to_datetime('20220103 14:00:00'))
        a.grid()
        a.legend()
        return f,a

    def optimize_clearsky_parameters(self,
                                     n_iterations = 4,
                                     min_clear_for_update_equivalent = 100,): #todo: Only for minute data. Self adjusting based on the time resolution of the data would be better.
        """Optimizes the clear sky parameters.
        Parameters
        ----------
        n_iterations : int, optional
            Number of iterations to run the optimization. Default is 4.
        min_clear_for_update_equivalent : int, optional
            Equivalent minimum number of clear sky points required for updating the parameters. This value is being adjusted based on the time resolution of the data.
            The number is the equivalent for minute data, which was the original resolution of the Radflux algorithm. Default is 100."""
        
        # adjust the number of clear sky points based on the time resolution of the data. 
        dt_in_m = np.median(self.dataset.datetime.values[1:]-self.dataset.datetime.values[:-1])/pd.to_timedelta(1,'m')
        min_clear_for_update = int(min_clear_for_update_equivalent * dt_in_m)  
        self.dataset.attrs['clear_sky_params_optimized'] = 'Failed'
        if self.verbose:
            print(('################\n'
                    'Original values\n'
                   f'nsw_exp:{self.get_attr('nsw_exp'):0.3f},\n' #normalized_total_shortwave_power_exponent
                   f'nsw_min: {self.get_attr('nsw_min'):0.1f},\n' 
                   f'normalized_total_shortwave_lower_limit_high_sun: {self.get_attr('normalized_total_shortwave_lower_limit_high_sun'):0.1f},\n' 
                   f'nsw_max: {self.get_attr('nsw_max'):0.1f},\n'
                   f'ndr_exp: {self.get_attr('ndr_exp'):0.3f}'))
        weight_by_mu0 = False
        for it in range(n_iterations):
            if self.verbose:
                print('################')
                print(f'Iteration {it+1}/{n_iterations}')
            if it == n_iterations - 1: # last iteration
                if self.verbose:
                    print('Set weight_by_mu0 = True for last iteration')
                weight_by_mu0 = True
                total = self.mask_clear_sky_radflux.sum()
                above_th = self.mask_clear_sky_radflux.where(self.mu0 > 0.6).sum()
                test_mu0_coverage_NSW_final_iteration = bool(above_th/total >= 0.45)
            else:
                # is mu0 larger than 80 of that at noon?
                test_mu0_coverage_NSW = bool(self.mu0.where(self.mask_clear_sky_radflux).max() > 0.8 * self.mu0.max())
            
            mu0_min = self.get_attr('mu0_min')
            test_mu0_coverage_diffuse_ratio = bool(self.mu0.where(self.mask_clear_sky_radflux & (self.mu0 > mu0_min)).min() < 0.4)


            # 1 check if sufficient clear sky points
            n_clear = int(self.mask_clear_sky_radflux.sum())
            if self.verbose:
                print('Number of clearsky (valid) points: ', n_clear)
            if n_clear < min_clear_for_update:
                if self.verbose:
                    print(f'Not enough clear sky points ({n_clear} < {min_clear_for_update}) -> skip optimazation, keep old params')
                self.dataset.attrs['clear_sky_params_optimized'] = 'Not enough clear sky points'
                return None
            
            #####
            # 2. Fit global with powerlaw -> Update Normalized shortwave (nsw) thresholds/exponent
             #TODO is this even a thing? Do we need to catch bad fits?
            # if fit is None:
            #     if self.verbose:
            #         print('Fit failed: keep old params')
            #     assert(False), (params, None, None, None)
            
            # params = atmcsk.fit_global_powerlaw_mu0(self.mu0, self.dataset.global_horizontal, self.mask_clear_sky_radflux, 
            #                                 mu0_min = self.get_attr('mu0_min'),
            #                                 min_points = min_clear_for_update) #todo: valid for minute data only. this should be adjustable, such that minute and second resolution data can be used as well. 
            res = fit_powerlaw_mu0(mu0 = self.mu0, 
                                      values = self.dataset.global_horizontal, 
                                      mask_clearsky= self.mask_clear_sky_radflux,
                                      mu0_min = self.get_attr('mu0_min'),
                                      min_points = min_clear_for_update,
                                      weight_by_mu0 = weight_by_mu0)
            params = res.fit_result
            self.tp_res_nsw = res

            # update the parameters in the dataset attributes
            self.dataset.attrs['nsw_exp'] = float(params.to_pandas().b)
            self.dataset.attrs['nsw_coeff'] = float(params.to_pandas().a)
            self.dataset.attrs['nsw_r2'] = float(params.to_pandas().r2)
            if it == n_iterations - 1:
                nsw_tol = 100
            else:
                nsw_tol = 150
            self.dataset.attrs['nsw_min'] = self.dataset.attrs['nsw_coeff'] - nsw_tol
            self.dataset.attrs['normalized_total_shortwave_lower_limit_high_sun'] = self.dataset.attrs['nsw_coeff'] - nsw_tol
            self.dataset.attrs['nsw_max'] = self.dataset.attrs['nsw_coeff'] + nsw_tol


                
            if isinstance(params, type(None)) and self.verbose:
                print('fit failed, probably not enough clear sky points')


            #########
            # 3. update NDR exponent

            res = fit_powerlaw_mu0(mu0 = self.mu0,
                                       values = self.dataset.diffuse_horizontal / self.dataset.global_horizontal,
                                       mask_clearsky = self.mask_clear_sky_radflux,
                                       mu0_min = self.get_attr('mu0_min'),
                                       min_points = min_clear_for_update,
                                       weight_by_1_over_mu0 = weight_by_mu0
                                       )
            self.tp_res_ndr = res
            params = res.fit_result

            # update the parameters in the dataset attributes
            self.dataset.attrs['ndr_exp'] = float(params.to_pandas().b)
            self.dataset.attrs['ndr_coeff'] = float(params.to_pandas().a)
            self.dataset.attrs['ndr_r2'] = float(params.to_pandas().r2)

            if self.verbose:
                print(( 'New values\n'
                       f'nsw_exp:{self.dataset.attrs['nsw_exp']:0.7f},\n'
                       f'nsw_min: {self.dataset.attrs['nsw_min']:0.2f},\n'                                    
                       f'normalized_total_shortwave_lower_limit_high_sun: {self.get_attr('normalized_total_shortwave_lower_limit_high_sun'):0.1f},\n' 
                       f'nsw_max: {self.dataset.attrs['nsw_max']:0.2f},\n'
                       f'ndr_exp: {self.dataset.attrs['ndr_exp']:0.7f}'))
            self.reset_all_masks() #this ensures that from now on all masks use the most up-to-date parameters.

        self.dataset.attrs['clear_sky_params_optimized'] = 'True'

        ##############
        ## extra, this does not really need to be here, as it does not add to the optimization, but it makes sence to have it here.
        res = fit_powerlaw_mu0(mu0 = self.mu0,
                            values = self.dataset.diffuse_horizontal,
                            mask_clearsky = self.mask_clear_sky_radflux,
                            mu0_min = self.get_attr('mu0_min'),
                            min_points = min_clear_for_update
                                )
        self.tp_res_diffuse = res
        params = res.fit_result

        self.dataset.attrs['normalized_diffuse_fit_coeff'] = float(params.to_pandas().a)
        self.dataset.attrs['normalized_diffuse_fit_exp'] = float(params.to_pandas().b)

        #### more tests
        # If b_diffr < −0.95 (too steep) → interpolate b_diffr.
        self.dataset.attrs['test_ndr_exp_too_steep'] = self.dataset.attrs['ndr_exp'] < -0.95
        #If b_diffr > −0.4 (too flat) → interpolate both a_diffr and b_diffr.
        self.dataset.attrs['test_ndr_exp_too_flat'] = self.dataset.attrs['ndr_exp'] > -0.4

        self.dataset.attrs['test_mu0_coverage_NSW'] = test_mu0_coverage_NSW
        self.dataset.attrs['test_mu0_coverage_NSW_final_iteration'] = test_mu0_coverage_NSW_final_iteration
        self.dataset.attrs['test_mu0_coverage_diffuse_ratio'] = test_mu0_coverage_diffuse_ratio

        if self.verbose:
            if not test_mu0_coverage_NSW:
                txtNSW = f'{test_mu0_coverage_NSW}, Discard total-SW a,b; interpolate them.'
            else:
                txtNSW = f'{test_mu0_coverage_NSW}, Keep total-SW a,b.'
            if not test_mu0_coverage_NSW_final_iteration:
                txtNSWfinal = f'{test_mu0_coverage_NSW_final_iteration}, Discard total-SW a,b; interpolate them.'
            else:
                txtNSWfinal = f'{test_mu0_coverage_NSW_final_iteration}, Keep total-SW a,b.'
            if not test_mu0_coverage_diffuse_ratio:
                txt_diffuse = f'{test_mu0_coverage_diffuse_ratio}, Discard only diffuse-ratio b; interpolate b. Keep diffuse-ratio a'
            else:
                txt_diffuse = f'{test_mu0_coverage_diffuse_ratio}'
            if self.dataset.attrs['test_ndr_exp_too_steep']:
                txt_ndr_too_steep = f'{self.dataset.attrs['test_ndr_exp_too_steep']}, interpolate b.'
            else: 
                txt_ndr_too_steep = f'{self.dataset.attrs['test_ndr_exp_too_steep']}, pass.'
            if self.dataset.attrs['test_ndr_exp_too_flat']:
                txt_ndr_too_flat = f'{self.dataset.attrs['test_ndr_exp_too_flat']}, interpolate a and b.'
            else:
                txt_ndr_too_flat = f'{self.dataset.attrs['test_ndr_exp_too_flat']}, pass.'
            print(('### mu0 coverage test results\n'
                    f'test_mu0_coverage_NSW: {txtNSW},\n'
                    f'test_mu0_coverage_NSW_final_iteration: {txtNSWfinal},\n'
                    f'test_mu0_coverage_diffuse_ratio: {txt_diffuse},\n'
                    f'test_ndr_exp_too_steep: {txt_ndr_too_steep},\n'
                    f'test_ndr_exp_too_flat: {txt_ndr_too_flat},\n'))
            
        f,aa = mpl.pyplot.subplots(3, sharex=True, gridspec_kw={'hspace':0})
        f.set_figheight(f.get_figheight() * 1.5)
        # 1. nsw optimization results
        a = aa[0]
        res = self.tp_res_nsw
        res.valid_points.plot(ax = a)
        res.fit.plot(ax = a)
        a.set_ylabel('Normalized global SW [W m-2]')
        a.set_xlabel('mu0 (cosine of solar zenith angle)')
        a.set_title('Normalized global power law fit')
        
        # 2. ndr optimization results
        a = aa[1]
        res = self.tp_res_ndr
        res.valid_points.plot(ax = a)
        res.fit.plot(ax = a)
        a.set_ylabel('Normalized diffuse ratio')
        a.set_xlabel('mu0 (cosine of solar zenith angle)')
        a.set_title('Normalized diffuse ratio power law fit')

        # 3. diffuse fit results
        a = aa[2]
        res = self.tp_res_diffuse
        res.valid_points.plot(ax = a)
        res.fit.plot(ax = a)
        a.set_ylabel('Normalized diffuse SW [W m-2]')
        a.set_xlabel('mu0 (cosine of solar zenith angle)')
        a.set_title('Normalized diffuse power law fit')
        return f,aa



    def apply_tilt_correction(self, sensor_response_time=None):
        # call irradiance properties to make sure they are initialized and data is present in the dataset
        # self.direct_normal_irradiation  
        self.sun_position
        ds = tiltcorrection.apply_tilt_correction(self.dataset, sensor_response_time=sensor_response_time)
        out = CombinedGlobalDiffuseDirect(ds)
        return out

    apply_tilt_correction.__doc__ = tiltcorrection.apply_tilt_correction.__doc__
