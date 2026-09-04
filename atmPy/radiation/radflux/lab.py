from __future__ import annotations

from dataclasses import asdict, dataclass
import pathlib as pl
import tomllib
from typing import Mapping
from importlib.resources import as_file, files
from atmPy.radiation.retrievals.broadband_shortwave_radiation import CombinedGlobalDiffuseDirect
import xarray as xr
import numpy as np
import sklearn
from atmPy.opt_imports import matplotlib as mpl
import warnings
import pandas as pd




path2defaul_sw_settings = (
    files("atmPy.radiation.radflux")
    / "resources"
    / "clear_sky_shortwave.example.toml"
)

SOLAR_CONSTANT = 1361.0


_ABBREVIATED_CLEAR_SKY_PARAMETER_NAMES = {
    "mu0_min": "maximum_solar_zenith_angle",
    "nsw_exp": "normalized_total_shortwave_power_exponent",
    "nsw_coeff": "normalized_total_shortwave_power_coefficient",
    "nsw_r2": "normalized_total_shortwave_coefficient_of_determination",
    "nsw_min": "normalized_total_shortwave_lower_limit_low_sun",
    "nsw_max": "normalized_total_shortwave_upper_limit",
    "diffuse_max_coeff": "maximum_diffuse_shortwave_irradiance",
    "diffuse_max_exp": "maximum_diffuse_shortwave_cosine_exponent",
    "normalized_diffuse_fit_coeff": "diffuse_shortwave_power_coefficient",
    "normalized_diffuse_fit_exp": "diffuse_shortwave_power_exponent",
    "max_dsw_dt": "shortwave_temporal_gradient_envelope_constant",
    "ndr_exp": "normalized_diffuse_ratio_power_exponent",
    "ndr_coeff": "normalized_diffuse_ratio_power_coefficient",
    "ndr_r2": "normalized_diffuse_ratio_coefficient_of_determination",
    "ndr_std_max": "normalized_diffuse_ratio_standard_deviation_limit",
    "ndr_window": "normalized_diffuse_ratio_standard_deviation_window_minutes",
    "ndr_std_max_estimated": "normalized_diffuse_ratio_standard_deviation_limit_estimate",
    "diffuse_max_coeff_estimated": "maximum_diffuse_shortwave_irradiance_estimate",
    "diffuse_max_exp_estimated": "maximum_diffuse_shortwave_cosine_exponent_estimate",
    "max_dsw_dt_estimated": "shortwave_temporal_gradient_envelope_constant_estimate",
    "clear_sky_params_optimized": "clear_sky_parameters_optimization_status",
}
_DESCRIPTIVE_TO_ABBREVIATED_CLEAR_SKY_PARAMETER_NAMES = {
    descriptive: abbreviated
    for abbreviated, descriptive in (
        _ABBREVIATED_CLEAR_SKY_PARAMETER_NAMES.items()
    )
}

_DESCRIPTIVE_CLEAR_SKY_PARAMETER_NAMES = (
    "maximum_solar_zenith_angle",
    "maximum_diffuse_shortwave_irradiance",
    "maximum_diffuse_shortwave_cosine_exponent",
    "normalized_diffuse_ratio_standard_deviation_limit",
    "normalized_diffuse_ratio_standard_deviation_window_minutes",
    "shortwave_temporal_gradient_envelope_constant",
    "minimum_total_shortwave_irradiance",
    "normalized_total_shortwave_intermediate_iteration_half_width",
    "normalized_total_shortwave_final_iteration_half_width",
    "normalized_total_shortwave_low_sun_cosine_boundary",
    "minimum_clear_sky_duration_for_daily_fit",
    "diffuse_ratio_exponent_lower_validity_limit",
    "diffuse_ratio_exponent_upper_validity_limit",
    "normalized_total_shortwave_power_exponent",
    "normalized_diffuse_ratio_power_exponent",
    "normalized_total_shortwave_lower_limit_low_sun",
    "normalized_total_shortwave_lower_limit_high_sun",
    "normalized_total_shortwave_upper_limit",
    "normalized_total_shortwave_power_coefficient",
    "normalized_total_shortwave_coefficient_of_determination",
    "normalized_diffuse_ratio_power_coefficient",
    "normalized_diffuse_ratio_coefficient_of_determination",
    "diffuse_shortwave_power_coefficient",
    "diffuse_shortwave_power_exponent",
)

_DESCRIPTIVE_CLEAR_SKY_ESTIMATE_NAMES = (
    "normalized_diffuse_ratio_standard_deviation_limit_estimate",
    "maximum_diffuse_shortwave_irradiance_estimate",
    "maximum_diffuse_shortwave_cosine_exponent_estimate",
    "shortwave_temporal_gradient_envelope_constant_estimate",
)


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
        mad = np.median(np.abs(residual - np.median(residual)))


        da =  xr.DataArray(
            np.array((A, b, r2, rmse, mae, median_abs_residual, outlier_fraction, mad), dtype=np.float64),
            dims=("fit_params_tcswd",),
            coords={"fit_params_tcswd": np.array(("a", "b", "r2","rmse","mae","median_abs_residual","outlier_fraction","mad"), dtype=object)},
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

@dataclass(frozen=True, slots=True)
class FixedParameters:
    maximum_solar_zenith_angle: float
    maximum_diffuse_shortwave_irradiance: float
    maximum_diffuse_shortwave_cosine_exponent: float
    normalized_diffuse_ratio_standard_deviation_limit: float
    normalized_diffuse_ratio_standard_deviation_window_minutes: int
    shortwave_temporal_gradient_envelope_constant: float
    minimum_total_shortwave_irradiance: float
    normalized_total_shortwave_intermediate_iteration_half_width: float
    normalized_total_shortwave_final_iteration_half_width: float
    normalized_total_shortwave_low_sun_cosine_boundary: float
    minimum_clear_sky_duration_for_daily_fit: int
    diffuse_ratio_exponent_lower_validity_limit: float
    diffuse_ratio_exponent_upper_validity_limit: float


@dataclass(frozen=True, slots=True)
class FitInitialValues:
    normalized_total_shortwave_power_exponent: float
    normalized_diffuse_ratio_power_exponent: float


@dataclass(frozen=True, slots=True)
class FirstIterationLimits:
    normalized_total_shortwave_lower_limit_low_sun: float
    normalized_total_shortwave_lower_limit_high_sun: float
    normalized_total_shortwave_upper_limit: float


@dataclass(frozen=True, slots=True)
class ClearSkyShortwaveSettings:
    fixed: FixedParameters
    initial: FitInitialValues
    first_iteration_limits: FirstIterationLimits
    legacy_names: Mapping[str, str]


def read_clear_sky_shortwave_settings(
    path: str | pl.Path = None,
) -> ClearSkyShortwaveSettings:
    """Read a clear-sky shortwave TOML settings file.

    Parameters
    ----------
    path : str | Path, optional
        Path to a settings (TOML) file. If None, the default settings file is used. Note, once loaded , values
        can not be changed. Make a copy of the default settings file 
        (atmPy.radiation.radflux.resources.clear_sky_shortwave.example.toml) and edit it to change values.
    """
    if path is None or path == "default":
        path = path2defaul_sw_settings

    path = pl.Path(path)

    with path.open("rb") as file:
        config = tomllib.load(file)

    shortwave = config["clear_sky_shortwave"]

    settings = ClearSkyShortwaveSettings(
        fixed=FixedParameters(**shortwave["fixed_parameters"]),
        initial=FitInitialValues(**shortwave["fit_initial_values"]),
        first_iteration_limits=FirstIterationLimits(
            **shortwave["first_iteration_limits"]
        ),
        legacy_names=dict(config.get("legacy_names", {})),
    )

    _validate_settings(settings)
    return settings


def _validate_settings(settings: ClearSkyShortwaveSettings) -> None:
    """Validate basic physical and algorithmic consistency."""
    fixed = settings.fixed
    initial = settings.initial
    limits = settings.first_iteration_limits

    if fixed.maximum_diffuse_shortwave_irradiance <= 0:
        raise ValueError(
            "maximum_diffuse_shortwave_irradiance must be positive."
        )

    if fixed.normalized_diffuse_ratio_standard_deviation_limit <= 0:
        raise ValueError(
            "normalized_diffuse_ratio_standard_deviation_limit "
            "must be positive."
        )

    if (
        fixed.normalized_diffuse_ratio_standard_deviation_window_minutes
        <= 0
    ):
        raise ValueError(
            "normalized_diffuse_ratio_standard_deviation_window_minutes "
            "must be positive."
        )

    if fixed.minimum_total_shortwave_irradiance <= 0:
        raise ValueError(
            "minimum_total_shortwave_irradiance must be positive."
        )

    if not 0 < fixed.normalized_total_shortwave_low_sun_cosine_boundary < 1:
        raise ValueError(
            "normalized_total_shortwave_low_sun_cosine_boundary "
            "must lie between 0 and 1."
        )

    if (
        fixed.normalized_total_shortwave_final_iteration_half_width
        > fixed.normalized_total_shortwave_intermediate_iteration_half_width
    ):
        raise ValueError(
            "The final-iteration NSW half-width should not exceed "
            "the intermediate-iteration half-width."
        )

    if (
        fixed.diffuse_ratio_exponent_lower_validity_limit
        >= fixed.diffuse_ratio_exponent_upper_validity_limit
    ):
        raise ValueError(
            "Diffuse-ratio exponent validity limits are reversed."
        )

    if not (
        limits.normalized_total_shortwave_lower_limit_low_sun
        < limits.normalized_total_shortwave_upper_limit
    ):
        raise ValueError("Invalid low-sun normalized-total-SW limits.")


    if initial.normalized_total_shortwave_power_exponent <= 0:
        raise ValueError(
            "normalized_total_shortwave_power_exponent must be positive."
        )


def _settings_parameters(
    settings: ClearSkyShortwaveSettings,
) -> dict[str, float | int]:
    return {
        **asdict(settings.fixed),
        **asdict(settings.initial),
        **asdict(settings.first_iteration_limits),
    }



class RadFlux(CombinedGlobalDiffuseDirect):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_attr(self, attr):
        """Return a dataset attribute, accepting legacy parameter names."""
        descriptive_name = _ABBREVIATED_CLEAR_SKY_PARAMETER_NAMES.get(
            attr, attr
        )
        attrs = self.dataset.attrs
        abbreviated_name = (
            _DESCRIPTIVE_TO_ABBREVIATED_CLEAR_SKY_PARAMETER_NAMES.get(
                descriptive_name
            )
        )

        if descriptive_name in attrs:
            value = attrs[descriptive_name]
        elif abbreviated_name in attrs:
            value = attrs[abbreviated_name]
            if descriptive_name == "maximum_solar_zenith_angle":
                value = np.rad2deg(np.arccos(value))
        else:
            if self.verbose:
                print(f'{descriptive_name} attribute is not set, using default')
            # self.dataset.attrs[attr] = default_config[attr]
            defaults = _settings_parameters(
                read_clear_sky_shortwave_settings()
            )
            if descriptive_name not in defaults:
                raise KeyError(
                    f"No clear-sky parameter named {descriptive_name!r}."
                )
            value = defaults[descriptive_name]
            attrs[descriptive_name] = value

        if attr == "mu0_min":
            value = np.cos(np.deg2rad(value))
        return value

    @property
    def clearsky_parameters(self):
        return self._get_clearsky_parameters()
    
    @clearsky_parameters.setter
    def clearsky_parameters(
        self,
        cs_params: Mapping | ClearSkyShortwaveSettings | str | pl.Path,
    ):
        """Sets the clearsky parameters in the dataset attributes. 
        Parameters
        ----------
        cs_params: Mapping, ClearSkyShortwaveSettings, str, or Path
            A mapping should use the descriptive clear-sky parameter names.
            Legacy abbreviated names are also accepted and converted.
            If str or Path, should be a path to a toml file containing the clearsky parameters.
            Note if cs_params ==  "default" the default clearsky parameters will be used.
        """
        self.reset_all_masks()
        if isinstance(cs_params, (str, pl.Path)):
            cs_params = read_clear_sky_shortwave_settings(cs_params)
            # params['nsw_coeff'] = cs_params.fixed.maximum_diffuse_shortwave_irradiance
            # params['nsw_r2'] = cs_params.fixed.maximum_diffuse_shortwave_cosine_exponent
            # params['normalized_diffuse_fit_coeff'] = cs_params.fixed.normalized_total_shortwave_final_iteration_half_width
            # params['normalized_diffuse_fit_exp'] = cs_params.fixed.normalized_total_shortwave_low_sun_cosine_boundary
            # minimum_clear_sky_duration_for_daily_fit
        if isinstance(cs_params, ClearSkyShortwaveSettings):
            cs_params = _settings_parameters(cs_params)
        # elif isinstance(cs_params, dict):
        elif not isinstance(cs_params, Mapping):
            raise TypeError(
                "clearsky_parameters must be a mapping, settings object, "
                "or TOML path."
            )

        params = {}
        for name, value in cs_params.items():
            descriptive_name = _ABBREVIATED_CLEAR_SKY_PARAMETER_NAMES.get(
                name, name
            )
            if name == "mu0_min":
                value = np.rad2deg(np.arccos(value))
            params[descriptive_name] = value

        for name, value in params.items():
            if name == 'normalized_diffuse_ratio_standard_deviation_window_minutes':
                if self.verbose:
                    print(
                        "Casting normalized diffuse ratio standard deviation "
                        "window to int"
                    )
                value = int(value)
            self.dataset.attrs[name] = value

        for abbreviated_name, descriptive_name in (
            _ABBREVIATED_CLEAR_SKY_PARAMETER_NAMES.items()
        ):
            if descriptive_name in params:
                self.dataset.attrs.pop(abbreviated_name, None)

    def _get_clearsky_parameters(self, include_estimates = True):
        parameter_names = list(_DESCRIPTIVE_CLEAR_SKY_PARAMETER_NAMES)
        if include_estimates:
            parameter_names.extend(_DESCRIPTIVE_CLEAR_SKY_ESTIMATE_NAMES)

        params = {}
        attrs = self.dataset.attrs
        for name in parameter_names:
            if name in attrs:
                params[name] = attrs[name]
                continue

            abbreviated_name = (
                _DESCRIPTIVE_TO_ABBREVIATED_CLEAR_SKY_PARAMETER_NAMES.get(name)
            )
            if abbreviated_name not in attrs:
                continue

            value = attrs[abbreviated_name]
            if abbreviated_name == "mu0_min":
                value = np.rad2deg(np.arccos(value))
            params[name] = value
        return params


    @property
    def mask_diffuse_magnitude(self):
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
        dif_max = self.get_attr('maximum_diffuse_shortwave_irradiance')
        dif_exp = self.get_attr(
            'maximum_diffuse_shortwave_cosine_exponent'
        )
        dif = self.dataset.diffuse_horizontal
        
        valid = (mu0 >= mu0_min) & dif.notnull()

        mu0_safe = mu0.where(mu0 > 0)
        dif_lim = dif_max * (mu0_safe ** dif_exp)

        test_mask = valid & (dif <= dif_lim)
        test_mask.name = "test_diffuse_mag"


        test_mask.attrs = {}
        test_mask.attrs["info"] = "Mask based on diffuse shortwave magnitude test. See Long & Ackerman (2000) and subsequent iterations for details."
        test_mask.attrs["unit"] = "1", 
        test_mask.attrs["long_name"] = "clear sky classification mask",
        test_mask.attrs["flag_values"] = '0, 1',
        test_mask.attrs["flag_meanings"] = "0: fails diffuse magnitude test (cloudy), 1: passes diffuse magnitude test (possible clear-sky)"
        out = {'mask': test_mask, 
               'auxiliary': {'dif_lim': dif_lim,
                             }}
        return out

    def plot_diffuse_magnitude_test(self, x = 'datetime'):
        self.mask_diffuse_magnitude # to initiate calculation
        if x == 'datetime':
            x = self.dataset.datetime
        elif x == 'mu0':
            x = self.mu0
        else:
            assert(False) 
        
        dif = self.dataset.diffuse_horizontal
        aux = self._diffuse_magnitude_test_auxiliary
        dif_lim = aux['dif_lim']

        f, aa = mpl.pyplot.subplots(2, sharex=True, height_ratios=[3,1], gridspec_kw={'hspace':0})
        #################
        a = aa[0]
        a.plot(x, dif, label = 'diffuse SW')
        a.plot(x, dif_lim, label = 'diffuse limit')
        a.set_ylim(0, dif_lim.max() * 1.1)
        a.set_ylabel('Diffuse SW [W m-2]')
        a.legend()
        a.set_title('Diffuse magnitude test')
        ##################
        a = aa[1]
        a.plot(x,self.mask_diffuse_magnitude,)
        a.set_ylabel('mask')
        #################
        a = aa[0]
        return f,aa

    @property
    def mask_global_irradiance_temporal_gradient(self,
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

    @property
    def mask_global_temporal_gradient(self):
        """Compatibility alias for the descriptive property name."""
        return self.mask_global_irradiance_temporal_gradient

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
        max_dsw_dt = self.get_attr(
            'shortwave_temporal_gradient_envelope_constant'
        )
        gsw = self.dataset.global_horizontal.copy(deep=True)
        # Differentiate wrt time; xarray returns W m-2 per nanosecond for datetime64,
        # so convert to per minute.
        dsw_dt_per_min = gsw.differentiate('datetime',
                                                 datetime_unit = 'm') 

        # Take absolute value and pad endpoints with NaNs
        dsw_dt_abs = np.abs(dsw_dt_per_min)
        dsw_dt_abs = dsw_dt_abs.reindex_like(gsw)  # align with original time axis

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
        self.mask_global_irradiance_temporal_gradient # just to initiate the calculation
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
        self.mask_global_irradiance_temporal_gradient.plot(ax = a,
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
        nsw_exp = self.get_attr(
            'normalized_total_shortwave_power_exponent'
        )
        # nsw_min = self.get_attr('nsw_min')
        nsw_min_low = self.get_attr(
            'normalized_total_shortwave_lower_limit_low_sun'
        )
        nsw_min_high = self.get_attr(
            'normalized_total_shortwave_lower_limit_high_sun'
        )
        nsw_max = self.get_attr('normalized_total_shortwave_upper_limit')
        mu0_boundary = self.get_attr(
            'normalized_total_shortwave_low_sun_cosine_boundary'
        )
        gsw = self.dataset.global_horizontal.copy(deep = True)
        mu0 = self.mu0
        valid = (mu0 >= mu0_min) & gsw.notnull()

        # Avoid division by zero
        mu0_safe = mu0.where(mu0 > 0)
        nsw = gsw / (mu0_safe ** nsw_exp)
        nsw = nsw.where(valid)

        # self.tp_nsw = nsw
        # self.tp_nsw_min_low = nsw_min_low
        # self.tp_nsw_min_high = nsw_min_high
        # self.tp_nsw_max = nsw_max

        nsw_min = xr.zeros_like(nsw, dtype=float )
        nsw_min = nsw_min.where(mu0 < mu0_boundary, other = nsw_min_high)
        nsw_min = nsw_min.where(mu0 > mu0_boundary, other = nsw_min_low)

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
        gsw = self.dataset.global_horizontal
        dsw = self.dataset.diffuse_horizontal
        mu0 = self.mu0
        mu0_min = self.get_attr('mu0_min')
        ndr_exp = self.get_attr(
            'normalized_diffuse_ratio_power_exponent'
        )
        ndr_std_max = self.get_attr(
            'normalized_diffuse_ratio_standard_deviation_limit'
        )
        window = self.get_attr(
            'normalized_diffuse_ratio_standard_deviation_window_minutes'
        )

        valid = (
            (mu0 >= mu0_min)
            & gsw.notnull()
            & dsw.notnull()
            & (gsw > 0)
        )
        
        mu0_safe = mu0.where(mu0 > 0)

        # Diffuse ratio
        dr = (dsw / gsw).where(valid)

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
        self.mask_clear_sky_shortwave_radflux # to initiate calculation
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
    def mask_clear_sky_shortwave_radflux(self) -> xr.DataArray:
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
            self.dataset['mask_clear_sky_shortwave_radflux'] = (self.mask_normalized_global_magnitude 
                                                                & self.mask_diffuse_magnitude 
                                                                & self.mask_global_irradiance_temporal_gradient
                                                                & self.mask_normalized_diffuse_ratio_variability)
            self.dataset.mask_clear_sky_shortwave_radflux.attrs = {}
            
            self.dataset.mask_clear_sky_shortwave_radflux.attrs["info"] = "Radflux clear sky mask according to Long & Ackerman (2000) and subsequent publication iterations."
            self.dataset.mask_clear_sky_shortwave_radflux.attrs["unit"] = "1", 
            self.dataset.mask_clear_sky_shortwave_radflux.attrs["long_name"] = "clear sky classification mask",
            self.dataset.mask_clear_sky_shortwave_radflux.attrs["flag_values"] = '0, 1',
            self.dataset.mask_clear_sky_shortwave_radflux.attrs["flag_meanings"] = "0: fails radflux clear-sky test (cloudy), 1: passes radflux clear-sky test (possible clear-sky)"

        optimization_status = self.dataset.attrs.get(
            'clear_sky_parameters_optimization_status'
        )
        if optimization_status is None:
            optimization_status = self.dataset.attrs.get(
                'clear_sky_params_optimized'
            )
        if optimization_status in (None, 'False'):
            warnings.warn('Clear-sky parameters have not been optimized! It is recommended to run optimize_clearsky_parameters().')
        return self.dataset['mask_clear_sky_shortwave_radflux']

    @property
    def mask_clear_sky_radflux(self) -> xr.DataArray:
        """Compatibility alias for the descriptive property name."""
        return self.mask_clear_sky_shortwave_radflux


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
        self.dataset.attrs[
            'clear_sky_parameters_optimization_status'
        ] = 'Failed'
        self.dataset.attrs.pop('clear_sky_params_optimized', None)
        if self.verbose:
            print(('################\n'
                    'Original values\n'
                   # normalized_total_shortwave_power_exponent
                   f'normalized_total_shortwave_power_exponent: '
                   f'{self.get_attr('normalized_total_shortwave_power_exponent'):0.3f},\n'
                   f'normalized_total_shortwave_lower_limit_low_sun: '
                   f'{self.get_attr('normalized_total_shortwave_lower_limit_low_sun'):0.1f},\n'
                   f'normalized_total_shortwave_lower_limit_high_sun: {self.get_attr('normalized_total_shortwave_lower_limit_high_sun'):0.1f},\n' 
                   f'normalized_total_shortwave_upper_limit: '
                   f'{self.get_attr('normalized_total_shortwave_upper_limit'):0.1f},\n'
                   f'normalized_diffuse_ratio_power_exponent: '
                   f'{self.get_attr('normalized_diffuse_ratio_power_exponent'):0.3f}'))
        weight_by_mu0 = False
        nsw_mu0_coverage = False
        nsw_final_mu0_coverage = False
        ndr_mu0_coverage = False
        for it in range(n_iterations):
            if self.verbose:
                print('################')
                print(f'Iteration {it+1}/{n_iterations}')
            if it == n_iterations - 1: # last iteration
                if self.verbose:
                    print('Set weight_by_mu0 = True for last iteration')
                weight_by_mu0 = True
                total = self.mask_clear_sky_shortwave_radflux.sum()
                above_th = self.mask_clear_sky_shortwave_radflux.where(self.mu0 > 0.6).sum()
                nsw_final_mu0_coverage = bool(above_th/total >= 0.45)
            else:
                # is mu0 larger than 80 of that at noon?
                nsw_mu0_coverage = bool(
                    self.mu0.where(self.mask_clear_sky_shortwave_radflux).max()
                    > 0.8 * self.mu0.max()
                )
            
            mu0_min = self.get_attr('mu0_min')
            ndr_mu0_coverage = bool(
                self.mu0.where(
                    self.mask_clear_sky_shortwave_radflux & (self.mu0 > mu0_min)
                ).min()
                < 0.4
            )


            # 1 check if sufficient clear sky points
            n_clear = int(self.mask_clear_sky_shortwave_radflux.sum())
            if self.verbose:
                print('Number of clearsky (valid) points: ', n_clear)
            if n_clear < min_clear_for_update:
                if self.verbose:
                    print(f'Not enough clear sky points ({n_clear} < {min_clear_for_update}) -> skip optimazation, keep old params')
                self.dataset.attrs[
                    'clear_sky_parameters_optimization_status'
                ] = 'Not enough clear sky points'
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
                                      mask_clearsky= self.mask_clear_sky_shortwave_radflux,
                                      mu0_min = self.get_attr('mu0_min'),
                                      min_points = min_clear_for_update,
                                      weight_by_mu0 = weight_by_mu0)
            params = res.fit_result
            self.normalized_total_shortwave_fit_results = res

            # update the parameters in the dataset attributes
            nsw_exp = float(params.sel(fit_params_tcswd='b'))
            nsw_coeff = float(params.sel(fit_params_tcswd='a'))
            nsw_r2 = float(params.sel(fit_params_tcswd='r2'))
            nsw_mad = float(params.sel(fit_params_tcswd='mad'))
            self.dataset.attrs[
                'normalized_total_shortwave_power_exponent'
            ] = nsw_exp
            self.dataset.attrs[
                'normalized_total_shortwave_power_coefficient'
            ] = nsw_coeff
            self.dataset.attrs[
                'normalized_total_shortwave_coefficient_of_determination'
            ] = nsw_r2
            self.dataset.attrs[
                'normalized_total_shortwave_median_absolute_deviation'
            ] = nsw_mad

            if it == n_iterations - 1:
                nsw_tol = self.get_attr(
                    'normalized_total_shortwave_final_iteration_half_width'
                )
            else:
                nsw_tol = self.get_attr(
                    'normalized_total_shortwave_intermediate_iteration_half_width'
                )
            self.dataset.attrs[
                'normalized_total_shortwave_lower_limit_low_sun'
            ] = nsw_coeff - nsw_tol
            self.dataset.attrs[
                'normalized_total_shortwave_lower_limit_high_sun'
            ] = nsw_coeff - nsw_tol
            self.dataset.attrs[
                'normalized_total_shortwave_upper_limit'
            ] = nsw_coeff + nsw_tol


                
            if isinstance(params, type(None)) and self.verbose:
                print('fit failed, probably not enough clear sky points')


            #########
            # 3. update NDR exponent

            res = fit_powerlaw_mu0(mu0 = self.mu0,
                                       values = self.dataset.diffuse_horizontal / self.dataset.global_horizontal,
                                       mask_clearsky = self.mask_clear_sky_shortwave_radflux,
                                       mu0_min = self.get_attr('mu0_min'),
                                       min_points = min_clear_for_update,
                                       weight_by_1_over_mu0 = weight_by_mu0
                                       )
            self.normalized_diffuse_ratio_fit_results = res
            params = res.fit_result

            # update the parameters in the dataset attributes
            ndr_exp = float(params.sel(fit_params_tcswd='b'))
            ndr_coeff = float(params.sel(fit_params_tcswd='a'))
            ndr_r2 = float(params.sel(fit_params_tcswd='r2'))
            ndr_mad = float(params.sel(fit_params_tcswd='mad'))
            self.dataset.attrs[
                'normalized_diffuse_ratio_power_exponent'
            ] = ndr_exp
            self.dataset.attrs[
                'normalized_diffuse_ratio_power_coefficient'
            ] = ndr_coeff
            self.dataset.attrs[
                'normalized_diffuse_ratio_coefficient_of_determination'
            ] = ndr_r2
            self.dataset.attrs[
                'normalized_diffuse_ratio_median_absolute_deviation'
            ] = ndr_mad

            if self.verbose:
                print(( 'New values\n'
                       f'normalized_total_shortwave_power_exponent: {nsw_exp:0.7f},\n'
                       f'normalized_total_shortwave_lower_limit_low_sun: '
                       f'{self.get_attr('normalized_total_shortwave_lower_limit_low_sun'):0.2f},\n'
                       f'normalized_total_shortwave_lower_limit_high_sun: {self.get_attr('normalized_total_shortwave_lower_limit_high_sun'):0.1f},\n' 
                       f'normalized_total_shortwave_upper_limit: '
                       f'{self.get_attr('normalized_total_shortwave_upper_limit'):0.2f},\n'
                       f'normalized_diffuse_ratio_power_exponent: {ndr_exp:0.7f}'))
            self.reset_all_masks() #this ensures that from now on all masks use the most up-to-date parameters.


        self.dataset.attrs[
            'clear_sky_parameters_optimization_status'
        ] = 'True'

        ##############
        ## extra, this does not really need to be here, as it does not add to the optimization, but it makes sence to have it here.
        res = fit_powerlaw_mu0(mu0 = self.mu0,
                            values = self.dataset.diffuse_horizontal,
                            mask_clearsky = self.mask_clear_sky_shortwave_radflux,
                            mu0_min = self.get_attr('mu0_min'),
                            min_points = min_clear_for_update
                                )
        self.diffuse_shortwave_fit_results = res
        params = res.fit_result

        self.dataset.attrs['diffuse_shortwave_power_coefficient'] = float(
            params.sel(fit_params_tcswd='a')
        )
        self.dataset.attrs['diffuse_shortwave_power_exponent'] = float(
            params.sel(fit_params_tcswd='b')
        )

        #### more tests
        # If b_diffr < −0.95 (too steep) → interpolate b_diffr.
        ndr_exp_min = self.get_attr(
            'diffuse_ratio_exponent_lower_validity_limit'
        )
        ndr_exp_max = self.get_attr(
            'diffuse_ratio_exponent_upper_validity_limit'
        )
        self.dataset.attrs[
            'diffuse_ratio_power_exponent_below_validity_limit'
        ] = ndr_exp < ndr_exp_min
        #If b_diffr > −0.4 (too flat) → interpolate both a_diffr and b_diffr.
        self.dataset.attrs[
            'diffuse_ratio_power_exponent_above_validity_limit'
        ] = ndr_exp > ndr_exp_max

        self.dataset.attrs[
            'normalized_total_shortwave_cosine_coverage_sufficient'
        ] = nsw_mu0_coverage
        self.dataset.attrs[
            'normalized_total_shortwave_final_iteration_cosine_coverage_sufficient'
        ] = nsw_final_mu0_coverage
        self.dataset.attrs[
            'diffuse_ratio_cosine_coverage_sufficient'
        ] = ndr_mu0_coverage

    



        if self.verbose:
            if not nsw_mu0_coverage:
                txt_nsw = f'{nsw_mu0_coverage}, Discard total-SW a,b; interpolate them.'
            else:
                txt_nsw = f'{nsw_mu0_coverage}, Keep total-SW a,b.'
            if not nsw_final_mu0_coverage:
                txt_nsw_final = f'{nsw_final_mu0_coverage}, Discard total-SW a,b; interpolate them.'
            else:
                txt_nsw_final = f'{nsw_final_mu0_coverage}, Keep total-SW a,b.'
            if not ndr_mu0_coverage:
                txt_diffuse = f'{ndr_mu0_coverage}, Discard only diffuse-ratio b; interpolate b. Keep diffuse-ratio a'
            else:
                txt_diffuse = f'{ndr_mu0_coverage}'
            if self.dataset.attrs['diffuse_ratio_power_exponent_below_validity_limit']:
                txt_ndr_too_steep = f'{self.dataset.attrs['diffuse_ratio_power_exponent_below_validity_limit']}, interpolate b.'
            else: 
                txt_ndr_too_steep = f'{self.dataset.attrs['diffuse_ratio_power_exponent_below_validity_limit']}, pass.'
            if self.dataset.attrs['diffuse_ratio_power_exponent_above_validity_limit']:
                txt_ndr_too_flat = f'{self.dataset.attrs['diffuse_ratio_power_exponent_above_validity_limit']}, interpolate a and b.'
            else:
                txt_ndr_too_flat = f'{self.dataset.attrs['diffuse_ratio_power_exponent_above_validity_limit']}, pass.'
            print((
                '### mu0 coverage test results\n'
                f'normalized_total_shortwave_cosine_coverage_sufficient: {txt_nsw},\n'
                f'normalized_total_shortwave_final_iteration_cosine_coverage_sufficient: {txt_nsw_final},\n'
                f'diffuse_ratio_cosine_coverage_sufficient: {txt_diffuse},\n'
                f'diffuse_ratio_power_exponent_below_validity_limit: {txt_ndr_too_steep},\n'
                f'diffuse_ratio_power_exponent_above_validity_limit: {txt_ndr_too_flat},\n'
            ))


        #todo: create an xarray dataset that contains the variables below (until line 1462). Determine values for the variables under conclusion based on the tests described above. When possible add attributes: unit(s), description, etc. This dataset is supposed to be my return object.
        #results
        normalized_total_shortwave_power_coefficient
        normalized_total_shortwave_power_exponent
        normalized_diffuse_ratio_power_coefficient
        normalized_diffuse_ratio_power_exponent

        # diagnostics
        n_clear
        n_clear_global_irradiance_termporal_gradiant
        n_clear_normalized_diffuse_ratio_variability
        n_clear_diffuse_magnitude
        n_clear_normalized_diffuse_ratio_variability

        mu0_coverage
        normalized_total_shortwave_median_absolute_deviation
        normalized_total_shortwave_coefficient_of_determination
        normalized_diffuse_ratio_median_absolute_deviation
        normalized_diffuse_ratio_coefficient_of_determination

        # tests
        diffuse_ratio_power_exponent_above_validity_limit
        diffuse_ratio_power_exponent_below_validity_limit
        normalized_total_shortwave_cosine_coverage_sufficient
        normalized_total_shortwave_final_iteration_cosine_coverage_sufficient
        diffuse_ratio_cosine_coverage_sufficient

        #conclusion
        normalized_total_shortwave_power_coefficient_is_valid
        normalized_total_shortwave_power_exponent_is_valid
        normalized_diffuse_ratio_power_coefficient_is_valid
        normalized_diffuse_ratio_power_exponent_is_valid




    def plot_clearsky_parameter_optimization_results(self):
        f,aa = mpl.pyplot.subplots(3, sharex=True, gridspec_kw={'hspace':0})
        f.set_figheight(f.get_figheight() * 1.5)
        # 1. nsw optimization results
        a = aa[0]
        res = self.normalized_total_shortwave_fit_results
        res.valid_points.plot(ax = a)
        res.fit.plot(ax = a)
        a.set_ylabel('Normalized global SW [W m-2]')
        a.set_xlabel('mu0 (cosine of solar zenith angle)')
        a.set_title('Normalized global power law fit')
        
        # 2. ndr optimization results
        a = aa[1]
        res = self.normalized_diffuse_ratio_fit_results
        res.valid_points.plot(ax = a)
        res.fit.plot(ax = a)
        a.set_ylabel('Normalized diffuse ratio')
        a.set_xlabel('mu0 (cosine of solar zenith angle)')
        a.set_title('Normalized diffuse ratio power law fit')

        # 3. diffuse fit results
        a = aa[2]
        res = self.diffuse_shortwave_fit_results
        res.valid_points.plot(ax = a)
        res.fit.plot(ax = a)
        a.set_ylabel('Normalized diffuse SW [W m-2]')
        a.set_xlabel('mu0 (cosine of solar zenith angle)')
        a.set_title('Normalized diffuse power law fit')
        return f,aa

    @property
    def clearsky_diffuse_horizontal(self):
        if 'clearsky_diffuse_horizontal' not in self.dataset:
            # params = self.clearsky_diffuse_irradiation_powerlow_fit_params
            params = self.clearsky_parameters
            assert(
                'diffuse_shortwave_power_coefficient' in params
                and 'diffuse_shortwave_power_exponent' in params
            ), (
                'Clear-sky diffuse shortwave power coefficient and exponent '
                'are not set. Either set them manually (e.g. from database) '
                'or run optimize_clearsky_parameters() first.'
            )
            a = params['diffuse_shortwave_power_coefficient']
            b = params['diffuse_shortwave_power_exponent']
            # A = params.to_pandas().a
            # b = params.to_pandas().b
            fit = a * np.power(self.mu0, b)
            self.dataset['clearsky_diffuse_horizontal'] = fit
            self.dataset.clearsky_diffuse_horizontal.attrs = {}
            self.dataset.clearsky_diffuse_horizontal.attrs['long_name'] = 'Clear sky diffuse horizontal irradiance (empirical power law fit)'
            self.dataset.clearsky_diffuse_horizontal.attrs['unit'] = 'W m-2'
        return self.dataset['clearsky_diffuse_horizontal']

    
    @property
    def clearsky_global_horizontal(self):
        if 'clearsky_global_horizontal' not in self.dataset:
            # params = self.clearsky_global_irradiation_powerlow_fit_params
            params = self.clearsky_parameters
            assert(
                'normalized_total_shortwave_power_coefficient' in params
                and 'normalized_total_shortwave_power_exponent' in params
            ), (
                'Clear-sky normalized total shortwave power coefficient and '
                'exponent are not set. Either set them manually (e.g. from '
                'database) or run optimize_clearsky_parameters() first.'
            )
            a = params['normalized_total_shortwave_power_coefficient']
            b = params['normalized_total_shortwave_power_exponent']
            fit = a * np.power(self.mu0, b)
            self.dataset['clearsky_global_horizontal'] = fit
            self.dataset.clearsky_global_horizontal.attrs = {}
            self.dataset.clearsky_global_horizontal.attrs['long_name'] = 'Clear sky global horizontal irradiance (empirical power law fit)'
            self.dataset.clearsky_global_horizontal.attrs['unit'] = 'W m-2'
        return self.dataset['clearsky_global_horizontal']
