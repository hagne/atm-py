from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib
from typing import Mapping
from importlib.resources import as_file, files

path2defaul_sw_settings = (
    files("atmPy.radiation.radflux")
    / "resources"
    / "clear_sky_shortwave.example.toml"
)

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
    normalized_total_shortwave_lower_limit: float
    normalized_total_shortwave_upper_limit: float


@dataclass(frozen=True, slots=True)
class ClearSkyShortwaveSettings:
    fixed: FixedParameters
    initial: FitInitialValues
    first_iteration_limits: FirstIterationLimits
    legacy_names: Mapping[str, str]


def read_clear_sky_shortwave_settings(
    path: str | Path = None,
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

    path = Path(path)

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
        limits.normalized_total_shortwave_lower_limit
        < limits.normalized_total_shortwave_upper_limit
    ):
        raise ValueError("Invalid low-sun normalized-total-SW limits.")


    if initial.normalized_total_shortwave_power_exponent <= 0:
        raise ValueError(
            "normalized_total_shortwave_power_exponent must be positive."
        )