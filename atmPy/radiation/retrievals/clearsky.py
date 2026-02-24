"""This is a collection of clear sky tests"""

import xarray as xr
import numpy as np
import pandas as pd

# TODO: come up with a better name
def normalized_diffuse_ratio_variability_test(
    global_irradiance: xr.DataArray,
    diffuse_irradiance: xr.DataArray,
    mu0: xr.DataArray,
    *,
    mu0_min: float,
    ndr_exp: float,
    ndr_std_max: float,
    window: int,
    # time_dim: str,
) -> xr.DataArray:
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
    ndr.name = "ndr"

    # Rolling std over the chosen window (centered to mimic a symmetric window)
    # window is in samples; ensure it's odd for symmetry if you care
    ndr_std = (
        ndr.rolling({'datetime': window}, center=True) #TODO variable name time will need to be adjusted
        .std()
    )

    test_mask = (ndr_std <= ndr_std_max) & ndr_std.notnull() & valid
    test_mask.name = "test_ndr_var"

    test_mask.attrs = {}
    test_mask.attrs["info"] = "Mask based on normalized diffuse ratio (NDR) variability test. See Long & Ackerman (2000) and subsequent iterations for details."
    test_mask.attrs["unit"] = "1", 
    test_mask.attrs["long_name"] = "clear sky classification mask",
    test_mask.attrs["flag_values"] = '0, 1',
    test_mask.attrs["flag_meanings"] = "0: fails NDR variability test (cloudy), 1: passes NDR variability test (possible clear-sky)"
    return test_mask

#TODO: the following should not just be for global, this should be applicable to any radiation, just the shresholds would change.
def global_irradiance_variability_test(self,
    global_irradiance: xr.DataArray,
    *,
    max_dsw_dt: float,
    # time_dim: str = "datetime",
) -> xr.DataArray:
    """
    Change-with-time test on global shortwave.

    Long & Ackerman compare the rate of change of surface SW to that of TOA
    SW; here we implement a generic magnitude limit on |d(sw_global)/dt|:

        |d sw_global / dt| <= max_dsw_dt

    where dt is computed from the time coordinate.

    Parameters
    ----------
    sw_global : xr.DataArray
        Downwelling global shortwave flux [W m-2].
    time_dim : str
        Name of the time dimension in `sw_global`.
    max_dsw_dt : float
        Maximum allowed |d(sw_global)/dt| in W m-2 per minute.

    Returns
    -------
    xr.DataArray
        Boolean mask where True indicates the point passes this test.
        Endpoints (first/last point) are treated as failing if derivative
        cannot be computed.
    """

    
    # Differentiate wrt time; xarray returns W m-2 per nanosecond for datetime64,
    # so convert to per minute.
    dsw_dt = global_irradiance.differentiate('datetime') 
    # Convert units: ns -> minutes
    # 1 minute = 60 s = 60 * 1e9 ns

    #TODO this will fail if the time is not in nanoseconds generalize
    ns_per_minute = np.float64(60 * 1e9)
    dsw_dt_per_min = dsw_dt * ns_per_minute

    # Take absolute value and pad endpoints with NaNs
    dsw_dt_abs = np.abs(dsw_dt_per_min)
    dsw_dt_abs = dsw_dt_abs.reindex_like(global_irradiance)  # align with original time axis

    test_mask = (dsw_dt_abs <= max_dsw_dt) & dsw_dt_abs.notnull()
    test_mask.name = "test_change_with_time"


    test_mask.attrs = {}
    test_mask.attrs["info"] = "Mask based on change-with-time test on global shortwave (global variability test). See Long & Ackerman (2000) and subsequent iterations for details."
    test_mask.attrs["unit"] = "1", 
    test_mask.attrs["long_name"] = "clear sky classification mask",
    test_mask.attrs["flag_values"] = '0, 1',
    test_mask.attrs["flag_meanings"] = "0: fails change-with-time test (cloudy), 1: passes change-with-time test (possible clear-sky)"

    return test_mask

# TODO: find better name for this function
def diffuse_magnitude_test(diffuse_irradiance: xr.DataArray,
                            mu0: xr.DataArray,
                            *,
                            mu0_min: float,
                            diffuse_max_coeff: float,
                            diffuse_max_exp: float,
) -> xr.DataArray:
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
    # mu0 = self.mu0
    # mu0_min = self.get_attr('mu0_min')
    # diffuse_max_coeff = self.get_attr('diffuse_max_coeff')
    # diffuse_max_exp = self.get_attr('diffuse_max_exp')
    # sw_diffuse = self.dataset.diffuse_horizontal
    
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

    return test_mask

# TODO come up with a better name
def normalized_global_magnitude_test(
    sw_global: xr.DataArray,
    mu0: xr.DataArray,
    *,
    mu0_min: float,
    nsw_exp: float,
    nsw_min: float,
    nsw_max: float,
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

    valid = (mu0 >= mu0_min) & sw_global.notnull()

    # Avoid division by zero
    mu0_safe = mu0.where(mu0 > 0)
    nsw = sw_global / (mu0_safe ** nsw_exp)
    nsw = nsw.where(valid)

    test_mask = valid & (nsw >= nsw_min) & (nsw <= nsw_max)
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
    return test_mask

def fit_global_powerlaw_mu0(
    mu0: xr.DataArray,
    global_irradiance: xr.DataArray, 
    mask_clearsky: xr.DataArray,
    *,
    mu0_min: float,
    min_points: int = 100,) -> xr.DataArray | None:
    """
    Fit a simple power law for `global_irradiance`:
    global_irradiance = A * mu0^b
    using a linear regression in log space:
    log(global_irradiance) = log(A) + b * log(mu0)
    over points where `mask_clearsky` is True, mu0 >= mu0_min, and
    global_irradiance > 0.

    Returns
    -------
    xr.DataArray or None
        Labeled output with:
        - tcswd_a: coefficient
        - tcswd_b: exponent
        None if not enough valid points.
    """
    cond = (
        mask_clearsky
        & (mu0 >= mu0_min)
        & mu0.notnull()
        & global_irradiance.notnull()
        & (global_irradiance > 0)
    )
    cond_vals = cond.values
    if int(cond_vals.sum()) < min_points:
        return None

    mu0_sel = mu0.values[cond_vals]
    y_sel = global_irradiance.values[cond_vals]

    # Flatten and drop NaNs
    valid = np.isfinite(mu0_sel) & np.isfinite(y_sel) & (mu0_sel > 0) & (y_sel > 0)
    if valid.sum() < min_points:
        return None

    x = np.log(mu0_sel[valid])
    z = np.log(y_sel[valid])

    # Simple least-squares fit: z = log(A) + b * x
    b, logA = np.polyfit(x, z, 1)
    A = np.exp(logA)

    da =  xr.DataArray(
        np.array((A, b), dtype=np.float64),
        dims=("fit_params_tcswd",),
        coords={"fit_params_tcswd": np.array(("a", "b"), dtype=object)},
        name="global_powerlaw_mu0_fit",
    )
    da.attrs['info'] = 'Fit result for global_irradiance = a * mu0^b under clearsky conditions.'

    return da

def fit_diffuse_global_ratio_mu0_powerlaw(
    mu0: xr.DataArray,
    diffuse_irradiance: xr.DataArray,
    global_irradiance: xr.DataArray, 
    mask_clearsky: xr.DataArray,
    *,
    mu0_min: float,
    min_points: int = 100,) -> xr.DataArray | None:
    """
    Fit a simple power law for `global_irradiance`:
    global_irradiance = A * mu0^b
    using a linear regression in log space:
    log(global_irradiance) = log(A) + b * log(mu0)
    over points where `mask_clearsky` is True, mu0 >= mu0_min, and
    global_irradiance > 0.

    Returns
    -------
    xr.DataArray or None
        Labeled output with:
        - tcswd_a: coefficient
        - tcswd_b: exponent
        None if not enough valid points.
    """
    cond = (
        mask_clearsky
        & (mu0 >= mu0_min)
        & mu0.notnull()
        & diffuse_irradiance.notnull()
        & global_irradiance.notnull()
        & (global_irradiance > 0)
    )
    cond_vals = cond.values
    if int(cond_vals.sum()) < min_points:
        assert(False), 'should this really happen?'

    mu0_sel = mu0.values[cond_vals]
    y_sel = diffuse_irradiance.values[cond_vals] / global_irradiance.values[cond_vals]

    # Flatten and drop NaNs
    valid = np.isfinite(mu0_sel) & np.isfinite(y_sel) & (mu0_sel > 0) & (y_sel > 0)
    if valid.sum() < min_points:
        return None

    x = np.log(mu0_sel[valid])
    z = np.log(y_sel[valid])

    # Simple least-squares fit: z = log(A) + b * x
    b, logA = np.polyfit(x, z, 1)
    A = np.exp(logA)

    da =  xr.DataArray(
        np.array((A, b), dtype=np.float64),
        dims=("fit_params_dgrcswd",),
        coords={"fit_params_dgrcswd": np.array(("a", "b"), dtype=object)},
        name="fit_diffuse_global_ratio_mu0_powerlaw",
    )
    da.attrs['info'] = 'Fit result for (diffuce_irradiance / global_irradiance) = a * mu0^b under clearsky conditions. here b is ndr_exp'
    # infer a ndr_std_max
    # Normalized diffuse ratio
    ndr = pd.DataFrame(y_sel * (mu0_sel ** (-b)))
    ndr_std = ndr.rolling(11, center=True).std()
    da.attrs['ndr_std_max_estimated'] = float(ndr_std.mean().iloc[0])
    return da

def fit_diffuse_mu0_powerlaw(
    mu0: xr.DataArray,
    diffuse_irradiance: xr.DataArray,
    mask_clearsky: xr.DataArray,
    *,
    mu0_min: float,
    min_points: int = 100,) -> xr.DataArray | None:
    """
    Fit a simple power law for `global_irradiance`:
    global_irradiance = A * mu0^b
    using a linear regression in log space:
    log(global_irradiance) = log(A) + b * log(mu0)
    over points where `mask_clearsky` is True, mu0 >= mu0_min, and
    global_irradiance > 0.

    Returns
    -------
    xr.DataArray or None
        Labeled output with:
        - tcswd_a: coefficient
        - tcswd_b: exponent
        None if not enough valid points.
    """
    cond = (
        mask_clearsky
        & (mu0 >= mu0_min)
        & mu0.notnull()
        & diffuse_irradiance.notnull()
        # & (diffuse_irradiance > 0)
    )
    cond_vals = cond.values
    if int(cond_vals.sum()) < min_points:
        assert(False), 'should this really happen?'

    mu0_sel = mu0.values[cond_vals]
    y_sel = diffuse_irradiance.values[cond_vals]

    # Flatten and drop NaNs
    valid = np.isfinite(mu0_sel) & np.isfinite(y_sel) & (mu0_sel > 0) & (y_sel > 0)
    if valid.sum() < min_points:
        assert(False), 'we need to handle this in some way'
        return None

    x = np.log(mu0_sel[valid])
    z = np.log(y_sel[valid])

    # Simple least-squares fit: z = log(A) + b * x
    b, logA = np.polyfit(x, z, 1)
    A = np.exp(logA)

    da =  xr.DataArray(
        np.array((A, b), dtype=np.float64),
        dims=("fit_params_dcswd",),
        coords={"fit_params_dcswd": np.array(("a", "b"), dtype=object)},
        name="fit_diffuse_mu0_powerlaw",
    )
    da.attrs['info'] = '''Fit result for diffuce_irradiance = a * mu0^b under clearsky conditions.
    Here a is the infered value for diffuse_max_coeff minus a margin of 10-20%. 
    b is the infered value for diffuse_max_exp and should be close to 0.5'''
    # da.attrs['infert_diffuse_max_coeff'] = A * 1.2
    # da.attrs['infert_diffuse_max_exp'] = b

    return da