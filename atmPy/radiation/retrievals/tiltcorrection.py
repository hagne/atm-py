import numpy as np
import xarray as xr
from scipy.signal import lfilter


def effective_mu(mut: xr.DataArray, tau: float) -> xr.DataArray:
    def _first_order_lowpass_1d(x, time, dt, tau):
        x = np.asarray(x, dtype=float)
        out = np.full_like(x, np.nan)
        finite = np.isfinite(x)
        if not finite.any():
            return out

        alpha = dt / (tau + dt)
        gap = np.diff(time) / np.timedelta64(1, "s")
        is_new_segment = finite.copy()
        is_new_segment[1:] &= (~finite[:-1]) | (gap > 3 * dt)
        is_end_segment = finite.copy()
        is_end_segment[:-1] &= (~finite[1:]) | (gap > 3 * dt)

        starts = np.flatnonzero(is_new_segment)
        stops = np.flatnonzero(is_end_segment) + 1
        for start, stop in zip(starts, stops):
            segment = x[start:stop]
            out[start:stop] = lfilter(
                [alpha],
                [1, -(1 - alpha)],
                segment,
                zi=[(1 - alpha) * segment[0]],
            )[0]
        return out

    dt = mut.datetime.diff("datetime").median().values / np.timedelta64(1, "s")
    # return mut
    mut = mut.clip(min=0)

    mu_eff = xr.apply_ufunc(
        _first_order_lowpass_1d,
        mut,
        mut.datetime,
        input_core_dims=[["datetime"], ["datetime"]],
        output_core_dims=[["datetime"]],
        kwargs={"dt": float(dt), "tau": tau},
        # vectorize=True,
        # dask="parallelized",
        # output_dtypes=[float],
    )
    # return mu_eff
    # mu_eff = mu_eff.transpose(*mut.dims)
    mu_eff.name = "mu_eff"
    mu_eff.attrs = {
        "long_name": "effective cosine of solar incidence angle",
        "units": "1",
        "description": (
            "Cosine of solar incidence angle filtered with a causal "
            "first-order low-pass response."
        ),
        "tau": tau,
        "tau_units": "s",
    }

    return mu_eff

def solar_incidence_angle(zenith, azimuth, pitch, roll, heading):
    """
    Calculate solar incidence angle on a tilted hemispheric detector.

    All angles are in radians.

    Conventions:
    - azimuth: clockwise from north
    - heading: clockwise from north
    - pitch: nose-up positive
    - roll: right-side-down positive
    """
    if isinstance(pitch, xr.DataArray):
        assert('positive' in pitch.attrs), 'pitch variable must have a "positive" attribute indicating the positive direction of pitch. Positive pitch should be nose-up.'
        assert(pitch.attrs['positive'] == 'nose-up'), f'Positive pitch should be nose-up. Current positive direction: {pitch.attrs["positive"]}'
    if isinstance(roll, xr.DataArray):
        assert('positive' in roll.attrs), 'roll variable must have a "positive" attribute indicating the positive direction of roll. Positive roll should be right-side-down.'
        assert(roll.attrs['positive'] == 'right-side-down'), f'Positive roll should be right-side-down. Current positive direction: {roll.attrs["positive"]}'

    # Sun vector in local ENU coordinates: x=east, y=north, z=up
    sx = np.sin(zenith) * np.sin(azimuth)
    sy = np.sin(zenith) * np.cos(azimuth)
    sz = np.cos(zenith)

    # Platform axes in ENU coordinates
    # Heading unit vector: platform forward direction
    fx = np.sin(heading)
    fy = np.cos(heading)
    fz = 0

    # Rightward unit vector
    rx = np.cos(heading)
    ry = -np.sin(heading)
    rz = 0

    # Up vector
    ux = 0
    uy = 0
    uz = 1

    # Detector normal after pitch and roll
    nx = (
        ux * np.cos(pitch) * np.cos(roll)
        - fx * np.sin(pitch)
        + rx * np.sin(roll)
    )
    ny = (
        uy * np.cos(pitch) * np.cos(roll)
        - fy * np.sin(pitch)
        + ry * np.sin(roll)
    )
    nz = (
        uz * np.cos(pitch) * np.cos(roll)
        - fz * np.sin(pitch)
        + rz * np.sin(roll)
    )

    # Normalize detector normal
    norm = np.sqrt(nx**2 + ny**2 + nz**2)
    nx /= norm
    ny /= norm
    nz /= norm

    cos_incidence = sx * nx + sy * ny + sz * nz
    cos_incidence = np.clip(cos_incidence, -1, 1)

    return np.arccos(cos_incidence)

def apply_tilt_correction(dataset: xr.Dataset, 
                          diffuse_global_ratio = None,
                          sensor_response_time = None,
                            # offset: bool = True,
                         ) -> xr.Dataset:
    """Apply the Long et al. (2010) downwelling shortwave tilt correction.

    Required dataset variables use exact names. All irradiance units must be
    ``W/m^2`` and all angle units must be ``radian``. The platform convention is
    positive pitch nose-down and positive roll right-wing-down.
    Parameters
    ----------
    dataset : xr.Dataset
    sensor_response_time : float, optional
        If provided, applies an additional correction to account for sensor response time. Response time should be the 95% response time of the sensor in seconds.
    """
        # if not isinstance(offset, bool):
    #     raise TypeError(f"offset must be a boolean, it is {type(offset)}")

    required_units = {
        "global_horizontal": "W/m^2",
        "diffuse_horizontal": "W/m^2",
        # "direct_horizontal": "W/m^2",
        "solar_zenith": "radian",
        "solar_azimuth": "radian",
        "platform_pitch": "radian",
        "platform_roll": "radian",
        "platform_heading": "radian",
            }
    missing_vars = [name for name in required_units if name not in dataset]
    if missing_vars:
        raise ValueError(f"The following required variables are missing from the dataset: {', '.join(missing_vars)}")
    wrong_units = [(name, dataset[name].attrs.get("units"), required_units[name]) for name in required_units if name in dataset and dataset[name].attrs.get("units") != required_units[name]]
    if wrong_units:
        raise ValueError(f"The following variables have incorrect units: {', '.join(f'{name}: {units} (expected: {expected_units})' for name, units, expected_units in wrong_units)}")
    # for name, units in required_units.items():
    #     if name not in dataset:
    #         raise ValueError(f"{name} variable is missing from dataset")
    #     if dataset[name].attrs.get("units") != units:
    #         raise ValueError(f"{name} units must be {units}. Current units: {dataset[name].attrs.get('units')}")

    # out = dataset.copy()
    mu0 = np.cos(dataset.solar_zenith)
    incangle = solar_incidence_angle(dataset.solar_zenith, dataset.solar_azimuth, dataset.platform_pitch, dataset.platform_roll, dataset.platform_heading)
    dataset["solar_incidence_angle"] = incangle
    dataset.solar_incidence_angle.attrs["units"] = "radian"
    dataset.solar_incidence_angle.attrs["long_name"] = "Solar incidence angle on tilted irradiance sensor"
    mut = np.cos(incangle)

    if sensor_response_time is not None:
        tau= -sensor_response_time/np.log(0.05)
        mut = effective_mu(mut, tau)

    Gt = dataset.global_horizontal

    ddt = dataset.diffuse_horizontal / Gt  # global_horizontal

    Diff = Gt*ddt
    G0 =  (-Gt*ddt*mu0 + Gt*ddt*mut + Gt*mu0)/mut
    N = (-Gt*ddt + Gt)/mut

    out = xr.Dataset()
    out["global_horizontal"] = G0
    out.global_horizontal.attrs.update(dataset.global_horizontal.attrs)
    tilt_correction_description = "Long et al. 2010"
    if sensor_response_time is not None:
        tilt_correction_description += f" with sensor response time correction (tau={tau:.1f} s)"
    out.global_horizontal.attrs["tilt_correction"] = tilt_correction_description
    out["diffuse_horizontal"] = Diff
    out.diffuse_horizontal.attrs.update(dataset.diffuse_horizontal.attrs)
    out["direct_normal"] = N
    # out.direct_normal.attrs.update(dataset.direct_horizontal.attrs)
    out.direct_normal.attrs["tilt_correction"] = tilt_correction_description


    out.attrs["tilt_corrected"] = "Long et al. 2010"

    out['mut'] = mut    
    out['mu0'] = mu0

    return out

from scipy.optimize import minimize

def minthis(ds, dr, dp, dh, sensor_response_time=0.2, poly_order=3, verbose=False):
    dst = ds.copy()
    dst['platform_roll'] = dst.platform_roll + dr
    dst['platform_pitch'] = dst.platform_pitch + dp
    dst['platform_heading'] = dst.platform_heading + dh

    # rin = atmbsr.CombinedGlobalDiffuseDirect(dst, verbose=verbose)
    out = apply_tilt_correction(dst, sensor_response_time=sensor_response_time)

    y = out.direct_normal
    dim = y.dims[0]
    x = y[dim]
    zenith_min_time = ds.solar_zenith.idxmin(dim=dim)
    if np.issubdtype(x.dtype, np.datetime64):
        xfit = (x - zenith_min_time) / np.timedelta64(1, 's')
    else:
        xfit = x - zenith_min_time

    yfit = y.values
    xfit = np.asarray(xfit.values, dtype=float)
    mask = np.isfinite(xfit) & np.isfinite(yfit)
    even_order = poly_order if poly_order % 2 == 0 else poly_order - 1
    powers = np.arange(0, even_order + 1, 2)
    design_matrix = np.column_stack([xfit[mask]**power for power in powers])
    coef, *_ = np.linalg.lstsq(design_matrix, yfit[mask], rcond=None)
    residual = yfit[mask] - design_matrix @ coef
    return float(np.std(residual))

def find_offset(ds, dr=None, dp=None, dh=None):
    """ Very experimental! Find installation offsets that minimize direct-normal irradiance variability.
    Note, this only works for moving platforms and mostly clearsky days. You will actuall
    need roll, pitch and heading to permanently change to get a result. 

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with platform_roll, platform_pitch, platform_heading,
        global_horizontal, and diffuse_horizontal variables.
    dr, dp, dh : float or None, optional
        Roll, pitch, and heading offsets in radians. If an offset is None, it is
        optimized. If an offset is provided, that value is held fixed and is not
        optimized.

    Returns
    -------
    pandas.Series
        Best offsets in radians and degrees, the minimized direct-normal
        standard deviation, and scipy optimization status fields.
    """
    # Search offsets in radians, but define/report them in degrees for readability.
    offset_names = ['dr', 'dp', 'dh']
    provided_offsets = {'dr': dr, 'dp': dp, 'dh': dh}
    free_offsets = [name for name in offset_names if provided_offsets[name] is None]
    bounds_deg = {name: (-15, 15) for name in offset_names}
    bounds = [tuple(np.deg2rad(v) for v in bounds_deg[name]) for name in free_offsets]
    optimization_trace = []

    def build_offsets(x):
        offsets = provided_offsets.copy()
        for name, value in zip(free_offsets, x):
            offsets[name] = value
        return offsets

    def objective(x):
        offsets = build_offsets(x)
        value = minthis(
            ds,
            offsets['dr'],
            offsets['dp'],
            offsets['dh'],
            sensor_response_time=0.2,
            verbose=False,
        )
        optimization_trace.append({
            'dr': offsets['dr'],
            'dp': offsets['dp'],
            'dh': offsets['dh'],
            'value': value,
        })
        return value

    if free_offsets:
        opt = minimize(
            objective,
            x0=np.zeros(len(free_offsets)),
            method='Powell',
            bounds=bounds,
            options={'xtol': 1e-4, 'ftol': 1e-4, 'disp': True},
        )
        best_offsets = build_offsets(opt.x)
        direct_normal_std = opt.fun
        success = opt.success
        message = opt.message
    else:
        best_offsets = provided_offsets
        direct_normal_std = objective([])
        success = True
        message = 'No offsets optimized; used provided values.'

    best = pd.Series({
        'dr_rad': best_offsets['dr'],
        'dp_rad': best_offsets['dp'],
        'dh_rad': best_offsets['dh'],
        'dr_deg': np.rad2deg(best_offsets['dr']),
        'dp_deg': np.rad2deg(best_offsets['dp']),
        'dh_deg': np.rad2deg(best_offsets['dh']),
        'direct_normal_std': direct_normal_std,
        'success': success,
        'message': message,
    })

    return best
