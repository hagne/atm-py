import numpy as np
import xarray as xr

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

def apply_tilt_correction(dataset: xr.Dataset, diffuse_global_ratio = None,
                            # offset: bool = True,
                         ) -> xr.Dataset:
    """Apply the Long et al. (2010) downwelling shortwave tilt correction.

    Required dataset variables use exact names. All irradiance units must be
    ``W/m^2`` and all angle units must be ``radian``. The platform convention is
    positive pitch nose-down and positive roll right-wing-down.

    
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
    mut = np.cos(incangle)
    Gt = dataset.global_horizontal

    ddt = dataset.diffuse_horizontal / Gt  # global_horizontal

    Diff = Gt*ddt
    G0 =  (-Gt*ddt*mu0 + Gt*ddt*mut + Gt*mu0)/mut
    N = (-Gt*ddt + Gt)/mut

    out = xr.Dataset()
    out["global_horizontal"] = G0
    out.global_horizontal.attrs.update(dataset.global_horizontal.attrs)
    out.global_horizontal.attrs["tilt_correction"] = "Long et al. 2010"
    out["diffuse_horizontal"] = Diff
    out.diffuse_horizontal.attrs.update(dataset.diffuse_horizontal.attrs)
    out["direct_normal"] = N
    out.direct_normal.attrs.update(dataset.direct_horizontal.attrs)
    out.direct_normal.attrs["tilt_correction"] = "Long et al. 2010" 


    out.attrs["tilt_corrected"] = "Long et al. 2010"
    return out
