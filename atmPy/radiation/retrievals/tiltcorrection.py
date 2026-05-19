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
    assert('positive' in pitch.attrs), 'pitch variable must have a "positive" attribute indicating the positive direction of pitch. Positive pitch should be nose-up.'
    assert(pitch.attrs['positive'] == 'nose-up'), f'Positive pitch should be nose-up. Current positive direction: {pitch.attrs["positive"]}'
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

def apply_tilt_correction(
    dataset: xr.Dataset,
    # offset: bool = True,
) -> xr.Dataset:
    """Apply the Long et al. (2010) downwelling shortwave tilt correction.

    Required dataset variables use exact names. All irradiance units must be
    ``W/m^2`` and all angle units must be ``radian``. The platform convention is
    positive pitch nose-down and positive roll right-wing-down.

    Let $G_t$ denote the measured total irradiance, $G_0$ the corrected total irradiance, 
    $N$ the direct normal component, $\mathrm{Diff}$ the diffuse component, $\mu_t$ the 
    cosine of the tilted solar zenith angle, $\mu_0$ the cosine of the horizontal solar 
    zenith angle, and $\mathrm{ddt}$ the diffuse-global ratio of tilted observations.

    Solving the system

    $$
    \begin{aligned}
    G_t &= N\mu_t + \mathrm{Diff} \\
    G_0 &= N\mu_0 + \mathrm{Diff} \\
    \mathrm{ddt} &= \frac{\mathrm{Diff}}{N}\mu_t
    \end{aligned}
    $$

    for the unknowns $G_0$, $N$, and $\mathrm{Diff}$, in terms of the remaining variables $G_t$, $\mu_t$, $\mu_0$, and $\mathrm{ddt}$, yields

    $$
    (G_0, N, \mathrm{Diff})
    =
    \left(
    \frac{G_t(\mathrm{ddt} + \mu_0\mu_t)}{\mathrm{ddt} + \mu_t^2},
    \frac{G_t\mu_t}{\mathrm{ddt} + \mu_t^2},
    \frac{G_t\mathrm{ddt}}{\mathrm{ddt} + \mu_t^2}
    \right).
    $$
    """
    print('bla', flush=True)
    # if not isinstance(offset, bool):
    #     raise TypeError(f"offset must be a boolean, it is {type(offset)}")

    required_units = {
        "global_horizontal": "W/m^2",
        "diffuse_horizontal": "W/m^2",
        "direct_horizontal": "W/m^2",
        "solar_zenith": "radian",
        "solar_azimuth": "radian",
        "platform_pitch": "radian",
        "platform_roll": "radian",
        "platform_heading": "radian",
            }

    for name, units in required_units.items():
        if name not in dataset:
            raise ValueError(f"{name} variable is missing from dataset")
        if dataset[name].attrs.get("units") != units:
            raise ValueError(f"{name} units must be {units}. Current units: {dataset[name].attrs.get('units')}")

    # out = dataset.copy()
    mu0 = np.cos(dataset.solar_zenith)
    incangle = solar_incidence_angle(dataset.solar_zenith, dataset.solar_azimuth, dataset.platform_pitch, dataset.platform_roll, dataset.platform_heading)
    mut = np.cos(incangle)
    ddt = dataset.diffuse_horizontal / dataset.global_horizontal
    Gt = dataset.global_horizontal

    Diff = Gt*ddt/(ddt + mut**2)
    G0 =  Gt*(ddt + mu0*mut)/(ddt + mut**2)
    N =   Gt*mut/(ddt + mut**2)

    out = xr.Dataset()
    out["global_horizontal"] = G0
    out.global_horizontal.attrs.update(dataset.global_horizontal.attrs)
    out.global_horizontal.attrs["tilt_correction"] = "Long et al. 2010"
    out["diffuse_horizontal"] = Diff
    out.diffuse_horizontal.attrs.update(dataset.diffuse_horizontal.attrs)
    out["direct_normal"] = N
    out.direct_normal.attrs.update(dataset.direct_horizontal.attrs)
    out.direct_normal.attrs["tilt_correction"] = "Long et al. 2010" 


    # sun_relative_azimuth = out.solar_azimuth - out.platform_heading
    # mu_tilt = (
    #     np.cos(out.platform_pitch) * np.cos(out.platform_roll) * mu0
    #     + np.sin(out.platform_pitch)
    #     * np.cos(out.platform_roll)
    #     * np.sin(out.solar_zenith)
    #     * np.cos(sun_relative_azimuth)
    #     - np.sin(out.platform_roll)
    #     * np.sin(out.solar_zenith)
    #     * np.sin(sun_relative_azimuth)
    # )

    # numerator = mu0 * out.direct_normal + out.diffuse_horizontal
    # denominator = mu_tilt * out.direct_normal + out.diffuse_horizontal
    # correction_factor = (numerator / denominator).where(denominator != 0)

    # if "global_horizontal_uncorrected_tilt" not in out:
    #     out["global_horizontal_uncorrected_tilt"] = out.global_horizontal
    # if "direct_normal_uncorrected_tilt" not in out:
    #     out["direct_normal_uncorrected_tilt"] = out.direct_normal
    # if "direct_horizontal" in out and "direct_horizontal_uncorrected_tilt" not in out:
    #     out["direct_horizontal_uncorrected_tilt"] = out.direct_horizontal

    # out["tilt_correction_factor"] = correction_factor
    # out["mu_tilt"] = mu_tilt
    # out["global_horizontal"] = out.global_horizontal * correction_factor
    # out["direct_horizontal"] = out.global_horizontal - out.diffuse_horizontal
    # out["direct_normal"] = out.direct_horizontal / mu0

    # out.global_horizontal.attrs.update(dataset.global_horizontal.attrs)
    # out.direct_horizontal.attrs["units"] = out.global_horizontal.attrs["units"]
    # out.direct_normal.attrs.update(dataset.direct_normal.attrs)
    # out.tilt_correction_factor.attrs.update(
    #     {"long_name": "Long et al. 2010 tilt correction factor"}
    # )
    # out.mu_tilt.attrs.update(
    #     {"long_name": "cosine of angle between sensor normal and direct solar beam"}
    # )
    # out.attrs["tilt_corrected"] = "Long et al. 2010"
    return out
