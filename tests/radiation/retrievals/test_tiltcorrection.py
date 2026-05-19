import unittest
import sys
import types

import numpy as np
import xarray as xr

sys.modules.setdefault("timezonefinder", types.ModuleType("timezonefinder"))

from atmPy.radiation.retrievals import tiltcorrection
from atmPy.radiation.retrievals.broadband_shortwave_radiation import CombinedGlobalDiffuseDirect


class TiltCorrectionTest(unittest.TestCase):
    def _dataset(self):
        time = xr.DataArray(
            np.array(["2024-06-01T12:00:00"], dtype="datetime64[ns]"),
            dims="datetime",
            name="datetime",
        )
        zenith = xr.DataArray([np.deg2rad(60.0)], dims="datetime", coords={"datetime": time})
        pitch = xr.DataArray(
            [np.deg2rad(60.0) - np.arccos(0.4)],
            dims="datetime",
            coords={"datetime": time},
        )
        direct_normal = xr.DataArray([800.0], dims="datetime", coords={"datetime": time})
        diffuse = xr.DataArray([100.0], dims="datetime", coords={"datetime": time})
        mu_t = xr.DataArray([0.4], dims="datetime", coords={"datetime": time})
        global_tilted = mu_t * direct_normal + diffuse
        ds = xr.Dataset(
            {
                "global_horizontal": global_tilted,
                "diffuse_horizontal": diffuse,
                "direct_normal": direct_normal,
                "solar_zenith": zenith,
                "solar_zenith_geometric": zenith,
                "solar_elevation_geometric": np.pi / 2 - zenith,
                "solar_elevation": np.pi / 2 - zenith,
                "solar_azimuth": xr.DataArray([0.0], dims="datetime", coords={"datetime": time}),
                "solar_equation_of_time": xr.DataArray([0.0], dims="datetime", coords={"datetime": time}),
                "solar_airmass": xr.DataArray([2.0], dims="datetime", coords={"datetime": time}),
                "solar_airmass_absolute": xr.DataArray([2.0], dims="datetime", coords={"datetime": time}),
                "solar_sun_earth_distance": xr.DataArray([1.0], dims="datetime", coords={"datetime": time}),
                "platform_pitch": pitch,
                "platform_roll": xr.DataArray([0.0], dims="datetime", coords={"datetime": time}),
                "platform_heading": xr.DataArray([0.0], dims="datetime", coords={"datetime": time}),
            },
            coords={"datetime": time},
        )
        for name in ["global_horizontal", "diffuse_horizontal", "direct_normal"]:
            ds[name].attrs["units"] = "W/m^2"
        for name in ["solar_zenith", "solar_azimuth", "platform_pitch", "platform_roll", "platform_heading"]:
            ds[name].attrs["units"] = "radian"
        return ds

    def test_apply_tilt_correction_uses_long_equation(self):
        ds = self._dataset()
        corrected = tiltcorrection.apply_tilt_correction(ds, tilt=True)

        self.assertAlmostEqual(float(corrected.global_horizontal.isel(datetime=0)), 500.0)
        self.assertAlmostEqual(float(corrected.direct_horizontal.isel(datetime=0)), 400.0)
        self.assertAlmostEqual(float(corrected.direct_normal.isel(datetime=0)), 800.0)

    def test_apply_tilt_correction_requires_expected_units(self):
        ds = self._dataset()
        ds.platform_pitch.attrs["units"] = "degree"

        with self.assertRaisesRegex(ValueError, "platform_pitch units must be radian"):
            tiltcorrection.apply_tilt_correction(ds, tilt=True)

    def test_combined_wrapper_replaces_shared_dataset(self):
        rad = CombinedGlobalDiffuseDirect(self._dataset())
        corrected = rad.apply_tilt_correction(tilt=True)

        self.assertAlmostEqual(float(rad.dataset.global_horizontal.isel(datetime=0)), 500.0)
        self.assertAlmostEqual(float(corrected.dataset.global_horizontal.isel(datetime=0)), 500.0)


if __name__ == "__main__":
    unittest.main()
