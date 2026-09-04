import sys
import types
import unittest

import numpy as np
import xarray as xr


sys.modules.setdefault("timezonefinder", types.ModuleType("timezonefinder"))

from atmPy.radiation.radflux.lab import RadFlux


class RadFluxClearSkyParameterNamesTest(unittest.TestCase):
    def setUp(self):
        ds = xr.Dataset(
            coords={
                "datetime": np.array(
                    ["2020-01-01T00:00"], dtype="datetime64[m]"
                )
            }
        )
        self.radflux = RadFlux(ds)

    def test_default_parameters_use_descriptive_public_names(self):
        self.radflux.clearsky_parameters = "default"

        params = self.radflux.clearsky_parameters

        self.assertEqual(params["maximum_solar_zenith_angle"], 90.0)
        self.assertEqual(
            params["normalized_total_shortwave_power_exponent"], 1.18
        )
        self.assertEqual(
            params[
                "normalized_diffuse_ratio_standard_deviation_window_minutes"
            ],
            11,
        )
        self.assertNotIn("mu0_min", params)
        self.assertNotIn("nsw_exp", params)
        self.assertNotIn("ndr_window", self.radflux.dataset.attrs)

    def test_abbreviated_input_is_normalized_at_public_boundary(self):
        self.radflux.clearsky_parameters = {
            "mu0_min": 0.5,
            "nsw_exp": 1.3,
            "ndr_window": 9.8,
        }

        params = self.radflux.clearsky_parameters

        self.assertAlmostEqual(params["maximum_solar_zenith_angle"], 60.0)
        self.assertEqual(
            params["normalized_total_shortwave_power_exponent"], 1.3
        )
        self.assertEqual(
            params[
                "normalized_diffuse_ratio_standard_deviation_window_minutes"
            ],
            9,
        )
        self.assertAlmostEqual(self.radflux.get_attr("mu0_min"), 0.5)
        self.assertEqual(
            self.radflux.get_attr(
                "normalized_total_shortwave_low_sun_cosine_boundary"
            ),
            0.2,
        )
        self.assertEqual(
            self.radflux.get_attr(
                "normalized_total_shortwave_power_exponent"
            ),
            1.3,
        )

    def test_existing_abbreviated_attributes_remain_readable(self):
        self.radflux.dataset.attrs.update(
            {
                "mu0_min": 0.25,
                "nsw_coeff": 1050.0,
                "nsw_exp": 1.2,
            }
        )

        params = self.radflux.clearsky_parameters

        self.assertAlmostEqual(
            params["maximum_solar_zenith_angle"],
            np.rad2deg(np.arccos(0.25)),
        )
        self.assertEqual(
            params["normalized_total_shortwave_power_coefficient"], 1050.0
        )
        self.assertEqual(
            self.radflux.get_attr(
                "normalized_total_shortwave_power_exponent"
            ),
            1.2,
        )

    def test_magnitude_tests_use_descriptive_dataset_attributes(self):
        mu0 = xr.DataArray(
            [0.1, 0.5, 0.9],
            dims="datetime",
            coords={
                "datetime": np.array(
                    [
                        "2020-01-01T00:00",
                        "2020-01-01T00:01",
                        "2020-01-01T00:02",
                    ],
                    dtype="datetime64[m]",
                )
            },
        )
        ds = xr.Dataset(
            {
                "mu0": mu0,
                "global_horizontal": 1000 * mu0**1.18,
                "diffuse_horizontal": 50 * mu0**0.5,
            }
        )
        radflux = RadFlux(ds)
        radflux.clearsky_parameters = "default"

        np.testing.assert_array_equal(
            radflux.mask_normalized_global_magnitude,
            [True, True, True],
        )
        np.testing.assert_array_equal(
            radflux.mask_diffuse_magnitude,
            [True, True, True],
        )


if __name__ == "__main__":
    unittest.main()
