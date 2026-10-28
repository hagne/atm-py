import unittest

import numpy as np
import pandas as pd
import xarray as xr

from atmPy.radiation.retrievals.spectral_irradiance import SolarIrradiation


class ApplyCalibrationLangleyTest(unittest.TestCase):
    def test_adds_relative_irradiance_uncertainty(self):
        datetime = pd.date_range('2026-01-01', periods=1)
        channel = [415, 500]
        ds = xr.Dataset(
            {
                'global_horizontal': (('datetime', 'channel'), [[50.0, 100.0]]),
                'diffuse_horizontal': (('datetime', 'channel'), [[10.0, 20.0]]),
                'direct_normal': (('datetime', 'channel'), [[40.0, 80.0]]),
                'toa_spectral_irradiance': ('channel', [1.0, 1.0]),
                'channel_wavelength': ('channel', channel),
            },
            coords={'datetime': datetime, 'channel': channel},
            attrs={
                'site': 'test',
                'site_name': 'Test site',
                'site_elevation': 0.0,
                'site_latitude': 0.0,
                'site_longitude': 0.0,
            },
        )
        sun_position_variables = [
            'zenith',
            'zenith_geometric',
            'elevation_geometric',
            'elevation',
            'azimuth',
            'equation_of_time',
            'airmass',
            'airmass_absolute',
            'sun_earth_distance',
        ]
        for variable in sun_position_variables:
            ds[variable] = ('datetime', [1.0])

        calibration = xr.Dataset(
            {
                'V0': ('channel', [100.0, 200.0]),
                'V0_uncertainty': ('channel', [2.0, 10.0]),
            },
            coords={'channel': channel},
        )

        calibrated = SolarIrradiation(ds).apply_calibration_langley(calibration).dataset

        np.testing.assert_allclose(calibrated.irradiance_uncertainty, [2.0, 5.0])
        self.assertEqual(calibrated.irradiance_uncertainty.dims, ('channel',))
        self.assertEqual(calibrated.irradiance_uncertainty.attrs['units'], '%')
        self.assertIn('95% confidence', calibrated.irradiance_uncertainty.attrs['description'])
        np.testing.assert_allclose(calibrated.aod_uncertainty, [0.02, 0.05])
        self.assertEqual(calibrated.aod_uncertainty.dims, ('channel',))
        self.assertEqual(calibrated.aod_uncertainty.attrs['units'], '1')
        self.assertEqual(calibrated.aod_uncertainty.attrs['reference_air_mass_factor'], 1.0)
        self.assertIn('aod_uncertainty / airmass', calibrated.aod_uncertainty.attrs['description'])


if __name__ == '__main__':
    unittest.main()
