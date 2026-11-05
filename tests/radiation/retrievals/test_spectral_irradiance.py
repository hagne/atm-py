import unittest
import warnings

import numpy as np
import pandas as pd
import xarray as xr

from atmPy.radiation.retrievals.spectral_irradiance import (
    DirectNormalIrradiation,
    SolarIrradiation,
    merge_by_datetime,
)


class MergeByDatetimeTest(unittest.TestCase):
    def setUp(self):
        self.target = xr.Dataset(
            coords={
                'datetime': np.array(
                    ['2026-01-01T00:00', '2026-01-01T00:30', '2026-01-01T01:00'],
                    dtype='datetime64[ns]',
                )
            }
        )

    def test_interpolates_and_adds_exact_variables(self):
        data = xr.Dataset(
            data_vars={
                'pressure': ('datetime', [1000.0, 990.0], {'units': 'hPa'}),
                'temperature': ('datetime', [0.0, 2.0], {'units': 'DegC'}),
            },
            coords={
                'datetime': np.array(
                    ['2026-01-01T00:00', '2026-01-01T01:00'],
                    dtype='datetime64[s]',
                )
            },
        )

        interpolated = merge_by_datetime(
            self.target,
            data,
            {'pressure': 'hPa', 'temperature': 'DegC'},
        )

        np.testing.assert_allclose(self.target.pressure, [1000.0, 995.0, 990.0])
        np.testing.assert_allclose(self.target.temperature, [0.0, 1.0, 2.0])
        self.assertEqual(interpolated.datetime.dtype, self.target.datetime.dtype)
        self.assertEqual(self.target.pressure.attrs['units'], 'hPa')

    def test_requires_exact_variable_and_unit(self):
        wrong_variable = xr.Dataset(
            {'air_pressure': ('datetime', [1000.0], {'units': 'hPa'})},
            coords={'datetime': self.target.datetime[:1]},
        )
        with self.assertRaisesRegex(KeyError, 'pressure'):
            merge_by_datetime(self.target, wrong_variable, {'pressure': 'hPa'})

        wrong_unit = wrong_variable.rename({'air_pressure': 'pressure'})
        wrong_unit.pressure.attrs['units'] = 'Pa'
        with self.assertRaisesRegex(ValueError, 'hPa'):
            merge_by_datetime(self.target, wrong_unit, {'pressure': 'hPa'})

    def test_rejects_data_without_temporal_overlap(self):
        data = xr.Dataset(
            {'pressure': ('datetime', [1000.0], {'units': 'hPa'})},
            coords={'datetime': [np.datetime64('2026-01-02')]},
        )

        with self.assertRaisesRegex(ValueError, 'no temporal overlap'):
            merge_by_datetime(self.target, data, {'pressure': 'hPa'})


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


class AngstromExponentTest(unittest.TestCase):
    def test_masks_infinite_aod_and_preserves_finite_result(self):
        datetime = pd.date_range('2026-01-01', periods=2)
        channel = np.array([415, 500, 670, 870, 1625])
        channel_wavelength = channel.astype(float)
        expected_alpha = 1.3
        expected_beta = 0.12
        finite_aod = expected_beta * (channel_wavelength * 1e-3) ** -expected_alpha
        aod = np.stack([finite_aod, finite_aod])
        aod[1, channel == 870] = np.inf
        ds = xr.Dataset(
            data_vars={
                'direct_normal': (('datetime', 'channel'), np.ones_like(aod)),
                'channel_wavelength': ('channel', channel_wavelength),
                'aod': (('datetime', 'channel'), aod),
                'aod_uncertainty': (
                    ('datetime', 'channel'),
                    np.full_like(aod, 0.01),
                ),
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

        direct_normal = DirectNormalIrradiation(ds)
        with warnings.catch_warnings():
            warnings.simplefilter('error', RuntimeWarning)
            angstrom_exponent = direct_normal.angstrom_exponent

        pair = angstrom_exponent.sel(angstrom_exponent_pair='415_870')
        self.assertAlmostEqual(pair.isel(datetime=0).item(), expected_alpha)
        self.assertTrue(np.isnan(pair.isel(datetime=1).item()))
        beta = direct_normal.angstrom_turbidity_coefficient.sel(
            angstrom_exponent_pair='415_870'
        )
        self.assertAlmostEqual(beta.isel(datetime=0).item(), expected_beta)


if __name__ == '__main__':
    unittest.main()
