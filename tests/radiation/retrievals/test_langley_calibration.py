import unittest

import numpy as np
import pandas as pd
import xarray as xr

from atmPy.radiation.retrievals.langley_calibration import Langley, fit_langley


class FitLangleyTest(unittest.TestCase):
    def _langley(self, signal, airmass):
        wavelength = np.array([500, 870, 940])
        ds = xr.Dataset(
            {'direct_normal': (('airmass', 'wavelength'), signal)},
            coords={'airmass': airmass, 'wavelength': wavelength},
        )
        return Langley(None, ds, airmass_limits=(2, 5))

    def test_robust_weighted_fit(self):
        airmass = np.array([5.0, 4.0, 3.0, 2.9, 2.8, 2.7, 2.6, 2.5])
        signal = 4.0 - 0.2 * airmass
        signal += np.array([0.01, -0.01, 0.015, 0.0, -0.01, 0.01, -0.005, 0.0])
        signal[4] += 1.0

        result = fit_langley(pd.Series(signal, index=airmass), robust=True, weighted=True)

        self.assertAlmostEqual(result['res_series']['slope'], -0.20334998395164727)
        self.assertAlmostEqual(result['res_series']['intercept'], 4.0151440866569414)

    def test_channel_noise_floor_controls_individual_removal(self):
        airmass = np.linspace(2.3, 4.6, 15)
        signal = 8.0 - airmass[:, np.newaxis] * np.array([0.1, 0.12, 0.2])
        signal[7, :2] += 0.15

        result = self._langley(signal, airmass).clean(
            threshold={870: 3, 500: 3},
            use_channels=[500, 870],
            clean_channels_individually=True,
            expected_instrument_noise={870: 0.1, 500: 0.01},
        )

        cleaned = result['langley'].langleys.loc[airmass[7]]
        self.assertTrue(np.isnan(cleaned[500]))
        self.assertTrue(np.isfinite(cleaned[[870, 940]]).all())

    def test_common_removal_requires_channel_agreement(self):
        airmass = np.linspace(2.3, 4.6, 21)
        signal = 8.0 - airmass[:, np.newaxis] * np.array([0.1, 0.12, 0.2])
        signal[10, :2] += 0.15
        signal[5, 0] += 0.15

        result = self._langley(signal, airmass).clean(
            threshold=3,
            use_channels=[500, 870],
        )

        self.assertNotIn(airmass[10], result['langley'].langleys.index)
        self.assertIn(airmass[5], result['langley'].langleys.index)
        self.assertTrue(result['removed_points'].loc[airmass[10]])
        self.assertFalse(result['removed_points'].loc[airmass[5]])


if __name__ == '__main__':
    unittest.main()
