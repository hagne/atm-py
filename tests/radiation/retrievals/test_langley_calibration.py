import unittest

import numpy as np
import pandas as pd
import xarray as xr

from atmPy.radiation.retrievals.langley_calibration import Langley, fit_langley


class FitLangleyTest(unittest.TestCase):
    def _langley(self, signal, airmass, datetime=None):
        channel = np.array([500, 870, 940])
        if datetime is None:
            datetime = pd.date_range('2026-05-02T12:00', periods=airmass.size, freq='1min')
        ds = xr.Dataset(
            {
                'direct_normal': (('airmass', 'channel'), signal),
                'datetime': ('airmass', datetime),
            },
            coords={'airmass': airmass, 'channel': channel},
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

    def test_supplied_weights_align_after_missing_values(self):
        airmass = np.array([2.2, 2.7, 3.1, 3.8, 4.6])
        signal = pd.Series([4.2, 3.9, np.nan, 3.1, 2.8], index=airmass)
        weights = pd.Series([4.0, 1.0, 3.0, 2.0, 5.0], index=airmass[[3, 0, 4, 1, 2]])

        result = fit_langley(signal, robust=False, weights=weights)

        retained = signal.dropna().sort_index(ascending=False)
        design = np.column_stack([np.ones(retained.size), retained.index.to_numpy(dtype=float)])
        aligned_weights = weights.reindex(retained.index).to_numpy(dtype=float)
        expected = np.linalg.solve(
            design.T @ (aligned_weights[:, np.newaxis] * design),
            design.T @ (aligned_weights * retained.to_numpy(dtype=float)),
        )
        self.assertAlmostEqual(result['res_series']['intercept'], expected[0])
        self.assertAlmostEqual(result['res_series']['slope'], expected[1])

    def test_expected_airmass_cdf_is_the_same_for_am_and_pm(self):
        airmass = np.array([2.1, 2.2, 2.4, 2.8, 3.6, 4.6])
        signal = 8.0 - airmass[:, np.newaxis] * np.array([0.1, 0.12, 0.2])
        datetime = pd.date_range('2026-05-02T12:00', periods=airmass.size, freq='1min')
        expected = np.linspace(0, 1, airmass.size)

        for timestamps in (datetime, datetime[::-1]):
            with self.subTest(datetime_order=timestamps[0]):
                langley = self._langley(signal, airmass, datetime=timestamps)
                np.testing.assert_allclose(langley.expected_airmass_cdf, expected)
                self.assertEqual(langley.expected_airmass_cdf.dims, ('airmass',))
                self.assertAlmostEqual(float(langley.airmass_weights.mean()), 1.0)

    def test_airmass_weights_are_inverse_expected_density(self):
        airmass = np.array([2.1, 2.3, 2.9, 4.1])
        signal = 8.0 - airmass[:, np.newaxis] * np.array([0.1, 0.12, 0.2])
        start = pd.Timestamp('2026-05-02T12:00')
        datetime = start + pd.to_timedelta([0, 1, 3, 6], unit='min')
        langley = self._langley(signal, airmass, datetime=datetime)

        expected_cdf = np.array([0.0, 1 / 6, 1 / 2, 1.0])
        expected_weights = np.array([1.2, 1.6, 2.16, 2.4])
        expected_weights /= expected_weights.mean()
        np.testing.assert_allclose(langley.expected_airmass_cdf, expected_cdf)
        np.testing.assert_allclose(langley.airmass_weights, expected_weights)

    def test_fixed_weights_and_coverage_after_masking(self):
        airmass = np.array([2.1, 2.2, 2.4, 2.8, 3.6, 4.6])
        signal = 8.0 - airmass[:, np.newaxis] * np.array([0.1, 0.12, 0.2])
        langley = self._langley(signal, airmass)
        full_weights = langley.airmass_weights.to_pandas()

        langley.langleys.loc[2.8] = np.nan
        langley.refresh()
        _ = langley.langley_fitres

        expected_weights = full_weights.reindex(langley.tp_langley.index)
        expected_weights /= expected_weights.mean()
        np.testing.assert_allclose(langley.tp_w, expected_weights)

        metrics = langley.airmass_coverage_metrics.sel(channel=500)
        self.assertEqual(metrics.retained_points.item(), 5)
        self.assertAlmostEqual(metrics.retained_fraction.item(), 5 / 6)
        self.assertAlmostEqual(metrics.cdf_span.item(), 1.0)
        self.assertAlmostEqual(metrics.maximum_cdf_gap.item(), 0.4)
        self.assertAlmostEqual(metrics.maximum_cdf_gap_ratio.item(), 2.0)
        retained_weights = full_weights.drop(index=2.8)
        effective_airmass = np.average(
            retained_weights.index.to_numpy(dtype=float),
            weights=retained_weights.to_numpy(dtype=float),
        )
        airmass_midpoint = (airmass.min() + airmass.max()) / 2
        airmass_half_range = (airmass.max() - airmass.min()) / 2
        self.assertAlmostEqual(metrics.effective_airmass.item(), effective_airmass)
        self.assertAlmostEqual(
            metrics.effective_airmass_offset.item(),
            (effective_airmass - airmass_midpoint) / airmass_half_range,
        )

    def test_coverage_gap_includes_missing_endpoints(self):
        airmass = np.array([2.1, 2.2, 2.4, 2.8, 3.6, 4.6])
        signal = 8.0 - airmass[:, np.newaxis] * np.array([0.1, 0.12, 0.2])
        langley = self._langley(signal, airmass)
        langley.langleys.loc[[2.1, 2.2], 500] = np.nan

        metrics = langley.airmass_coverage_metrics.sel(channel=500)
        self.assertAlmostEqual(metrics.maximum_cdf_gap.item(), 0.4)
        self.assertAlmostEqual(metrics.maximum_cdf_gap_ratio.item(), 2.0)
        self.assertGreater(metrics.effective_airmass_offset.item(), 0)

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
