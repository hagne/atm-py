import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from atmPy.radiation.retrievals.langley_calibration import Langley
from atmPy.radiation.retrievals.langley_db import LangleyResultDatabase


class LangleyResultDatabaseTest(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.database_path = Path(self.temporary_directory.name) / 'langley.sqlite'
        self.database = LangleyResultDatabase(
            self.database_path,
            create_if_not_exist=True,
            version='1.2.3',
        )

        airmass = np.array([2.1, 2.3, 2.7, 3.2, 3.8, 4.6])
        channels = np.array([500, 870, 940])
        signal = 8.0 - airmass[:, np.newaxis] * np.array([0.1, 0.12, 0.2])
        signal[2, 500 == channels] = np.nan
        dataset = xr.Dataset(
            {
                'direct_normal': (('airmass', 'channel'), signal),
                'datetime': (
                    'airmass',
                    pd.date_range('2026-05-02T12:00', periods=airmass.size, freq='1min'),
                ),
            },
            coords={'airmass': airmass, 'channel': channels},
        )
        self.langley = Langley(
            None,
            dataset,
            airmass_limits=(2, 5),
            when='am',
        )

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_write_langley_stores_fit_coverage_and_processing_results(self):
        self.database.write_langley(
            '2026-05-02',
            '/data/input.nc',
            self.langley,
            processing_date='2026-05-03T01:02:03',
            processing_server='test-host',
        )

        stored = self.database.dump_langley_results()
        self.assertEqual(stored.shape[0], 3)
        result = stored.loc[('2026-05-02', 'am', 500)]
        fitresults = self.langley.langley_fitres.loc[500]
        coverage = self.langley.airmass_coverage_metrics.sel(channel=500)

        self.assertAlmostEqual(result.intercept, fitresults.intercept)
        self.assertAlmostEqual(result.slope_stderr, fitresults.slope_stderr)
        self.assertEqual(result.retained_points, coverage.retained_points.item())
        self.assertAlmostEqual(result.maximum_cdf_gap, coverage.maximum_cdf_gap.item())
        self.assertTrue({
            'expected_points',
            'airmass_midpoint',
            'airmass_half_range',
        }.isdisjoint(stored.columns))
        self.assertAlmostEqual(
            result.residual_correlation_determinant,
            self.langley.langley_residual_correlation_prop['determinant'],
        )
        self.assertEqual(result.input_file, '/data/input.nc')
        self.assertEqual(result.processed_at, '2026-05-03T01:02:03')
        self.assertEqual(result.process_version, '1.2.3')
        self.assertEqual(result.processing_server, 'test-host')

    def test_write_langley_results_updates_existing_channel_rows(self):
        write_kwargs = {
            'date': '2026-05-02',
            'path2file': '/data/input.nc',
            'when': 'pm',
            'langley_fitresults': self.langley.langley_fitres.copy(),
            'airmass_coverage_metrics': self.langley.airmass_coverage_metrics,
            'processing_date': '2026-05-03T01:02:03',
            'processing_server': 'test-host',
        }
        self.database.write_langley_results(**write_kwargs)
        write_kwargs['langley_fitresults'].loc[500, 'intercept'] = 9.25
        write_kwargs['path2file'] = '/data/reprocessed.nc'
        self.database.write_langley_results(**write_kwargs)

        stored = self.database.dump_langley_results()
        self.assertEqual(stored.shape[0], 3)
        result = stored.loc[('2026-05-02', 'pm', 500)]
        self.assertEqual(result.intercept, 9.25)
        self.assertEqual(result.input_file, '/data/reprocessed.nc')

    def test_write_rejects_results_with_different_channels(self):
        metrics = self.langley.airmass_coverage_metrics.sel(channel=[500, 870])

        with self.assertRaisesRegex(ValueError, 'same channels'):
            self.database.write_langley_results(
                '2026-05-02',
                '/data/input.nc',
                'am',
                self.langley.langley_fitres,
                metrics,
            )

    def test_write_stores_channel_specific_coverage_values(self):
        metrics = self.langley.airmass_coverage_metrics.copy()
        metrics['retained_points'] = ('channel', [4, 5, 6])
        metrics['effective_airmass'] = ('channel', [3.1, 3.2, 3.3])

        self.database.write_langley_results(
            '2026-05-02',
            '/data/input.nc',
            'am',
            self.langley.langley_fitres,
            metrics,
        )

        stored = self.database.dump_langley_results()
        self.assertEqual(stored.loc[('2026-05-02', 'am', 500)].retained_points, 4)
        self.assertEqual(stored.loc[('2026-05-02', 'am', 940)].effective_airmass, 3.3)

    def test_add_empty_entry_stores_provenance_without_results(self):
        self.database.add_empty_entry(
            '2026-05-02',
            '/data/input.nc',
            'am',
            500,
            processing_date='2026-05-03T01:02:03',
            processing_server='test-host',
        )

        result = self.database.dump_langley_results().loc[
            ('2026-05-02', 'am', 500)
        ]
        self.assertEqual(result.input_file, '/data/input.nc')
        self.assertEqual(result.processed_at, '2026-05-03T01:02:03')
        self.assertEqual(result.process_version, '1.2.3')
        self.assertEqual(result.processing_server, 'test-host')
        self.assertTrue(result.drop([
            'input_file',
            'processed_at',
            'process_version',
            'processing_server',
        ]).isna().all())


if __name__ == '__main__':
    unittest.main()
