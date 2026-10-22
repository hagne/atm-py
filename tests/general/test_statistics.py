import pathlib
import shutil
import subprocess
import unittest
from unittest import mock

import numpy as np
import pandas as pd
import scipy
import scipy.stats

from atmPy.general import statistics


def _r_mgcv_available():
    rscript = shutil.which("Rscript")
    if rscript is None:
        return False

    try:
        result = subprocess.run(
            [
                rscript,
                "--vanilla",
                "-e",
                "quit(status = ifelse(requireNamespace('mgcv', quietly = TRUE), 0, 1))",
            ],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


R_MGCV_AVAILABLE = _r_mgcv_available()


class ParentTimeSeries:
    def __init__(self, data):
        self.data = data


class GamClimatologyRTest(unittest.TestCase):
    def _model(self, datetime, values, *, grid_size=41):
        data = pd.DataFrame({"value": values}, index=datetime)
        parent = ParentTimeSeries(data)
        model = statistics.Statistics(parent).gamcl_r
        model.trend_nsplines = 10
        model.seasonality_nsplines = 12
        model.prediction_grid_size = grid_size
        return model

    def test_statistics_exposes_r_climatology_with_robust_parameter(self):
        datetime = pd.date_range("2020-01-01", periods=12, freq="30D")
        model = self._model(datetime, np.arange(datetime.size, dtype=float))

        self.assertIsInstance(model, statistics.GamClimatologyR)
        self.assertFalse(model.robust)
        model.robust = True
        self.assertTrue(model.robust)
        with self.assertRaises(NotImplementedError):
            model.interactions_lam = 1

    def test_missing_rscript_has_actionable_error(self):
        datetime = pd.date_range("2020-01-01", periods=12, freq="30D")
        model = self._model(datetime, np.arange(datetime.size, dtype=float))
        model.rscript_path = "/path/that/does/not/contain/Rscript"

        with self.assertRaisesRegex(RuntimeError, "r-base.*r-mgcv"):
            _ = model.prediction

    def test_rscript_protocol_and_output_mapping_without_r(self):
        datetime = pd.date_range("2020-01-01", periods=12, freq="30D")
        values = np.arange(datetime.size, dtype=float)
        model = self._model(datetime, values, grid_size=3)
        model.robust = True
        captured = {}

        def fake_rscript(command, **kwargs):
            captured["command"] = command
            captured["kwargs"] = kwargs
            captured["input"] = pd.read_csv(command[3])
            output_dir = pathlib.Path(command[4])

            fitted = values + 0.25
            pd.DataFrame(
                {
                    "row_id": np.arange(values.size),
                    "prediction": fitted,
                    "standard_error": np.full(values.size, 0.1),
                    "lower": fitted - 0.5,
                    "upper": fitted + 0.5,
                    "mean_confidence": np.full(values.size, 0.2),
                }
            ).to_csv(output_dir / "prediction.csv", index=False)
            pd.DataFrame(
                {
                    "days_since_start": [0.0, 165.0, 330.0],
                    "effect": [0.0, 0.5, 1.0],
                    "standard_error": [0.1, 0.1, 0.1],
                    "lower": [-0.2, 0.3, 0.8],
                    "upper": [0.2, 0.7, 1.2],
                }
            ).to_csv(output_dir / "trend.csv", index=False)
            pd.DataFrame(
                {
                    "doy": [0.0, 182.5, 365.0],
                    "effect": [-0.5, 0.0, 0.5],
                    "standard_error": [0.1, 0.1, 0.1],
                    "lower": [-0.7, -0.2, 0.3],
                    "upper": [-0.3, 0.2, 0.7],
                }
            ).to_csv(output_dir / "seasonality.csv", index=False)
            pd.DataFrame(
                {
                    "key": [
                        "backend",
                        "family",
                        "robust",
                        "method",
                        "formula",
                        "intercept",
                        "nobs",
                        "effective_degrees_of_freedom",
                        "deviance_explained",
                        "residual_scale",
                        "robust_degrees_of_freedom",
                        "residual_degrees_of_freedom",
                        "smoothing_parameters",
                        "converged",
                        "outer_convergence",
                        "r_version",
                        "mgcv_version",
                    ],
                    "value": [
                        "mgcv",
                        "Scaled t(4,1)",
                        "TRUE",
                        "REML",
                        "y ~ s(days_since_start) + s(doy)",
                        "1.5",
                        "12",
                        "5.0",
                        "0.8",
                        "1.0",
                        "4.0",
                        "7.0",
                        "0.2,0.3",
                        "TRUE",
                        "full convergence",
                        "4.5.0",
                        "1.9-3",
                    ],
                }
            ).to_csv(output_dir / "metadata.tsv", sep="\t", index=False)
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

        with (
            mock.patch.object(
                model,
                "_resolve_rscript",
                return_value="/fake/bin/Rscript",
            ),
            mock.patch(
                "atmPy.general.statistics.subprocess.run",
                side_effect=fake_rscript,
            ) as run,
        ):
            prediction = model.prediction

        run.assert_called_once()
        command = captured["command"]
        self.assertEqual(command[:3], ["/fake/bin/Rscript", "--vanilla", "-"])
        self.assertEqual(command[5:10], ["TRUE", "spline", "10", "12", "3"])
        self.assertEqual(command[11:13], ["FALSE", "FALSE"])
        self.assertIn("requireNamespace(\"mgcv\", quietly = TRUE)", captured["kwargs"]["input"])
        self.assertTrue(captured["kwargs"]["capture_output"])
        self.assertTrue(captured["kwargs"]["text"])
        self.assertFalse(captured["kwargs"]["check"])
        self.assertEqual(
            captured["input"].columns.tolist(),
            ["row_id", "y", "days_since_start", "doy"],
        )
        np.testing.assert_array_equal(captured["input"].row_id, np.arange(12))
        np.testing.assert_allclose(prediction.prediction, values + 0.25)
        np.testing.assert_allclose(
            prediction.prediction_mean_confidence,
            0.2,
        )
        self.assertEqual(
            prediction.prediction_mean_confidence.dims,
            ("datetime_full",),
        )
        self.assertEqual(
            prediction.prediction_mean_confidence.attrs["description"],
            "Pointwise symmetric confidence half-width for the fitted mean "
            "based on mgcv's unconditional covariance; excludes residual "
            "observation variability. Add and subtract this value from prediction "
            "to construct the interval.",
        )
        self.assertEqual(
            prediction.prediction_mean_confidence.attrs["confidence_level"],
            model.prediction_confidence,
        )
        self.assertEqual(
            prediction.prediction_mean_confidence.attrs["interval_construction"],
            "lower = prediction - prediction_mean_confidence; "
            "upper = prediction + prediction_mean_confidence",
        )
        np.testing.assert_allclose(prediction.residual, -0.25)
        np.testing.assert_allclose(prediction.partial_datetime, [0.0, 0.5, 1.0])
        np.testing.assert_allclose(prediction.partial_doy, [-0.5, 0.0, 0.5])
        self.assertEqual(prediction.attrs["family"], "Scaled t(4,1)")
        self.assertEqual(prediction.attrs["robust"], 1)
        self.assertEqual(prediction.attrs["outer_convergence"], "full convergence")
        self.assertEqual(prediction.attrs["r_version"], "4.5.0")
        self.assertEqual(
            prediction.attrs["prediction_interval"],
            "approximate pointwise response interval",
        )

    @unittest.skipUnless(R_MGCV_AVAILABLE, "Rscript with mgcv is not available")
    def test_prediction_schema_and_missing_observation_alignment(self):
        datetime = pd.date_range("2018-01-01", periods=4 * 37, freq="10D")
        elapsed = (datetime - datetime[0]) / pd.Timedelta(days=1)
        doy = datetime.day_of_year.to_numpy() - 1
        values = (
            3.0
            + 0.0005 * elapsed.to_numpy()
            + 0.7 * np.sin(2 * np.pi * doy / 365.2425)
        )
        values[17] = np.nan
        model = self._model(datetime, values)

        prediction = model.prediction

        required_variables = {
            "partial_datetime",
            "partial_datetime_confidence",
            "partial_doy",
            "partial_doy_confidence",
            "prediction",
            "prediction_confidence",
            "prediction_mean_confidence",
            "residual",
        }
        self.assertTrue(required_variables.issubset(prediction.data_vars))
        self.assertEqual(prediction.partial_datetime.dims, ("datetime",))
        self.assertEqual(
            prediction.partial_datetime_confidence.dims,
            ("datetime", "boundary"),
        )
        self.assertEqual(prediction.partial_doy.dims, ("doy",))
        self.assertEqual(
            prediction.partial_doy_confidence.dims,
            ("doy", "boundary"),
        )
        self.assertEqual(prediction.prediction.dims, ("datetime_full",))
        self.assertEqual(
            prediction.prediction_confidence.dims,
            ("datetime_full", "boundary"),
        )
        self.assertEqual(
            prediction.prediction_mean_confidence.dims,
            ("datetime_full",),
        )
        self.assertEqual(
            prediction.prediction_mean_confidence.attrs["description"],
            "Pointwise symmetric confidence half-width for the fitted mean "
            "based on mgcv's unconditional covariance; excludes residual "
            "observation variability. Add and subtract this value from prediction "
            "to construct the interval.",
        )
        self.assertEqual(
            prediction.prediction_mean_confidence.attrs["confidence_level"],
            model.prediction_confidence,
        )
        self.assertEqual(
            prediction.prediction_mean_confidence.attrs["interval_construction"],
            "lower = prediction - prediction_mean_confidence; "
            "upper = prediction + prediction_mean_confidence",
        )
        self.assertEqual(prediction.residual.dims, ("datetime_full",))
        self.assertEqual(
            prediction.boundary.to_numpy().tolist(),
            ["bottom", "top"],
        )
        self.assertEqual(prediction.sizes["datetime"], model.prediction_grid_size)
        self.assertEqual(prediction.sizes["doy"], model.prediction_grid_size)
        np.testing.assert_array_equal(
            prediction.datetime_full.to_numpy(),
            datetime.to_numpy(),
        )
        self.assertTrue(np.isfinite(prediction.prediction).all())
        self.assertTrue(np.isfinite(prediction.attrs["intercept"]))
        self.assertEqual(prediction.attrs["backend"], "mgcv")
        self.assertEqual(prediction.attrs["family"], "gaussian")
        self.assertEqual(prediction.attrs["robust"], 0)
        self.assertEqual(prediction.attrs["converged"], 1)
        self.assertEqual(prediction.attrs["outer_convergence"], "full convergence")
        self.assertTrue(prediction.attrs["r_version"])
        self.assertTrue(prediction.attrs["mgcv_version"])

        expected_residual = values - prediction.prediction.to_numpy()
        np.testing.assert_allclose(
            prediction.residual,
            expected_residual,
            equal_nan=True,
        )
        self.assertTrue(np.isnan(prediction.residual.isel(datetime_full=17)))

        autocorrelation = model.autocorrelation
        self.assertEqual(float(autocorrelation.sel(lag=0, acf_params="acf")), 1.0)
        self.assertGreaterEqual(autocorrelation.attrs["lagmax_above_ci"], 0)

        for estimate_name, confidence_name in (
            ("partial_datetime", "partial_datetime_confidence"),
            ("partial_doy", "partial_doy_confidence"),
            ("prediction", "prediction_confidence"),
        ):
            estimate = prediction[estimate_name]
            confidence = prediction[confidence_name]
            self.assertTrue(
                bool((confidence.sel(boundary="bottom") <= estimate).all())
            )
            self.assertTrue(
                bool((estimate <= confidence.sel(boundary="top")).all())
            )

        critical_value = scipy.stats.norm.ppf(
            (1 + model.prediction_confidence) / 2
        )
        np.testing.assert_allclose(
            prediction.prediction_mean_confidence,
            critical_value * prediction.prediction_standard_error,
        )
        response_interval_width = (
            prediction.prediction_confidence.sel(boundary="top")
            - prediction.prediction_confidence.sel(boundary="bottom")
        )
        self.assertTrue(
            bool(
                (
                    2 * prediction.prediction_mean_confidence
                    < response_interval_width
                ).all()
            )
        )

    @unittest.skipUnless(R_MGCV_AVAILABLE, "Rscript with mgcv is not available")
    def test_robust_fit_resists_large_outliers(self):
        datetime = pd.date_range("2017-01-01", periods=6 * 52, freq="7D")
        elapsed = ((datetime - datetime[0]) / pd.Timedelta(days=1)).to_numpy()
        doy = datetime.day_of_year.to_numpy() - 1
        truth = (
            2.0
            + 0.0007 * elapsed
            + 0.8 * np.sin(2 * np.pi * doy / 365.2425)
        )
        contaminated = truth + 0.05 * np.sin(np.arange(datetime.size) * 1.7)
        contaminated[9::21] += 8.0
        model = self._model(datetime, contaminated)

        gaussian_prediction = model.prediction.prediction.to_numpy().copy()
        gaussian_result = model.prediction
        model.robust = True
        robust_result = model.prediction
        robust_prediction = robust_result.prediction.to_numpy()

        gaussian_rmse = np.sqrt(np.mean((gaussian_prediction - truth) ** 2))
        robust_rmse = np.sqrt(np.mean((robust_prediction - truth) ** 2))
        self.assertIsNot(gaussian_result, robust_result)
        self.assertTrue(robust_result.attrs["family"].lower().startswith("scaled t"))
        self.assertEqual(robust_result.attrs["robust"], 1)
        self.assertGreater(
            robust_result.attrs["robust_degrees_of_freedom"],
            2,
        )
        self.assertLess(robust_rmse, 0.1)
        self.assertLess(robust_rmse, gaussian_rmse * 0.25)

    @unittest.skipUnless(R_MGCV_AVAILABLE, "Rscript with mgcv is not available")
    def test_interaction_only_uses_tensor_basis(self):
        datetime = pd.date_range("2018-01-01", periods=96, freq="14D")
        elapsed = ((datetime - datetime[0]) / pd.Timedelta(days=1)).to_numpy()
        doy = datetime.day_of_year.to_numpy() - 1
        values = 2 + 0.2 * np.sin(2 * np.pi * doy / 365.2425) * (
            1 + elapsed / elapsed.max()
        )
        model = self._model(datetime, values, grid_size=7)
        model.interactions_only = True
        model.interactions_nsplines = (5, 6)
        model.trend_nsplines = 1
        model.seasonality_nsplines = 2

        prediction = model.prediction

        self.assertNotIn("partial_datetime", prediction)
        self.assertNotIn("partial_doy", prediction)
        self.assertEqual(prediction.interactions.dims, ("datetime", "doy"))
        self.assertEqual(prediction.interactions.shape, (7, 7))
        self.assertIn("te(days_since_start, doy", prediction.attrs["formula"])


if __name__ == "__main__":
    unittest.main()
