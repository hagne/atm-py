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


class GamClimatologyTest(unittest.TestCase):
    def test_prediction_confidence_is_symmetric_half_width(self):
        class FakeGam:
            coef_ = np.array([1.5])

            @staticmethod
            def predict(features):
                return 2 + features[:, 0]

            @classmethod
            def prediction_intervals(cls, features, width):
                prediction = cls.predict(features)
                confidence = np.linspace(0.25, 0.75, features.shape[0])
                return np.column_stack(
                    (prediction - confidence, prediction + confidence)
                )

        datetime = pd.date_range("2020-01-01", periods=4, freq="30D")
        model = statistics.GamClimatology()
        model.interactions_only = True
        model._start_date = datetime[0]
        model._data = pd.DataFrame(
            {
                "y": [2.0, 3.0, 4.0, 5.0],
                "x_dsincestart": [0.0, 30.0, 60.0, 90.0],
                "x_doy": [0.0, 30.0, 60.0, 90.0],
            },
            index=datetime,
        )
        model._fit_res = FakeGam()

        prediction = model.prediction

        np.testing.assert_allclose(
            prediction.prediction_confidence,
            np.linspace(0.25, 0.75, datetime.size),
        )
        self.assertEqual(
            prediction.prediction_confidence.dims,
            ("datetime_full",),
        )
        self.assertEqual(
            prediction.prediction_confidence.attrs["confidence_level"],
            model.prediction_confidence,
        )


class GamClimatologyRTest(unittest.TestCase):
    def _model(self, datetime, values, *, grid_size=41):
        data = pd.DataFrame({"value": values}, index=datetime)
        parent = ParentTimeSeries(data)
        model = statistics.Statistics(parent).gamcl_r
        model.trend_nsplines = 10
        model.seasonality_nsplines = 12
        model.prediction_grid_size = grid_size
        return model

    def _serial_model(self, datetime=None):
        if datetime is None:
            datetime = pd.date_range("2019-01-01", periods=160, freq="7D")
        elapsed = ((datetime - datetime[0]) / pd.Timedelta(days=1)).to_numpy()
        doy = datetime.day_of_year.to_numpy() - 1
        innovations = 0.04 * np.sin(np.arange(datetime.size) * 0.37)
        noise = np.empty(datetime.size)
        noise[0] = innovations[0]
        for index in range(1, datetime.size):
            noise[index] = 0.65 * noise[index - 1] + innovations[index]
        values = (
            2.0
            + 0.0007 * elapsed
            + 0.5 * np.sin(2 * np.pi * doy / 365.2425)
            + noise
        )
        model = self._model(datetime, values, grid_size=21)
        model.trend_nsplines = 8
        model.seasonality_nsplines = 10
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

    def test_serial_correlation_options_validate_and_invalidate_fit(self):
        datetime = pd.date_range("2020-01-01", periods=12, freq="30D")
        model = self._model(datetime, np.arange(datetime.size, dtype=float))

        self.assertFalse(model.ncv)
        self.assertFalse(model.ar1)
        self.assertFalse(model.newey_west)
        self.assertIsNone(model.ncv_window_days)
        self.assertIsNone(model.ar1_rho)
        self.assertIsNone(model.newey_west_lag)

        for attribute in ("ncv", "ar1", "newey_west"):
            with self.subTest(attribute=attribute):
                with self.assertRaises(TypeError):
                    setattr(model, attribute, 1)

        model._prediction = object()
        model._fit_res = object()
        model._fit_metadata = object()
        model._autocorrelation = object()
        model.ncv = True
        self.assertIsNone(model._prediction)
        self.assertIsNone(model._fit_res)
        self.assertIsNone(model._fit_metadata)
        self.assertIsNone(model._autocorrelation)
        with self.assertRaises(ValueError):
            model.ar1 = True
        with self.assertRaises(ValueError):
            model.newey_west = True

        model.ncv = False
        model.ar1 = True
        with self.assertRaises(ValueError):
            model.ncv = True
        with self.assertRaises(ValueError):
            model.newey_west = True

        model.ar1 = False
        model.newey_west = True
        with self.assertRaises(ValueError):
            model.ncv = True
        with self.assertRaises(ValueError):
            model.ar1 = True

        for attribute, valid_value, invalid_values in (
            ("ncv_window_days", 14.5, (0, -1, np.inf, "14")),
            ("ar1_rho", 0.4, (-1, 1, np.nan, "0.4")),
            ("newey_west_lag", 0, (-1, 2.5, "3")),
        ):
            with self.subTest(attribute=attribute):
                setattr(model, attribute, valid_value)
                self.assertEqual(getattr(model, attribute), valid_value)
                model._prediction = object()
                setattr(model, attribute, None)
                self.assertIsNone(getattr(model, attribute))
                self.assertIsNone(model._prediction)
                for invalid_value in invalid_values:
                    with self.assertRaises((TypeError, ValueError)):
                        setattr(model, attribute, invalid_value)

    def test_missing_rscript_has_actionable_error(self):
        datetime = pd.date_range("2020-01-01", periods=12, freq="30D")
        model = self._model(datetime, np.arange(datetime.size, dtype=float))
        model.rscript_path = "/path/that/does/not/contain/Rscript"

        with self.assertRaisesRegex(RuntimeError, "r-base.*r-mgcv"):
            _ = model.prediction

    def test_serial_correlation_options_are_forwarded_to_rscript(self):
        datetime = pd.date_range("2020-01-01", periods=12, freq="30D")
        values = np.arange(datetime.size, dtype=float)
        cases = (
            (None, None, "none", "NA"),
            ("ncv", None, "ncv", "NA"),
            ("ar1", 0.6, "ar1", "0.6"),
            ("newey_west", 0, "newey_west", "0.0"),
        )

        for method, parameter, expected_method, expected_parameter in cases:
            with self.subTest(method=method, parameter=parameter):
                model = self._model(datetime, values)
                if method is not None:
                    setattr(model, method, True)
                    if parameter is not None:
                        parameter_name = {
                            "ncv": "ncv_window_days",
                            "ar1": "ar1_rho",
                            "newey_west": "newey_west_lag",
                        }[method]
                        setattr(model, parameter_name, parameter)
                failed_process = subprocess.CompletedProcess(
                    [], 1, stdout="", stderr="expected test stop"
                )
                with (
                    mock.patch.object(
                        model,
                        "_resolve_rscript",
                        return_value="/fake/bin/Rscript",
                    ),
                    mock.patch(
                        "atmPy.general.statistics.subprocess.run",
                        return_value=failed_process,
                    ) as run,
                    self.assertRaisesRegex(RuntimeError, "expected test stop"),
                ):
                    _ = model.prediction

                command = run.call_args.args[0]
                self.assertEqual(
                    command[-2:],
                    [expected_method, expected_parameter],
                )

    def test_rscript_protocol_and_output_mapping_without_r(self):
        datetime = pd.date_range("2020-01-01", periods=12, freq="30D")
        values = np.arange(datetime.size, dtype=float)
        model = self._model(datetime, values, grid_size=3)
        model.robust = True
        model.newey_west_lag = 3
        model.newey_west = True
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
                    "autocorrelation_residual": values - fitted,
                    "confidence": np.full(values.size, 0.5),
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
                        "serial_correlation_method",
                        "serial_correlation_modelled",
                        "uncertainty_adjusted_for_serial_correlation",
                        "covariance_description",
                        "ncv_window_days",
                        "ncv_window_estimated",
                        "ar1_rho",
                        "ar1_rho_estimated",
                        "serial_correlation_cadence_days",
                        "serial_correlation_sections",
                        "newey_west_lag",
                        "newey_west_lag_estimated",
                        "newey_west_kernel",
                        "newey_west_score_centered",
                        "newey_west_small_sample_correction",
                        "ar1_approximation",
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
                        "newey_west",
                        "FALSE",
                        "TRUE",
                        "a coefficient-level Newey-West HAC covariance",
                        "",
                        "FALSE",
                        "",
                        "FALSE",
                        "30",
                        "1",
                        "3",
                        "FALSE",
                        "Bartlett",
                        "TRUE",
                        "none",
                        "none",
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
        self.assertEqual(command[15:], ["newey_west", "3.0"])
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
        np.testing.assert_allclose(prediction.prediction_confidence, 0.5)
        self.assertEqual(
            prediction.prediction_confidence.dims,
            ("datetime_full",),
        )
        self.assertEqual(
            prediction.prediction_confidence.attrs["description"],
            "Pointwise symmetric prediction half-width for an individual "
            "response; includes residual observation variability. Add and "
            "subtract this value from prediction to construct the interval.",
        )
        self.assertEqual(
            prediction.prediction_confidence.attrs["confidence_level"],
            model.prediction_confidence,
        )
        self.assertEqual(
            prediction.prediction_confidence.attrs["interval_construction"],
            "lower = prediction - prediction_confidence; "
            "upper = prediction + prediction_confidence",
        )
        self.assertEqual(
            prediction.prediction_mean_confidence.dims,
            ("datetime_full",),
        )
        self.assertEqual(
            prediction.prediction_mean_confidence.attrs["description"],
            "Pointwise symmetric confidence half-width for the fitted mean "
            "based on a coefficient-level Newey-West HAC covariance; excludes "
            "residual observation variability. Add and subtract this value from "
            "prediction to construct the interval.",
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
        np.testing.assert_allclose(prediction.residual_for_autocorrelation, -0.25)
        np.testing.assert_allclose(prediction.partial_datetime, [0.0, 0.5, 1.0])
        np.testing.assert_allclose(prediction.partial_doy, [-0.5, 0.0, 0.5])
        self.assertEqual(prediction.attrs["family"], "Scaled t(4,1)")
        self.assertEqual(prediction.attrs["robust"], 1)
        self.assertEqual(prediction.attrs["outer_convergence"], "full convergence")
        self.assertEqual(prediction.attrs["r_version"], "4.5.0")
        self.assertEqual(prediction.attrs["serial_correlation_method"], "newey_west")
        self.assertEqual(prediction.attrs["serial_correlation_modelled"], 0)
        self.assertEqual(
            prediction.attrs["uncertainty_adjusted_for_serial_correlation"],
            1,
        )
        self.assertEqual(prediction.attrs["newey_west_lag"], 3)
        self.assertEqual(prediction.attrs["newey_west_lag_estimated"], 0)
        self.assertEqual(prediction.attrs["newey_west_kernel"], "Bartlett")
        self.assertEqual(prediction.attrs["newey_west_score_centered"], 1)
        self.assertEqual(
            prediction.attrs["newey_west_small_sample_correction"],
            "none",
        )
        self.assertEqual(
            prediction.attrs["prediction_interval"],
            "approximate pointwise response interval",
        )

        axis = mock.Mock()
        axis.get_figure.return_value = object()
        line = mock.Mock()
        line.get_color.return_value = "C0"
        axis.get_lines.return_value = [line]
        with mock.patch(
            "xarray.plot.accessor.DataArrayPlotAccessor.__call__"
        ):
            model.plot_prediction(
                ax=axis,
                show_confidence=True,
                show_original_data=False,
            )
        fill_args = axis.fill_between.call_args.args
        np.testing.assert_allclose(fill_args[1], prediction.prediction + 0.5)
        np.testing.assert_allclose(fill_args[2], prediction.prediction - 0.5)

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
            ("datetime_full",),
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
        response_critical_value = scipy.stats.t.ppf(
            (1 + model.prediction_confidence) / 2,
            prediction.attrs["residual_degrees_of_freedom"],
        )
        expected_prediction_confidence = response_critical_value * np.sqrt(
            prediction.prediction_standard_error**2
            + prediction.attrs["residual_scale"] ** 2
        )
        np.testing.assert_allclose(
            prediction.prediction_confidence,
            expected_prediction_confidence,
        )
        self.assertTrue(
            bool(
                (
                    prediction.prediction_mean_confidence
                    < prediction.prediction_confidence
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
    def test_robust_ncv_uses_explicit_time_window(self):
        model = self._serial_model()
        model.robust = True
        model.ncv_window_days = 21
        model.ncv = True

        prediction = model.prediction

        self.assertTrue(prediction.attrs["family"].lower().startswith("scaled t"))
        self.assertEqual(prediction.attrs["serial_correlation_method"], "ncv")
        self.assertEqual(prediction.attrs["serial_correlation_modelled"], 0)
        self.assertEqual(
            prediction.attrs["uncertainty_adjusted_for_serial_correlation"],
            1,
        )
        self.assertEqual(prediction.attrs["ncv_window_days"], 21)
        self.assertEqual(prediction.attrs["ncv_window_estimated"], 0)
        self.assertTrue(bool(np.isfinite(prediction.prediction_standard_error).all()))

        critical_value = scipy.stats.norm.ppf(
            (1 + model.prediction_confidence) / 2
        )
        np.testing.assert_allclose(
            prediction.prediction_mean_confidence,
            critical_value * prediction.prediction_standard_error,
        )

    @unittest.skipUnless(R_MGCV_AVAILABLE, "Rscript with mgcv is not available")
    def test_ar1_fit_reports_structure_and_finite_intervals(self):
        model = self._serial_model()
        model.ar1_rho = 0.5
        model.ar1 = True

        prediction = model.prediction

        self.assertEqual(prediction.attrs["serial_correlation_method"], "ar1")
        self.assertEqual(prediction.attrs["serial_correlation_modelled"], 1)
        self.assertEqual(prediction.attrs["ar1_rho"], 0.5)
        self.assertEqual(prediction.attrs["ar1_rho_estimated"], 0)
        self.assertEqual(prediction.attrs["serial_correlation_cadence_days"], 7)
        self.assertEqual(prediction.attrs["serial_correlation_sections"], 1)
        autocorrelation = model.autocorrelation
        self.assertEqual(
            autocorrelation.attrs["residual_source"],
            "residual_for_autocorrelation",
        )
        self.assertTrue(bool(np.isfinite(prediction.prediction_confidence).all()))
        self.assertTrue(
            bool(
                (
                    prediction.prediction_confidence
                    > prediction.prediction_mean_confidence
                ).all()
            )
        )

    @unittest.skipUnless(R_MGCV_AVAILABLE, "Rscript with mgcv is not available")
    def test_newey_west_uses_point_dependent_covariance(self):
        baseline = self._serial_model().prediction
        model = self._serial_model()
        model.newey_west_lag = 4
        model.newey_west = True

        prediction = model.prediction

        np.testing.assert_allclose(prediction.prediction, baseline.prediction)
        adjustment = (
            prediction.prediction_standard_error
            / baseline.prediction_standard_error
        )
        self.assertGreater(float(adjustment.max() - adjustment.min()), 0.01)
        self.assertEqual(
            prediction.attrs["serial_correlation_method"],
            "newey_west",
        )
        self.assertEqual(prediction.attrs["newey_west_lag"], 4)
        self.assertEqual(prediction.attrs["newey_west_lag_estimated"], 0)
        self.assertEqual(prediction.attrs["newey_west_kernel"], "Bartlett")
        np.testing.assert_allclose(
            model.confidence_interval_adjusted_newey_wet,
            prediction.partial_datetime_confidence,
        )
        np.testing.assert_allclose(
            model.confidence_interval_adjusted_newey_west,
            prediction.partial_datetime_confidence,
        )

    @unittest.skipUnless(R_MGCV_AVAILABLE, "Rscript with mgcv is not available")
    def test_ar1_rejects_irregular_timestamps(self):
        datetime = pd.date_range(
            "2019-01-01", periods=160, freq="7D"
        ).to_numpy(copy=True)
        datetime[50] += np.timedelta64(2, "D")
        model = self._serial_model(pd.DatetimeIndex(datetime))
        model.ar1_rho = 0.5
        model.ar1 = True

        with self.assertRaisesRegex(RuntimeError, "requires regularly spaced data"):
            _ = model.prediction

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
        self.assertEqual(
            prediction.interactions_standard_error.dims,
            ("datetime", "doy"),
        )
        self.assertEqual(
            prediction.interactions_confidence.dims,
            ("datetime", "doy", "boundary"),
        )
        self.assertTrue(
            bool(
                (
                    prediction.interactions_confidence.sel(boundary="bottom")
                    <= prediction.interactions
                ).all()
            )
        )
        self.assertTrue(
            bool(
                (
                    prediction.interactions
                    <= prediction.interactions_confidence.sel(boundary="top")
                ).all()
            )
        )
        self.assertIn("te(days_since_start, doy", prediction.attrs["formula"])


if __name__ == "__main__":
    unittest.main()
