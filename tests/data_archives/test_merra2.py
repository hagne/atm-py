import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import xarray as xr

timezonefinder_stub = types.ModuleType("timezonefinder")


class TimezoneFinder:
    def timezone_at(self, lng, lat):
        if lng > -100:
            return "America/Chicago"
        return "America/Denver"


timezonefinder_stub.TimezoneFinder = TimezoneFinder
sys.modules.setdefault("timezonefinder", timezonefinder_stub)

from atmPy.data_archives.nasa import merra2


class EarthaccessStub:
    def __init__(self):
        self.login_kwargs = None
        self.login_environment = None
        self.search_kwargs = None
        self.dataset_search_kwargs = None
        self.download_args = None

    def login(self, **kwargs):
        self.login_kwargs = kwargs
        self.login_environment = {
            "EARTHDATA_TOKEN": os.environ.get("EARTHDATA_TOKEN"),
            "EARTHDATA_USERNAME": os.environ.get("EARTHDATA_USERNAME"),
            "EARTHDATA_PASSWORD": os.environ.get("EARTHDATA_PASSWORD"),
        }

    def search_datasets(self, **kwargs):
        self.dataset_search_kwargs = kwargs
        return [CollectionResultStub()]

    def search_data(self, **kwargs):
        self.search_kwargs = kwargs
        return ["granule-1", "granule-2"]

    def download(self, granules, download_dir):
        self.download_args = (granules, download_dir)
        return [download_dir / "day-1.nc4", download_dir / "day-2.nc4"]


class CollectionResultStub:
    def concept_id(self):
        return "C123456789-GES_DISC"


class DownloadFutureStub:
    def __init__(self, path):
        self.path = path

    def result(self):
        return self.path


class HarmonyProcessingFailedError(Exception):
    pass


class HarmonyClientStub:
    def __init__(self, auth_kwargs, processing_error=None):
        self.auth_kwargs = auth_kwargs
        self.processing_error = processing_error
        self.submitted_request = None
        self.wait_args = None
        self.download_kwargs = None

    def submit(self, request):
        self.submitted_request = request
        return "harmony-job-id"

    def wait_for_processing(self, job_id, show_progress=False):
        self.wait_args = (job_id, show_progress)
        if self.processing_error is not None:
            raise self.processing_error

    def download_all(self, job_id, **kwargs):
        self.download_kwargs = {"job_id": job_id, **kwargs}
        directory = Path(kwargs["directory"])
        return [DownloadFutureStub(directory / "subset.nc4")]


class HarmonyStub:
    client = types.SimpleNamespace(
        ProcessingFailedException=HarmonyProcessingFailedError
    )

    def __init__(self, processing_error=None):
        self.processing_error = processing_error
        self.client_instance = None

    def Client(self, **kwargs):
        self.client_instance = HarmonyClientStub(
            kwargs,
            processing_error=self.processing_error,
        )
        return self.client_instance

    @staticmethod
    def Collection(concept_id):
        return types.SimpleNamespace(id=concept_id)

    @staticmethod
    def BBox(west, south, east, north):
        return types.SimpleNamespace(
            west=west,
            south=south,
            east=east,
            north=north,
        )

    @staticmethod
    def Dimension(*, name, min, max):
        return types.SimpleNamespace(name=name, min=min, max=max)

    @staticmethod
    def Request(**kwargs):
        return types.SimpleNamespace(**kwargs)


class MissingCredentialsError(Exception):
    pass


class MissingCredentialsEarthaccessStub:
    exceptions = types.SimpleNamespace(
        LoginStrategyUnavailable=MissingCredentialsError
    )

    def login(self, **kwargs):
        raise MissingCredentialsError


class EulaNotAcceptedError(Exception):
    pass


class EulaNotAcceptedEarthaccessStub(EarthaccessStub):
    exceptions = types.SimpleNamespace(
        EulaNotAccepted=EulaNotAcceptedError,
    )

    def download(self, granules, download_dir):
        raise EulaNotAcceptedError


def make_dataset(units="Dobsons"):
    datetime = np.array(
        [
            "2024-01-01T00:30:00",
            "2024-01-01T01:30:00",
            "2024-01-01T02:30:00",
        ],
        dtype="datetime64[ns]",
    )
    ds = xr.Dataset(
        {
            "TO3": (
                ("time", "lat", "lon"),
                np.array([[[280.0]], [[281.0]], [[282.0]]]),
                {"units": units, "long_name": "total_column_ozone"},
            ),
            "T2M": (
                ("time", "lat", "lon"),
                np.array([[[270.0]], [[271.0]], [[272.0]]]),
            ),
        },
        coords={"time": datetime, "lat": [40.0], "lon": [-105.0]},
        attrs={"source": "MERRA2"},
    )
    return ds


class DownloadTotalColumnOzoneTests(unittest.TestCase):
    def test_downloads_expected_product_and_returns_requested_values(self):
        earthaccess = EarthaccessStub()
        source = make_dataset()

        with tempfile.TemporaryDirectory() as tmpdir:
            with (
                mock.patch.dict(os.environ, {}, clear=True),
                mock.patch.object(merra2, "earthaccess", earthaccess),
                mock.patch.object(
                    merra2.xr, "open_mfdataset", return_value=source
                ) as open_mfdataset,
            ):
                ds = merra2.download_total_column_ozone(
                    "2024-01-01T01:00:00Z",
                    "2024-01-01T02:00:00Z",
                    download_dir=tmpdir,
                    chunks={"time": 24},
                    settings_path=Path(tmpdir) / "earthdata.ini",
                    backend="earthaccess",
                )

        self.assertEqual(earthaccess.login_kwargs, {"strategy": "netrc"})
        self.assertEqual(earthaccess.search_kwargs["short_name"], "M2T1NXSLV")
        self.assertEqual(earthaccess.search_kwargs["version"], "5.12.4")
        self.assertEqual(
            earthaccess.search_kwargs["temporal"],
            ("2024-01-01T01:00:00Z", "2024-01-01T02:00:00Z"),
        )
        self.assertEqual(earthaccess.download_args[1], Path(tmpdir))
        open_mfdataset.assert_called_once_with(
            [Path(tmpdir) / "day-1.nc4", Path(tmpdir) / "day-2.nc4"],
            combine="by_coords",
            chunks={"time": 24},
        )
        self.assertEqual(list(ds.data_vars), ["TO3"])
        self.assertEqual(ds.TO3.attrs["units"], "Dobsons")
        self.assertEqual(ds.sizes["datetime"], 1)
        self.assertAlmostEqual(float(ds.TO3.isel(datetime=0, lat=0, lon=0)), 281.0)

    def test_rejects_unexpected_ozone_units(self):
        earthaccess = EarthaccessStub()

        with tempfile.TemporaryDirectory() as tmpdir:
            with (
                mock.patch.object(merra2, "earthaccess", earthaccess),
                mock.patch.object(
                    merra2.xr,
                    "open_mfdataset",
                    return_value=make_dataset(units="kg m-2"),
                ),
            ):
                with self.assertRaisesRegex(ValueError, "TO3 units must be 'Dobsons'"):
                    merra2.download_total_column_ozone(
                        "2024-01-01",
                        "2024-01-02",
                        download_dir=tmpdir,
                        settings_path=Path(tmpdir) / "earthdata.ini",
                        backend="earthaccess",
                    )

    def test_rejects_reversed_time_period(self):
        with self.assertRaisesRegex(ValueError, "end must be on or after start"):
            merra2.download_total_column_ozone("2024-01-02", "2024-01-01")

    def test_date_only_end_includes_the_entire_day(self):
        client = merra2.Merra2(backend="earthaccess")
        with (
            mock.patch.object(client, "login"),
            mock.patch.object(client, "_download_granules"),
            mock.patch.object(
                merra2.xr,
                "open_mfdataset",
                return_value=make_dataset(),
            ),
        ):
            client.paths = [Path("day.nc4")]
            ds = client.download(
                "2024-01-01",
                "2024-01-01",
                variables="TO3",
            )

        self.assertEqual(ds.sizes["datetime"], 3)

    def test_uses_environment_credentials_without_interactive_prompt(self):
        earthaccess = EarthaccessStub()
        environment = {
            "EARTHDATA_USERNAME": "test-user",
            "EARTHDATA_PASSWORD": "test-password",
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            with (
                mock.patch.dict(os.environ, environment, clear=True),
                mock.patch.object(merra2, "earthaccess", earthaccess),
                mock.patch.object(
                    merra2.xr, "open_mfdataset", return_value=make_dataset()
                ),
            ):
                merra2.download_total_column_ozone(
                    "2024-01-01",
                    "2024-01-02",
                    download_dir=tmpdir,
                    settings_path=Path(tmpdir) / "earthdata.ini",
                    backend="earthaccess",
                )

        self.assertEqual(earthaccess.login_kwargs, {"strategy": "environment"})

    def test_missing_credentials_raise_clear_error_without_prompt(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            settings_path = Path(tmpdir) / "earthdata.ini"
            with (
                mock.patch.dict(os.environ, {}, clear=True),
                mock.patch.object(
                    merra2, "earthaccess", MissingCredentialsEarthaccessStub()
                ),
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "A settings template was created"
                ):
                    merra2.Merra2(
                        settings_path=settings_path,
                        auth_strategy="auto",
                    ).login()

            self.assertEqual(
                settings_path.read_text(),
                "[earthdata]\n"
                "username = YOUR_USERNAME\n"
                "password = YOUR_PASSWORD\n",
            )

    def test_unaccepted_eula_raises_actionable_error(self):
        earthaccess = EulaNotAcceptedEarthaccessStub()
        with tempfile.TemporaryDirectory() as tmpdir:
            settings_path = Path(tmpdir) / "earthdata.ini"
            with (
                mock.patch.dict(
                    os.environ,
                    {
                        "EARTHDATA_USERNAME": "test-user",
                        "EARTHDATA_PASSWORD": "test-password",
                    },
                    clear=True,
                ),
                mock.patch.object(merra2, "earthaccess", earthaccess),
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    "users/test-user/unaccepted_eulas",
                ):
                    merra2.download_total_column_ozone(
                        "2024-01-01",
                        "2024-01-02",
                        download_dir=tmpdir,
                        settings_path=settings_path,
                        backend="earthaccess",
                    )

    def test_harmony_subsets_generic_collection_and_variables(self):
        earthaccess = EarthaccessStub()
        harmony = HarmonyStub()
        source = make_dataset()

        with tempfile.TemporaryDirectory() as tmpdir:
            settings_path = Path(tmpdir) / "earthdata.ini"
            settings_path.write_text(
                "[earthdata]\n"
                "username = test-user\n"
                "password = test-password\n"
            )
            client = merra2.Merra2(
                settings_path=settings_path,
                download_dir=tmpdir,
                short_name="M2I3NPASM",
                version="5.12.4",
            )
            with (
                mock.patch.dict(os.environ, {}, clear=True),
                mock.patch.object(merra2, "earthaccess", earthaccess),
                mock.patch.object(merra2, "harmony", harmony),
                mock.patch.object(
                    merra2.xr,
                    "open_mfdataset",
                    return_value=source,
                ),
            ):
                ds = client.download(
                    "2024-01-01T01:00:00Z",
                    "2024-01-01T02:00:00Z",
                    variables=["T2M", "TO3"],
                    bbox=(-110, 35, -100, 45),
                    dimensions={"lev": (850, 1000)},
                    show_progress=False,
                )

        request = harmony.client_instance.submitted_request
        self.assertEqual(
            earthaccess.dataset_search_kwargs,
            {"short_name": "M2I3NPASM", "version": "5.12.4"},
        )
        self.assertEqual(
            harmony.client_instance.auth_kwargs,
            {"auth": ("test-user", "test-password")},
        )
        self.assertEqual(request.collection.id, "C123456789-GES_DISC")
        self.assertEqual(request.variables, ["T2M", "TO3"])
        self.assertEqual(
            (
                request.spatial.west,
                request.spatial.south,
                request.spatial.east,
                request.spatial.north,
            ),
            (-110, 35, -100, 45),
        )
        self.assertEqual(
            (
                request.dimensions[0].name,
                request.dimensions[0].min,
                request.dimensions[0].max,
            ),
            ("lev", 850, 1000),
        )
        self.assertEqual(
            harmony.client_instance.wait_args,
            ("harmony-job-id", False),
        )
        self.assertEqual(list(ds.data_vars), ["T2M", "TO3"])
        self.assertEqual(ds.sizes["datetime"], 1)
        self.assertTrue(
            Path(
                harmony.client_instance.download_kwargs["directory"]
            ).is_relative_to(Path(tmpdir) / "harmony")
        )

    def test_harmony_403_explains_service_failure_and_full_file_fallback(self):
        earthaccess = EarthaccessStub()
        harmony = HarmonyStub(
            processing_error=HarmonyProcessingFailedError(
                "Subsetter failed with 403 error retrieving source granule"
            )
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            client = merra2.Merra2(
                settings_path=Path(tmpdir) / "earthdata.ini",
                download_dir=tmpdir,
                auth_strategy="netrc",
            )
            with (
                mock.patch.object(merra2, "earthaccess", earthaccess),
                mock.patch.object(merra2, "harmony", harmony),
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    "backend='earthaccess'",
                ):
                    client.download(
                        "2026-05-13",
                        "2026-05-13",
                        variables="TO3",
                        show_progress=False,
                    )


class EarthdataSettingsTests(unittest.TestCase):
    def test_creates_settings_file_with_secure_permissions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "earthdata.ini"

            returned_path = merra2.create_earthdata_settings_file(path)

            self.assertEqual(returned_path, path)
            self.assertEqual(
                path.read_text(),
                "[earthdata]\n"
                "username = YOUR_USERNAME\n"
                "password = YOUR_PASSWORD\n",
            )
            self.assertEqual(path.stat().st_mode & 0o777, 0o600)

    def test_auto_authentication_uses_completed_settings_file(self):
        earthaccess = EarthaccessStub()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "earthdata.ini"
            path.write_text(
                "[earthdata]\n"
                "username = test-user\n"
                "password = test-password\n"
            )

            with (
                mock.patch.dict(os.environ, {}, clear=True),
                mock.patch.object(merra2, "earthaccess", earthaccess),
            ):
                client = merra2.Merra2(
                    settings_path=path,
                    auth_strategy="auto",
                )
                client.login()

        self.assertEqual(earthaccess.login_kwargs, {"strategy": "environment"})
        self.assertEqual(client.username, "test-user")
        self.assertEqual(client.password, "test-password")
        self.assertEqual(
            earthaccess.login_environment,
            {
                "EARTHDATA_TOKEN": None,
                "EARTHDATA_USERNAME": "test-user",
                "EARTHDATA_PASSWORD": "test-password",
            },
        )

    def test_settings_authentication_does_not_leave_credentials_in_environment(self):
        earthaccess = EarthaccessStub()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "earthdata.ini"
            path.write_text(
                "[earthdata]\n"
                "username = test-user\n"
                "password = test-password\n"
            )
            existing_environment = {
                "EARTHDATA_TOKEN": "existing-token",
                "EARTHDATA_USERNAME": "existing-user",
                "EARTHDATA_PASSWORD": "existing-password",
            }

            with (
                mock.patch.dict(os.environ, existing_environment, clear=True),
                mock.patch.object(merra2, "earthaccess", earthaccess),
            ):
                client = merra2.Merra2(
                    settings_path=path,
                    auth_strategy="settings",
                )
                client.login()

                self.assertEqual(
                    {name: os.environ.get(name) for name in existing_environment},
                    existing_environment,
                )

        self.assertEqual(
            earthaccess.login_environment,
            {
                "EARTHDATA_TOKEN": None,
                "EARTHDATA_USERNAME": "test-user",
                "EARTHDATA_PASSWORD": "test-password",
            },
        )
        self.assertEqual(earthaccess.login_kwargs, {"strategy": "environment"})

    def test_class_stores_product_configuration_and_result(self):
        earthaccess = EarthaccessStub()
        source = make_dataset()
        with tempfile.TemporaryDirectory() as tmpdir:
            expected_paths = [
                Path(tmpdir) / "day-1.nc4",
                Path(tmpdir) / "day-2.nc4",
            ]
            client = merra2.Merra2(
                settings_path=Path(tmpdir) / "earthdata.ini",
                download_dir=tmpdir,
                chunks={"time": 24},
                auth_strategy="netrc",
                backend="earthaccess",
            )
            with (
                mock.patch.object(merra2, "earthaccess", earthaccess),
                mock.patch.object(
                    merra2.xr,
                    "open_mfdataset",
                    return_value=source,
                ),
            ):
                ds = client.download_total_column_ozone(
                    "2024-01-01T01:00:00",
                    "2024-01-01T02:00:00",
                )

        self.assertEqual(client.short_name, "M2T1NXSLV")
        self.assertEqual(client.ozone_variable, "TO3")
        self.assertIs(client.dataset, ds)
        self.assertEqual(client.paths, expected_paths)


if __name__ == "__main__":
    unittest.main()
