import sys
import tempfile
import types
import unittest
from datetime import datetime, timezone
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

from atmPy.data_archives.nasa import geos_fp


def make_dataset():
    time = np.array(
        [
            "2026-08-12T00:30:00",
            "2026-08-12T01:30:00",
            "2026-08-12T02:30:00",
        ],
        dtype="datetime64[ns]",
    )
    lat = [39.0, 40.0, 41.0]
    lon = [-106.0, -105.0, -104.0]
    lev = [1000.0, 850.0, 700.0]
    return xr.Dataset(
        data_vars={
            "t2m": (
                ("time", "lat", "lon"),
                270.0 + np.arange(27).reshape(3, 3, 3),
                {"units": "K"},
            ),
            "rh": (
                ("time", "lev", "lat", "lon"),
                np.arange(81, dtype=float).reshape(3, 3, 3, 3),
                {"units": "1"},
            ),
        },
        coords={"time": time, "lev": lev, "lat": lat, "lon": lon},
        attrs={"source": "synthetic GEOS-FP"},
    )


class GeosFpUrlTests(unittest.TestCase):
    def test_builds_assimilation_latest_and_cycle_urls(self):
        base = "https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg"

        cases = (
            (
                {"stream": "assim"},
                f"{base}/assim/tavg1_2d_slv_Nx",
            ),
            (
                {"stream": "fcast"},
                f"{base}/fcast/tavg1_2d_slv_Nx.latest",
            ),
            (
                {
                    "stream": "fcast",
                    "cycle": datetime(2026, 8, 12, 6),
                },
                f"{base}/fcast/tavg1_2d_slv_Nx/"
                "tavg1_2d_slv_Nx.20260812_06",
            ),
            (
                {"stream": "seamless"},
                f"{base}/seamless/tavg1_2d_slv_Nx.latest",
            ),
            (
                {
                    "stream": "seamless",
                    "cycle": "2026-08-12T06:00:00Z",
                },
                f"{base}/seamless/tavg1_2d_slv_Nx/"
                "tavg1_2d_slv_Nx.20260812_06",
            ),
        )

        for kwargs, expected in cases:
            with self.subTest(**kwargs):
                client = geos_fp.GeosFp(
                    collection="tavg1_2d_slv_Nx",
                    **kwargs,
                )
                self.assertEqual(client.opendap_url, expected)

    def test_cycle_is_converted_to_utc(self):
        client = geos_fp.GeosFp(
            stream="fcast",
            cycle=datetime(
                2026,
                8,
                12,
                0,
                tzinfo=timezone.utc,
            ),
        )

        self.assertTrue(client.opendap_url.endswith(".20260812_00"))

    def test_opens_the_public_dataset_with_the_netcdf4_backend(self):
        source = make_dataset()
        client = geos_fp.GeosFp(chunks={"time": 24})

        with (
            mock.patch.object(
                geos_fp.netcdf4,
                "module_available",
                True,
            ),
            mock.patch.object(
                geos_fp.xr,
                "open_dataset",
                return_value=source,
            ) as open_dataset,
        ):
            returned = client._open_remote_dataset()

        self.assertIs(returned, source)
        open_dataset.assert_called_once_with(
            client.opendap_url,
            engine="netcdf4",
            chunks={"time": 24},
        )

    def test_open_error_preserves_network_detail_and_suggests_ipv4(self):
        client = geos_fp.GeosFp()

        with (
            mock.patch.object(
                geos_fp.netcdf4,
                "module_available",
                True,
            ),
            mock.patch.object(
                geos_fp.xr,
                "open_dataset",
                side_effect=OSError("curl timeout"),
            ),
            self.assertRaisesRegex(
                OSError,
                "curl timeout.*GeosFp\\(force_ipv4=True\\)",
            ),
        ):
            client._open_remote_dataset()

    def test_ipv4_connection_resolves_only_ipv4_addresses(self):
        ssl_context = mock.Mock()
        upstream_socket = mock.Mock()
        address = ("169.154.151.145", 443)

        with (
            mock.patch.object(
                geos_fp.socket,
                "getaddrinfo",
                return_value=[
                    (
                        geos_fp.socket.AF_INET,
                        geos_fp.socket.SOCK_STREAM,
                        6,
                        "",
                        address,
                    )
                ],
            ) as getaddrinfo,
            mock.patch.object(
                geos_fp.socket,
                "socket",
                return_value=upstream_socket,
            ),
        ):
            connection = geos_fp._IPv4HTTPSConnection(
                "opendap.nccs.nasa.gov",
                443,
                timeout=30,
                context=ssl_context,
            )
            connection.connect()

        getaddrinfo.assert_called_once_with(
            "opendap.nccs.nasa.gov",
            443,
            geos_fp.socket.AF_INET,
            geos_fp.socket.SOCK_STREAM,
        )
        upstream_socket.connect.assert_called_once_with(address)
        ssl_context.wrap_socket.assert_called_once_with(
            upstream_socket,
            server_hostname="opendap.nccs.nasa.gov",
        )


class GeosFpDownloadTests(unittest.TestCase):
    def test_force_ipv4_routes_opendap_through_loopback(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            client = geos_fp.GeosFp(
                download_dir=tmpdir,
                chunks=None,
                force_ipv4=True,
            )
            local_url = (
                "http://127.0.0.1:54321/dods/GEOS-5/fp/0.25_deg/assim/"
                "tavg1_2d_slv_Nx"
            )
            with (
                mock.patch.object(
                    geos_fp,
                    "_nccs_ipv4_proxy",
                ) as ipv4_proxy,
                mock.patch.object(
                    client,
                    "_open_remote_dataset",
                    return_value=make_dataset(),
                ) as open_remote,
            ):
                ipv4_proxy.return_value.__enter__.return_value = local_url
                path = client.download_file(
                    "2026-08-12",
                    "2026-08-12",
                    variables="t2m",
                )

            ipv4_proxy.assert_called_once_with(client.opendap_url)
            open_remote.assert_called_once_with(local_url)
            self.assertTrue(path.is_file())

    def test_download_returns_requested_numeric_spatial_and_vertical_subset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            client = geos_fp.GeosFp(
                collection="tavg1_2d_slv_Nx",
                stream="assim",
                download_dir=tmpdir,
                chunks=None,
            )
            with mock.patch.object(
                client,
                "_open_remote_dataset",
                return_value=make_dataset(),
            ) as open_remote:
                ds = client.download(
                    "2026-08-12T01:00:00Z",
                    "2026-08-12T02:00:00Z",
                    variables=["t2m", "rh"],
                    bbox=(-105.25, 39.5, -103.75, 40.5),
                    dimensions={"lev": (800.0, 900.0)},
                )

            try:
                open_remote.assert_called_once_with()
                self.assertEqual(list(ds.data_vars), ["t2m", "rh"])
                self.assertEqual(
                    ds.sizes,
                    {"datetime": 1, "lat": 1, "lon": 2, "lev": 1},
                )
                np.testing.assert_equal(
                    ds.datetime.values[0],
                    np.datetime64("2026-08-12T01:30:00", "ns"),
                )
                self.assertEqual(ds.lat.item(), 40.0)
                self.assertEqual(ds.lon.to_index().tolist(), [-105.0, -104.0])
                self.assertEqual(ds.lev.item(), 850.0)
                self.assertAlmostEqual(
                    float(
                        ds.t2m.sel(
                            datetime="2026-08-12T01:30:00",
                            lat=40,
                            lon=-105,
                        )
                    ),
                    283.0,
                )
                self.assertAlmostEqual(
                    float(
                        ds.rh.sel(
                            datetime="2026-08-12T01:30:00",
                            lev=850,
                            lat=40,
                            lon=-105,
                        )
                    ),
                    40.0,
                )
                self.assertEqual(ds.t2m.attrs["units"], "K")
                self.assertEqual(
                    ds.attrs["geos_fp_collection"],
                    "tavg1_2d_slv_Nx",
                )
                self.assertEqual(ds.attrs["geos_fp_stream"], "assim")
                self.assertEqual(
                    ds.attrs["geos_fp_opendap_url"],
                    client.opendap_url,
                )
            finally:
                ds.close()

    def test_download_file_reuses_cache_unless_overwrite_is_requested(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            client = geos_fp.GeosFp(
                download_dir=tmpdir,
                chunks=None,
            )
            with mock.patch.object(
                client,
                "_open_remote_dataset",
                side_effect=make_dataset,
            ) as open_remote:
                first = client.download_file(
                    "2026-08-12",
                    "2026-08-12",
                    variables="t2m",
                )
                self.assertEqual(open_remote.call_count, 1)

                second = client.download_file(
                    "2026-08-12",
                    "2026-08-12",
                    variables="t2m",
                )
                self.assertEqual(open_remote.call_count, 1)

                overwritten = client.download_file(
                    "2026-08-12",
                    "2026-08-12",
                    variables="t2m",
                    overwrite=True,
                )
                self.assertEqual(open_remote.call_count, 2)

                files = client.download_files(
                    "2026-08-12",
                    "2026-08-12",
                    variables="t2m",
                )

            self.assertIsInstance(first, Path)
            self.assertEqual(first, second)
            self.assertEqual(first, overwritten)
            self.assertEqual(files, [first])
            self.assertTrue(first.is_file())
            with xr.open_dataset(first) as cached:
                self.assertEqual(list(cached.data_vars), ["t2m"])
                self.assertIn("datetime", cached.coords)
                self.assertEqual(cached.sizes["datetime"], 3)

    def test_repeated_download_closes_previous_cached_dataset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            client = geos_fp.GeosFp(download_dir=tmpdir, chunks=None)
            previous = mock.Mock(spec=xr.Dataset)
            client.dataset = previous

            with (
                mock.patch.object(
                    client,
                    "download_file",
                    return_value=Path(tmpdir) / "subset.nc4",
                ),
                mock.patch.object(
                    geos_fp.netcdf4,
                    "module_available",
                    True,
                ),
                mock.patch.object(
                    geos_fp.xr,
                    "open_dataset",
                    return_value=make_dataset(),
                ),
            ):
                client.download("2026-08-12", "2026-08-12")

            previous.close.assert_called_once_with()

    def test_rejects_invalid_requests_before_opening_remote_dataset(self):
        client = geos_fp.GeosFp(chunks=None)
        with mock.patch.object(client, "_open_remote_dataset") as open_remote:
            with self.assertRaisesRegex(
                ValueError,
                "end must be on or after start",
            ):
                client.download("2026-08-13", "2026-08-12")

            with self.assertRaisesRegex(ValueError, "bbox"):
                client.download(
                    "2026-08-12",
                    "2026-08-12",
                    bbox=(-100.0, 40.0, -110.0, 45.0),
                )

            with self.assertRaisesRegex(ValueError, "lev"):
                client.download(
                    "2026-08-12",
                    "2026-08-12",
                    dimensions={"lev": (1000.0, 850.0)},
                )

            open_remote.assert_not_called()

    def test_missing_variable_uses_exact_lowercase_dataset_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            client = geos_fp.GeosFp(
                download_dir=tmpdir,
                chunks=None,
            )
            for variable in ("missing", "T2M"):
                with self.subTest(variable=variable):
                    with (
                        mock.patch.object(
                            client,
                            "_open_remote_dataset",
                            return_value=make_dataset(),
                        ),
                        self.assertRaisesRegex(KeyError, variable),
                    ):
                        client.download_file(
                            "2026-08-12",
                            "2026-08-12",
                            variables=variable,
                        )


class GeosFpValidationTests(unittest.TestCase):
    def test_rejects_invalid_stream_and_assimilation_cycle(self):
        with self.assertRaisesRegex(ValueError, "stream"):
            geos_fp.GeosFp(stream="archive")

        with self.assertRaisesRegex(ValueError, "cycle"):
            geos_fp.GeosFp(
                stream="assim",
                cycle="2026-08-12T06:00:00Z",
            )

    def test_maps_remote_time_chunks_to_cached_datetime(self):
        client = geos_fp.GeosFp(chunks={"time": 24, "lat": 10})

        self.assertEqual(
            client._local_chunks(),
            {"datetime": 24, "lat": 10},
        )


if __name__ == "__main__":
    unittest.main()
