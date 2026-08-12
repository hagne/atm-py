"""Download GEOS Forward Processing data from NCCS OPeNDAP.

GEOS-FP is distributed by NASA's Global Modeling and Assimilation Office
through a public GrADS Data Server at the NASA Center for Climate Simulation.
Unlike the Earthdata products handled by :mod:`.earthdata`, this service does
not use CMR collection metadata or Earthdata authentication.  It exposes large
virtual datasets that can be subset remotely before they are saved locally.
"""

import contextlib
import hashlib
import http.client
import http.server
import json
import re
import socket
import ssl
import tempfile
import threading
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Literal
from urllib.parse import urlsplit, urlunsplit

import xarray as xr

from atmPy import opt_imports


netcdf4 = opt_imports.OptionalImport("netCDF4")

GeosFpStream = Literal["assim", "fcast", "seamless"]

_HOP_BY_HOP_HEADERS = {
    "connection",
    "keep-alive",
    "proxy-authenticate",
    "proxy-authorization",
    "te",
    "trailer",
    "transfer-encoding",
    "upgrade",
}


class _IPv4HTTPSConnection(http.client.HTTPSConnection):
    """HTTPS connection whose upstream socket is restricted to IPv4."""

    def connect(self):
        last_error = None
        addresses = socket.getaddrinfo(
            self.host,
            self.port,
            socket.AF_INET,
            socket.SOCK_STREAM,
        )
        for family, socket_type, protocol, _, address in addresses:
            upstream_socket = socket.socket(family, socket_type, protocol)
            upstream_socket.settimeout(self.timeout)
            try:
                if self.source_address:
                    upstream_socket.bind(self.source_address)
                upstream_socket.connect(address)
                self.sock = self._context.wrap_socket(
                    upstream_socket,
                    server_hostname=self.host,
                )
                return
            except OSError as error:
                last_error = error
                upstream_socket.close()

        if last_error is not None:
            raise last_error
        raise OSError(f"No IPv4 address found for {self.host}")


class _IPv4ProxyServer(http.server.ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, server_address, handler, target_url: str):
        super().__init__(server_address, handler)
        target = urlsplit(target_url)
        self.target_host = target.hostname
        self.target_port = target.port or 443
        self.ssl_context = ssl.create_default_context()


class _IPv4ProxyHandler(http.server.BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def do_GET(self):
        response_started = False
        connection = _IPv4HTTPSConnection(
            self.server.target_host,
            self.server.target_port,
            timeout=120,
            context=self.server.ssl_context,
        )
        try:
            headers = {
                name: value
                for name, value in self.headers.items()
                if name.lower() not in _HOP_BY_HOP_HEADERS | {"host"}
            }
            headers["Host"] = self.server.target_host
            headers["Connection"] = "close"
            connection.request("GET", self.path, headers=headers)
            response = connection.getresponse()

            self.send_response_only(response.status, response.reason)
            response_started = True
            for name, value in response.getheaders():
                if name.lower() not in _HOP_BY_HOP_HEADERS:
                    self.send_header(name, value)
            self.send_header("Connection", "close")
            self.end_headers()
            while data := response.read(1024 * 1024):
                self.wfile.write(data)
        except OSError as error:
            if not response_started:
                self.send_error(502, explain=str(error))
        finally:
            connection.close()
            self.close_connection = True

    def log_message(self, format, *args):
        pass


@contextlib.contextmanager
def _nccs_ipv4_proxy(opendap_url: str):
    """Expose an NCCS URL on loopback while connecting upstream over IPv4."""
    target = urlsplit(opendap_url)
    server = _IPv4ProxyServer(
        ("127.0.0.1", 0),
        _IPv4ProxyHandler,
        opendap_url,
    )
    thread = threading.Thread(
        target=server.serve_forever,
        kwargs={"poll_interval": 0.05},
        daemon=True,
    )
    thread.start()
    try:
        local_url = urlunsplit(
            (
                "http",
                f"127.0.0.1:{server.server_port}",
                target.path,
                target.query,
                target.fragment,
            )
        )
        yield local_url
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


class GeosFp:
    """Access GEOS-FP near-real-time data through NCCS OPeNDAP.

    Parameters
    ----------
    collection
        GEOS-FP file collection, for example ``"tavg1_2d_slv_Nx"`` for
        hourly, time-averaged single-level diagnostics or
        ``"inst3_3d_asm_Np"`` for three-hourly pressure-level assimilated
        meteorology.  NCCS OPeNDAP variable names are lowercase.
    stream
        ``"assim"`` accesses the continuously growing assimilation history.
        ``"fcast"`` accesses a forecast cycle, and ``"seamless"`` accesses
        the NCCS analysis-plus-forecast dataset available for some collections.
    cycle
        Forecast or seamless cycle as a datetime, ISO datetime string, or
        ``"YYYYMMDD_HH"``.  By default, the server's ``latest`` alias is used.
        Assimilation datasets do not take a cycle.
    download_dir
        Directory for cached NetCDF subsets.  The default is
        ``~/.cache/atmPy/geos_fp/<collection>``.
    chunks
        Chunking passed to xarray when opening the remote and cached datasets.
        The default, ``None``, avoids requiring dask.  A ``time`` key in a
        chunk dictionary is mapped to ``datetime`` for the cached dataset.
    force_ipv4
        Route OPeNDAP requests through an ephemeral loopback proxy whose
        upstream connection is restricted to IPv4.  This works around networks
        on which the NCCS IPv6 endpoint accepts a connection but stalls during
        TLS negotiation.  TLS certificate verification remains enabled between
        the proxy and NCCS.

    Notes
    -----
    The official collection catalogs are available at ``catalog_url``.  The
    ``latest`` forecast and seamless aliases change as GEOS-FP advances.  Pass
    ``overwrite=True`` to :meth:`download` or :meth:`download_file` to refresh
    a cached request made through one of these aliases.

    Following the GEOS-5.43.0 transition on 2026-02-26, NCCS OPeNDAP cannot
    serve older assimilation data for ``inst3_3d_aer_Nv``,
    ``tavg3_2d_adg_Nx``, or ``tavg3_2d_aer_Nx``.  Forecast data are unaffected.
    """

    base_url = "https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg"
    project_url = (
        "https://gmao.gsfc.nasa.gov/gmao-products/"
        "geos-near-real-time-data-products/"
    )
    documentation_url = "https://www.nccs.nasa.gov/opendap-servers/"

    def __init__(
        self,
        *,
        collection: str = "tavg1_2d_slv_Nx",
        stream: GeosFpStream = "assim",
        cycle: str | date | datetime | None = None,
        download_dir: str | Path | None = None,
        chunks: str | dict[str, int] | None = None,
        force_ipv4: bool = False,
    ):
        if not isinstance(collection, str) or not re.fullmatch(
            r"[A-Za-z0-9_-]+",
            collection,
        ):
            raise ValueError(
                "collection must be a GEOS-FP collection name containing "
                "only letters, numbers, underscores, or hyphens"
            )
        if stream not in ("assim", "fcast", "seamless"):
            raise ValueError("stream must be 'assim', 'fcast', or 'seamless'")
        if stream == "assim" and cycle is not None:
            raise ValueError("cycle is only valid for fcast or seamless streams")

        self.collection = collection
        self.stream = stream
        self.cycle = self._normalize_cycle(cycle)
        if download_dir is None:
            download_dir = (
                Path.home() / ".cache" / "atmPy" / "geos_fp" / collection
            )
        self.download_dir = Path(download_dir).expanduser()
        self.chunks = chunks
        self.force_ipv4 = force_ipv4

        self.path: Path | None = None
        self.paths: list[Path] | None = None
        self.dataset: xr.Dataset | None = None

    @staticmethod
    def _normalize_cycle(value: str | date | datetime | None) -> str | None:
        if value is None or value == "latest":
            return None
        if isinstance(value, str):
            try:
                cycle_datetime = datetime.strptime(value, "%Y%m%d_%H")
            except ValueError:
                try:
                    cycle_datetime = datetime.fromisoformat(
                        value.replace("Z", "+00:00")
                    )
                except ValueError as error:
                    raise ValueError(
                        "cycle must be an ISO datetime or 'YYYYMMDD_HH'"
                    ) from error
        elif isinstance(value, date) and not isinstance(value, datetime):
            cycle_datetime = datetime.combine(value, datetime.min.time())
        elif isinstance(value, datetime):
            cycle_datetime = value
        else:
            raise TypeError("cycle must be a string, date, datetime, or None")

        if cycle_datetime.tzinfo is not None:
            cycle_datetime = cycle_datetime.astimezone(timezone.utc).replace(
                tzinfo=None
            )
        if any(
            (
                cycle_datetime.minute,
                cycle_datetime.second,
                cycle_datetime.microsecond,
            )
        ):
            raise ValueError("cycle must fall on a whole UTC hour")
        return cycle_datetime.strftime("%Y%m%d_%H")

    @property
    def catalog_url(self) -> str:
        """Return the live NCCS catalog URL for the selected stream."""
        return f"{self.base_url}/{self.stream}"

    @property
    def opendap_url(self) -> str:
        """Return the NCCS OPeNDAP dataset URL selected by this instance."""
        base = f"{self.catalog_url}/{self.collection}"
        if self.stream == "assim":
            return base
        if self.cycle is None:
            return f"{base}.latest"
        return f"{base}/{self.collection}.{self.cycle}"

    @staticmethod
    def _as_utc_naive(
        value: str | date | datetime,
        *,
        end_of_day: bool = False,
    ) -> datetime:
        if isinstance(value, str):
            try:
                value = date.fromisoformat(value)
            except ValueError:
                try:
                    value = datetime.fromisoformat(value.replace("Z", "+00:00"))
                except ValueError as error:
                    raise ValueError(
                        "start and end must be ISO dates or datetimes"
                    ) from error

        if isinstance(value, date) and not isinstance(value, datetime):
            time = datetime.max.time() if end_of_day else datetime.min.time()
            value = datetime.combine(value, time)
        if not isinstance(value, datetime):
            raise TypeError("start and end must be strings, dates, or datetimes")
        if value.tzinfo is not None:
            value = value.astimezone(timezone.utc).replace(tzinfo=None)
        return value

    @staticmethod
    def _coordinate_slice(
        ds: xr.Dataset,
        coordinate: str,
        lower,
        upper,
    ) -> slice:
        if coordinate not in ds.coords:
            raise ValueError(
                "GEOS-FP dataset does not contain the requested "
                f"{coordinate!r} coordinate"
            )
        index = ds.indexes.get(coordinate)
        if index is not None and index.is_monotonic_decreasing:
            return slice(upper, lower)
        return slice(lower, upper)

    def _validate_download_request(
        self,
        start: str | date | datetime,
        end: str | date | datetime,
        *,
        variables: str | list[str] | tuple[str, ...] | None,
        bbox: tuple[float, float, float, float] | None,
        dimensions: dict[str, tuple[float, float]] | None,
    ) -> tuple[
        datetime,
        datetime,
        list[str] | None,
        tuple[float, float, float, float] | None,
        dict[str, tuple[float, float]] | None,
    ]:
        start_datetime = self._as_utc_naive(start)
        end_datetime = self._as_utc_naive(end, end_of_day=True)
        if end_datetime < start_datetime:
            raise ValueError("end must be on or after start")

        if variables is None:
            variable_names = None
        elif isinstance(variables, str):
            variable_names = [variables]
        else:
            variable_names = list(variables)
        if variable_names is not None and (
            not variable_names
            or not all(
                isinstance(variable, str) and variable
                for variable in variable_names
            )
        ):
            raise ValueError("variables must contain at least one variable name")
        if variable_names is not None and len(variable_names) != len(
            set(variable_names)
        ):
            raise ValueError("variables cannot contain duplicate names")

        normalized_bbox = None
        if bbox is not None:
            if len(bbox) != 4:
                raise ValueError("bbox must contain west, south, east, north")
            try:
                normalized_bbox = tuple(float(value) for value in bbox)
            except (TypeError, ValueError) as error:
                raise ValueError("bbox bounds must be numbers") from error
            west, south, east, north = normalized_bbox
            if not (-180 <= west < east <= 180):
                raise ValueError(
                    "bbox west and east must satisfy "
                    "-180 <= west < east <= 180"
                )
            if not (-90 <= south < north <= 90):
                raise ValueError(
                    "bbox south and north must satisfy "
                    "-90 <= south < north <= 90"
                )

        normalized_dimensions = None
        if dimensions:
            normalized_dimensions = {}
            for name, bounds in dimensions.items():
                if not isinstance(name, str) or not name:
                    raise ValueError("dimension names must be non-empty strings")
                if name in ("time", "datetime", "lat", "lon"):
                    raise ValueError(
                        f"Use start/end or bbox instead of dimensions[{name!r}]"
                    )
                if len(bounds) != 2:
                    raise ValueError(
                        f"dimension {name!r} must have (minimum, maximum) bounds"
                    )
                try:
                    lower, upper = (float(value) for value in bounds)
                except (TypeError, ValueError) as error:
                    raise ValueError(
                        f"dimension {name!r} bounds must be numbers"
                    ) from error
                if lower > upper:
                    raise ValueError(
                        f"dimension {name!r} must have (minimum, maximum) bounds"
                    )
                normalized_dimensions[name] = (lower, upper)

        return (
            start_datetime,
            end_datetime,
            variable_names,
            normalized_bbox,
            normalized_dimensions,
        )

    def _cache_path(
        self,
        *,
        start: datetime,
        end: datetime,
        variables: list[str] | None,
        bbox: tuple[float, float, float, float] | None,
        dimensions: dict[str, tuple[float, float]] | None,
    ) -> Path:
        request = {
            "schema": 1,
            "opendap_url": self.opendap_url,
            "start": start.isoformat(),
            "end": end.isoformat(),
            "variables": variables,
            "bbox": bbox,
            "dimensions": dimensions,
        }
        request_json = json.dumps(request, sort_keys=True)
        request_id = hashlib.sha256(request_json.encode()).hexdigest()[:12]
        period = f"{start:%Y%m%dT%H%M%S}_{end:%Y%m%dT%H%M%S}"
        return self.download_dir / self.stream / f"{period}_{request_id}.nc4"

    @staticmethod
    def _require_netcdf4() -> None:
        if not netcdf4.module_available:
            raise ImportError(
                "NCCS OPeNDAP access requires netCDF4. Install the NASA "
                'dependencies with: pip install -e ".[nasa]"'
            )

    def _open_remote_dataset(
        self,
        opendap_url: str | None = None,
    ) -> xr.Dataset:
        self._require_netcdf4()
        if opendap_url is None:
            opendap_url = self.opendap_url
        try:
            return xr.open_dataset(
                opendap_url,
                engine="netcdf4",
                chunks=self.chunks,
            )
        except OSError as error:
            ipv4_hint = ""
            if not self.force_ipv4:
                ipv4_hint = (
                    " If this is a curl timeout on a host with working IPv4, "
                    "retry with GeosFp(force_ipv4=True)."
                )
            raise OSError(
                f"Could not open GEOS-FP collection {self.collection!r} from "
                f"{self.opendap_url}. Check the collection, stream, and cycle "
                f"against {self.catalog_url}, and check NCCS/network "
                f"connectivity. netCDF4 reported: {error}.{ipv4_hint}"
            ) from error

    @contextlib.contextmanager
    def _remote_dataset(self):
        if not self.force_ipv4:
            with self._open_remote_dataset() as remote:
                yield remote
            return

        with _nccs_ipv4_proxy(self.opendap_url) as local_url:
            with self._open_remote_dataset(local_url) as remote:
                yield remote

    def _local_chunks(self) -> str | dict[str, int] | None:
        if not isinstance(self.chunks, dict) or "time" not in self.chunks:
            return self.chunks
        chunks = dict(self.chunks)
        chunks["datetime"] = chunks.pop("time")
        return chunks

    @staticmethod
    def _ensure_nonempty(ds: xr.Dataset, coordinate: str) -> None:
        if coordinate in ds.dims and ds.sizes[coordinate] == 0:
            raise FileNotFoundError(
                f"The requested GEOS-FP subset contains no {coordinate} values"
            )

    def _write_subset(
        self,
        path: Path,
        *,
        start: datetime,
        end: datetime,
        variables: list[str] | None,
        bbox: tuple[float, float, float, float] | None,
        dimensions: dict[str, tuple[float, float]] | None,
    ) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary_file = tempfile.NamedTemporaryFile(
            prefix=f".{path.stem}-",
            suffix=path.suffix,
            dir=path.parent,
            delete=False,
        )
        temporary_path = Path(temporary_file.name)
        temporary_file.close()

        try:
            with self._remote_dataset() as remote:
                if "time" not in remote.coords:
                    raise ValueError(
                        "GEOS-FP dataset does not contain the expected time "
                        "coordinate"
                    )
                if variables is not None:
                    missing_variables = [
                        variable
                        for variable in variables
                        if variable not in remote.data_vars
                    ]
                    if missing_variables:
                        raise KeyError(
                            f"Variables not found in {self.collection}: "
                            f"{', '.join(missing_variables)}"
                        )
                    subset = remote[variables]
                else:
                    subset = remote

                subset = subset.sel(
                    time=self._coordinate_slice(remote, "time", start, end)
                )
                self._ensure_nonempty(subset, "time")

                if bbox is not None:
                    west, south, east, north = bbox
                    subset = subset.sel(
                        lat=self._coordinate_slice(
                            subset,
                            "lat",
                            south,
                            north,
                        ),
                        lon=self._coordinate_slice(
                            subset,
                            "lon",
                            west,
                            east,
                        ),
                    )
                    self._ensure_nonempty(subset, "lat")
                    self._ensure_nonempty(subset, "lon")

                if dimensions:
                    subset = subset.sel(
                        {
                            name: self._coordinate_slice(subset, name, *bounds)
                            for name, bounds in dimensions.items()
                        }
                    )
                    for name in dimensions:
                        self._ensure_nonempty(subset, name)

                subset = subset.rename(time="datetime")
                subset.attrs.update(
                    {
                        "geos_fp_collection": self.collection,
                        "geos_fp_stream": self.stream,
                        "geos_fp_opendap_url": self.opendap_url,
                    }
                )
                if self.cycle is not None:
                    subset.attrs["geos_fp_cycle"] = self.cycle
                subset.to_netcdf(temporary_path, engine="netcdf4")

            temporary_path.replace(path)
        finally:
            if temporary_path.exists():
                temporary_path.unlink()

    def download_file(
        self,
        start: str | date | datetime,
        end: str | date | datetime,
        *,
        variables: str | list[str] | tuple[str, ...] | None = None,
        bbox: tuple[float, float, float, float] | None = None,
        dimensions: dict[str, tuple[float, float]] | None = None,
        overwrite: bool = False,
    ) -> Path:
        """Download an OPeNDAP subset and return its cached NetCDF path.

        ``bbox`` is ordered west, south, east, north.  ``dimensions`` maps
        other exact coordinate names, such as ``lev``, to inclusive minimum
        and maximum values.  Remote selection occurs before values are
        transferred from NCCS.
        """
        (
            start_datetime,
            end_datetime,
            variable_names,
            normalized_bbox,
            normalized_dimensions,
        ) = self._validate_download_request(
            start,
            end,
            variables=variables,
            bbox=bbox,
            dimensions=dimensions,
        )
        path = self._cache_path(
            start=start_datetime,
            end=end_datetime,
            variables=variable_names,
            bbox=normalized_bbox,
            dimensions=normalized_dimensions,
        )
        if overwrite or not path.exists():
            self._write_subset(
                path,
                start=start_datetime,
                end=end_datetime,
                variables=variable_names,
                bbox=normalized_bbox,
                dimensions=normalized_dimensions,
            )

        self.path = path
        self.paths = [path]
        return path

    def download_files(
        self,
        start: str | date | datetime,
        end: str | date | datetime,
        *,
        variables: str | list[str] | tuple[str, ...] | None = None,
        bbox: tuple[float, float, float, float] | None = None,
        dimensions: dict[str, tuple[float, float]] | None = None,
        overwrite: bool = False,
    ) -> list[Path]:
        """Download an OPeNDAP subset and return it as a one-item path list."""
        return [
            self.download_file(
                start,
                end,
                variables=variables,
                bbox=bbox,
                dimensions=dimensions,
                overwrite=overwrite,
            )
        ]

    def download(
        self,
        start: str | date | datetime,
        end: str | date | datetime,
        *,
        variables: str | list[str] | tuple[str, ...] | None = None,
        bbox: tuple[float, float, float, float] | None = None,
        dimensions: dict[str, tuple[float, float]] | None = None,
        overwrite: bool = False,
    ) -> xr.Dataset:
        """Download a GEOS-FP subset and open the cached xarray dataset."""
        path = self.download_file(
            start,
            end,
            variables=variables,
            bbox=bbox,
            dimensions=dimensions,
            overwrite=overwrite,
        )
        self._require_netcdf4()
        if self.dataset is not None:
            self.dataset.close()
        self.dataset = xr.open_dataset(
            path,
            engine="netcdf4",
            chunks=self._local_chunks(),
        )
        return self.dataset
