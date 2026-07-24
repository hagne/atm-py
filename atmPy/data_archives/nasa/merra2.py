"""Download MERRA-2 data from NASA Earthdata."""

import configparser
import hashlib
import json
import os
import warnings
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Literal

import xarray as xr

from atmPy import opt_imports


earthaccess = opt_imports.OptionalImport("earthaccess")
harmony = opt_imports.OptionalImport("harmony")

AuthStrategy = Literal[
    "auto",
    "settings",
    "netrc",
    "environment",
    "interactive",
]
DownloadBackend = Literal["harmony", "earthaccess"]


class Merra2:
    """Access and subset a MERRA-2 collection using NASA Earthdata."""

    short_name = "M2T1NXSLV"
    version = "5.12.4"
    ozone_variable = "TO3"
    ozone_units = "Dobsons"
    earthdata_eula_url = (
        "https://urs.earthdata.nasa.gov/users/{username}/unaccepted_eulas"
    )
    settings_template = """\
[earthdata]
username = YOUR_USERNAME
password = YOUR_PASSWORD
"""

    def __init__(
        self,
        *,
        settings_path: str | Path | None = None,
        download_dir: str | Path | None = None,
        chunks: str | dict[str, int] | None = "auto",
        auth_strategy: AuthStrategy = "settings",
        short_name: str = "M2T1NXSLV",
        version: str = "5.12.4",
        backend: DownloadBackend = "harmony",
    ):
        if backend not in ("harmony", "earthaccess"):
            raise ValueError("backend must be 'harmony' or 'earthaccess'")

        self.short_name = short_name
        self.version = version
        self.backend = backend

        if settings_path is None:
            settings_path = Path.home() / ".config" / "atmPy" / "earthdata.ini"
        if download_dir is None:
            download_dir = (
                Path.home() / ".cache" / "atmPy" / "merra2" / self.short_name
            )

        self.settings_path = Path(settings_path).expanduser()
        self.download_dir = Path(download_dir).expanduser()
        self.chunks = chunks
        self.auth_strategy = auth_strategy

        self.username: str | None = None
        self.password: str | None = None
        self.token: str | None = None
        self.auth = None
        self.collection = None
        self.concept_id: str | None = None
        self.harmony_client = None
        self.harmony_request = None
        self.harmony_job_id: str | None = None
        self.granules = None
        self.paths: list[Path] | None = None
        self.dataset: xr.Dataset | None = None

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
                value = datetime.fromisoformat(value.replace("Z", "+00:00"))

        if isinstance(value, date) and not isinstance(value, datetime):
            time = datetime.max.time() if end_of_day else datetime.min.time()
            value = datetime.combine(value, time)

        if not isinstance(value, datetime):
            raise TypeError("start and end must be strings, dates, or datetimes")
        if value.tzinfo is not None:
            value = value.astimezone(timezone.utc).replace(tzinfo=None)
        return value

    def create_settings_file(self) -> Path:
        """Create the Earthdata settings template without overwriting a file."""
        self.settings_path.parent.mkdir(parents=True, exist_ok=True)
        if not self.settings_path.exists():
            self.settings_path.write_text(self.settings_template)
        self.settings_path.chmod(0o600)
        return self.settings_path

    def read_settings(self, *, required: bool = True) -> tuple[str, str] | None:
        """Read credentials from the Earthdata settings file."""
        if not self.settings_path.exists():
            if required:
                self.create_settings_file()
                raise RuntimeError(
                    "NASA Earthdata credentials are not configured. Edit "
                    f"{self.settings_path} and replace YOUR_USERNAME and "
                    "YOUR_PASSWORD, then retry."
                )
            return None

        settings = configparser.ConfigParser(interpolation=None)
        try:
            settings.read(self.settings_path)
            username = settings["earthdata"]["username"].strip()
            password = settings["earthdata"]["password"].strip()
        except (configparser.Error, KeyError) as error:
            raise ValueError(
                f"{self.settings_path} must contain an [earthdata] section "
                "with username and password"
            ) from error

        if username == "YOUR_USERNAME" or password == "YOUR_PASSWORD":
            if required:
                raise RuntimeError(
                    f"Replace YOUR_USERNAME and YOUR_PASSWORD in "
                    f"{self.settings_path}, then retry."
                )
            return None
        if not username or not password:
            raise ValueError(
                f"username and password in {self.settings_path} cannot be empty"
            )

        self.username = username
        self.password = password
        return username, password

    @staticmethod
    def _has_environment_credentials() -> bool:
        return bool(
            os.environ.get("EARTHDATA_TOKEN")
            or (
                os.environ.get("EARTHDATA_USERNAME")
                and os.environ.get("EARTHDATA_PASSWORD")
            )
        )

    def _login_from_settings(self) -> None:
        self.read_settings()
        credential_names = (
            "EARTHDATA_TOKEN",
            "EARTHDATA_USERNAME",
            "EARTHDATA_PASSWORD",
        )
        previous_environment = {
            name: os.environ.get(name) for name in credential_names
        }

        # earthaccess has no public login(username=..., password=...) method.
        # Use its environment strategy, then restore the process environment.
        os.environ.pop("EARTHDATA_TOKEN", None)
        os.environ["EARTHDATA_USERNAME"] = self.username
        os.environ["EARTHDATA_PASSWORD"] = self.password
        try:
            self.auth = earthaccess.login(strategy="environment")
        finally:
            for name, value in previous_environment.items():
                if value is None:
                    os.environ.pop(name, None)
                else:
                    os.environ[name] = value

    def login(self) -> None:
        """Authenticate with NASA Earthdata using the configured strategy."""
        strategy = self.auth_strategy
        if strategy == "auto":
            if self._has_environment_credentials():
                strategy = "environment"
            elif self.read_settings(required=False):
                strategy = "settings"
            else:
                strategy = "netrc"

        try:
            if strategy == "settings":
                self._login_from_settings()
            else:
                self.auth = earthaccess.login(strategy=strategy)
                if strategy == "environment":
                    self.username = os.environ.get("EARTHDATA_USERNAME")
                    self.password = os.environ.get("EARTHDATA_PASSWORD")
                    self.token = os.environ.get("EARTHDATA_TOKEN")
                elif self.auth is not None:
                    self.username = getattr(self.auth, "username", None)
        except earthaccess.exceptions.LoginStrategyUnavailable as error:
            self.create_settings_file()
            raise RuntimeError(
                "NASA Earthdata credentials were not found.\n\n"
                f"A settings template was created at {self.settings_path}.\n"
                "Edit it to contain:\n"
                "[earthdata]\n"
                "username = YOUR_USERNAME\n"
                "password = YOUR_PASSWORD\n\n"
                "Then rerun the download."
            ) from error

    @property
    def eula_url(self) -> str:
        username = self.username or "earthaccess"
        return self.earthdata_eula_url.format(username=username)

    def _find_collection(self) -> None:
        collections = earthaccess.search_datasets(
            short_name=self.short_name,
            version=self.version,
        )
        if not collections:
            raise FileNotFoundError(
                f"No MERRA-2 collection found for {self.short_name} "
                f"version {self.version}"
            )

        self.collection = collections[0]
        concept_id = self.collection.concept_id
        if callable(concept_id):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", FutureWarning)
                concept_id = concept_id()
        self.concept_id = concept_id

    def _harmony_download_dir(
        self,
        *,
        variables: list[str],
        start: datetime,
        end: datetime,
        bbox: tuple[float, float, float, float] | None,
        dimensions: dict[str, tuple[float, float]] | None,
    ) -> Path:
        request_description = {
            "short_name": self.short_name,
            "version": self.version,
            "variables": sorted(variables),
            "start": start.isoformat(),
            "end": end.isoformat(),
            "bbox": bbox,
            "dimensions": dimensions,
        }
        request_json = json.dumps(request_description, sort_keys=True)
        request_id = hashlib.sha256(request_json.encode()).hexdigest()[:12]
        return self.download_dir / "harmony" / request_id

    def _download_harmony(
        self,
        *,
        variables: list[str],
        start: datetime,
        end: datetime,
        bbox: tuple[float, float, float, float] | None,
        dimensions: dict[str, tuple[float, float]] | None,
        show_progress: bool,
    ) -> None:
        try:
            client_class = harmony.Client
        except ImportError as error:
            raise ImportError(
                "Harmony subsetting requires harmony-py. Install the NASA "
                'dependencies with: pip install -e ".[nasa]"'
            ) from error

        self._find_collection()

        client_kwargs = {}
        if self.token:
            client_kwargs["token"] = self.token
        elif self.username and self.password:
            client_kwargs["auth"] = (self.username, self.password)
        self.harmony_client = client_class(**client_kwargs)

        spatial = harmony.BBox(*bbox) if bbox is not None else None
        harmony_dimensions = None
        if dimensions:
            harmony_dimensions = [
                harmony.Dimension(name=name, min=bounds[0], max=bounds[1])
                for name, bounds in dimensions.items()
            ]

        self.harmony_request = harmony.Request(
            collection=harmony.Collection(self.concept_id),
            spatial=spatial,
            temporal={"start": start, "stop": end},
            dimensions=harmony_dimensions,
            variables=variables,
        )
        self.harmony_job_id = self.harmony_client.submit(self.harmony_request)
        try:
            self.harmony_client.wait_for_processing(
                self.harmony_job_id,
                show_progress=show_progress,
            )
        except harmony.client.ProcessingFailedException as error:
            message = str(error)
            if "403 error retrieving" not in message:
                raise
            raise RuntimeError(
                "NASA Harmony accepted the subset request, but its OPeNDAP "
                "subsetter could not access the source granule (HTTP 403). "
                "This is not an xarray error.\n\n"
                "First confirm that all agreements listed here are accepted:\n"
                f"{self.eula_url}\n\n"
                "If they are already accepted, NASA's Harmony/OPeNDAP service "
                "is unavailable for this granule. Retry later, or download the "
                "complete original granule by passing backend='earthaccess'. "
                "atmPy does not switch automatically because a daily MERRA-2 "
                "file can be hundreds of MB.\n\n"
                f"Harmony job: {self.harmony_job_id}"
            ) from error

        subset_dir = self._harmony_download_dir(
            variables=variables,
            start=start,
            end=end,
            bbox=bbox,
            dimensions=dimensions,
        )
        subset_dir.mkdir(parents=True, exist_ok=True)
        downloads = self.harmony_client.download_all(
            self.harmony_job_id,
            directory=str(subset_dir),
            overwrite=False,
        )
        self.paths = [Path(download.result()) for download in downloads]

    def _download_granules(
        self,
        *,
        start: datetime,
        end: datetime,
    ) -> None:
        temporal = (
            f"{start.isoformat()}Z",
            f"{end.isoformat()}Z",
        )
        self.granules = earthaccess.search_data(
            short_name=self.short_name,
            version=self.version,
            temporal=temporal,
        )
        if not self.granules:
            raise FileNotFoundError(
                f"No {self.short_name} granules found between "
                f"{temporal[0]} and {temporal[1]}"
            )

        try:
            self.paths = earthaccess.download(
                self.granules,
                self.download_dir,
            )
        except earthaccess.exceptions.EulaNotAccepted as error:
            raise RuntimeError(
                "NASA Earthdata login succeeded, but your account has not "
                "accepted a required GES DISC data-use agreement. Open this "
                "page in a web browser, sign in, and accept the listed "
                f"agreement:\n{self.eula_url}\n"
                "Then rerun the download. atmPy cannot accept a legal "
                "agreement on your behalf."
            ) from error

    @staticmethod
    def _coordinate_slice(
        ds: xr.Dataset,
        coordinate: str,
        lower: float,
        upper: float,
    ) -> slice:
        if coordinate not in ds.coords:
            raise ValueError(
                f"MERRA-2 dataset does not contain the requested "
                f"{coordinate!r} coordinate"
            )
        index = ds.indexes.get(coordinate)
        if index is not None and index.is_monotonic_decreasing:
            return slice(upper, lower)
        return slice(lower, upper)

    def download(
        self,
        start: str | date | datetime,
        end: str | date | datetime,
        *,
        variables: str | list[str] | tuple[str, ...],
        bbox: tuple[float, float, float, float] | None = None,
        dimensions: dict[str, tuple[float, float]] | None = None,
        backend: DownloadBackend | None = None,
        login: bool = True,
        show_progress: bool = True,
    ) -> xr.Dataset:
        """Download a variable subset from a MERRA-2 collection.

        ``bbox`` must be ordered as west, south, east, north. ``dimensions``
        maps dimension names to inclusive minimum and maximum values. Harmony
        performs these subsets before download. Use ``backend="earthaccess"``
        to download complete granules and subset them locally instead.
        """
        start_datetime = self._as_utc_naive(start)
        end_datetime = self._as_utc_naive(end, end_of_day=True)
        if end_datetime < start_datetime:
            raise ValueError("end must be on or after start")

        if isinstance(variables, str):
            variable_names = [variables]
        else:
            variable_names = list(variables)
        if not variable_names or not all(
            isinstance(variable, str) and variable
            for variable in variable_names
        ):
            raise ValueError("variables must contain at least one variable name")

        if bbox is not None:
            if len(bbox) != 4:
                raise ValueError("bbox must contain west, south, east, north")
            west, south, east, north = bbox
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

        if dimensions:
            for name, bounds in dimensions.items():
                if len(bounds) != 2 or bounds[0] > bounds[1]:
                    raise ValueError(
                        f"dimension {name!r} must have (minimum, maximum) bounds"
                    )

        selected_backend = backend or self.backend
        if selected_backend not in ("harmony", "earthaccess"):
            raise ValueError("backend must be 'harmony' or 'earthaccess'")

        self.download_dir.mkdir(parents=True, exist_ok=True)
        if login:
            self.login()

        if selected_backend == "harmony":
            self._download_harmony(
                variables=variable_names,
                start=start_datetime,
                end=end_datetime,
                bbox=bbox,
                dimensions=dimensions,
                show_progress=show_progress,
            )
        else:
            self._download_granules(
                start=start_datetime,
                end=end_datetime,
            )
        if not self.paths:
            raise RuntimeError("NASA Earthdata did not return any downloaded files")

        ds = xr.open_mfdataset(
            self.paths,
            combine="by_coords",
            chunks=self.chunks,
        )
        missing_variables = [
            variable for variable in variable_names if variable not in ds
        ]
        if missing_variables:
            ds.close()
            raise KeyError(
                f"Variables not found in {self.short_name}: "
                f"{', '.join(missing_variables)}"
            )
        ds = ds[variable_names]
        if "time" not in ds.coords:
            ds.close()
            raise ValueError(
                "MERRA-2 dataset does not contain the expected time coordinate"
            )
        ds = ds.rename(time="datetime")

        if selected_backend == "earthaccess":
            if bbox is not None:
                west, south, east, north = bbox
                ds = ds.sel(
                    lat=self._coordinate_slice(ds, "lat", south, north),
                    lon=self._coordinate_slice(ds, "lon", west, east),
                )
            if dimensions:
                ds = ds.sel(
                    {
                        name: self._coordinate_slice(ds, name, *bounds)
                        for name, bounds in dimensions.items()
                    }
                )

        self.dataset = ds.sel(datetime=slice(start_datetime, end_datetime))
        return self.dataset

    def download_total_column_ozone(
        self,
        start: str | date | datetime,
        end: str | date | datetime,
        *,
        bbox: tuple[float, float, float, float] | None = None,
        backend: DownloadBackend | None = None,
        login: bool = True,
        show_progress: bool = True,
    ) -> xr.Dataset:
        """Download hourly total-column ozone for an inclusive time period."""
        ds = self.download(
            start,
            end,
            variables=self.ozone_variable,
            bbox=bbox,
            backend=backend,
            login=login,
            show_progress=show_progress,
        )
        units = ds[self.ozone_variable].attrs.get("units")
        if units != self.ozone_units:
            ds.close()
            self.dataset = None
            raise ValueError(
                f"{self.ozone_variable} units must be "
                f"{self.ozone_units!r}, found {units!r}"
            )

        return ds


def create_earthdata_settings_file(path: str | Path | None = None) -> Path:
    """Create the settings file used by :class:`Merra2`."""
    return Merra2(settings_path=path).create_settings_file()


def download_total_column_ozone(
    start: str | date | datetime,
    end: str | date | datetime,
    *,
    download_dir: str | Path | None = None,
    chunks: str | dict[str, int] | None = "auto",
    login: bool = True,
    auth_strategy: AuthStrategy = "auto",
    settings_path: str | Path | None = None,
    bbox: tuple[float, float, float, float] | None = None,
    backend: DownloadBackend = "harmony",
    show_progress: bool = True,
) -> xr.Dataset:
    """Compatibility wrapper for :meth:`Merra2.download_total_column_ozone`."""
    merra2 = Merra2(
        settings_path=settings_path,
        download_dir=download_dir,
        chunks=chunks,
        auth_strategy=auth_strategy,
        backend=backend,
    )
    return merra2.download_total_column_ozone(
        start,
        end,
        bbox=bbox,
        login=login,
        show_progress=show_progress,
    )
