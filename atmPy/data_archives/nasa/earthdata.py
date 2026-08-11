"""Common access to data products distributed through NASA Earthdata."""

import configparser
import hashlib
import json
import os
import re
import warnings
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Literal

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


class Earthdata:
    """Authenticate with NASA Earthdata and retrieve collection files.

    This class handles the common Earthdata interfaces: authentication, CMR
    collection resolution, Harmony subsetting, and original-granule downloads.
    Product-specific subclasses are responsible for interpreting the downloaded
    files and applying their scientific coordinate and metadata conventions.

    Parameters
    ----------
    short_name
        NASA collection short name. This may be omitted when the instance is
        used only to manage Earthdata credentials.
    provider
        Optional CMR collection provider, for example ``"GES_DISC"``.
    version
        Collection version. By default, CMR's newest current version is used.
    backend
        ``"harmony"`` requests remote variable and spatial subsetting.
        ``"earthaccess"`` downloads original granules.
    """

    current_collection_states = frozenset({"ACTIVE", "COMPLETE"})
    require_current_version = False
    cache_namespace = "earthdata"
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
        auth_strategy: AuthStrategy = "settings",
        short_name: str | None = None,
        provider: str | None = None,
        version: str | None = None,
        backend: DownloadBackend = "harmony",
    ):
        if backend not in ("harmony", "earthaccess"):
            raise ValueError("backend must be 'harmony' or 'earthaccess'")
        if short_name is not None and (
            not isinstance(short_name, str) or not short_name.strip()
        ):
            raise ValueError("short_name must be a NASA collection short name")
        if provider is not None and (
            not isinstance(provider, str) or not provider.strip()
        ):
            raise ValueError("provider must be a non-empty string or None")
        if version is not None and (
            not isinstance(version, str) or not version.strip()
        ):
            raise ValueError("version must be a non-empty string or None")

        self.short_name = short_name.strip() if short_name is not None else None
        self.provider = provider.strip() if provider is not None else None
        self.requested_version = version.strip() if version is not None else None
        self.version = self.requested_version
        self.backend = backend

        if settings_path is None:
            settings_path = Path.home() / ".config" / "atmPy" / "earthdata.ini"
        if download_dir is None:
            download_dir = Path.home() / ".cache" / "atmPy" / self.cache_namespace
            if self.short_name is not None:
                download_dir /= self.short_name

        self.settings_path = Path(settings_path).expanduser()
        self.download_dir = Path(download_dir).expanduser()
        self.auth_strategy = auth_strategy

        self.username: str | None = None
        self.password: str | None = None
        self.token: str | None = None
        self.auth = None
        self.collection = None
        self.concept_id: str | None = None
        self.collection_status: str | None = None
        self.product_title: str | None = None
        self.native_format: str | None = None
        self.dataset_url: str | None = None
        self.data_url: str | None = None
        self.product_metadata: dict[str, str | bool] | None = None
        self.harmony_client = None
        self.harmony_request = None
        self.harmony_job_id: str | None = None
        self.granules = None
        self.paths: list[Path] | None = None

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

    @staticmethod
    def _collection_umm(collection) -> dict:
        return collection.get("umm", {})

    @classmethod
    def _collection_matches(cls, collection) -> bool:
        return True

    @classmethod
    def _version_key(cls, collection) -> tuple[int, ...]:
        version = cls._collection_umm(collection).get("Version", "")
        return tuple(int(part) for part in re.findall(r"\d+", version))

    @classmethod
    def _current_collection(cls, collections: list):
        candidates = [
            collection
            for collection in collections
            if cls._collection_matches(collection)
            and cls._collection_umm(collection).get("CollectionProgress")
            in cls.current_collection_states
        ]
        if not candidates:
            return None

        return max(
            candidates,
            key=lambda collection: (
                cls._version_key(collection),
                cls._collection_umm(collection).get("CollectionProgress")
                == "ACTIVE",
                collection.get("meta", {}).get("revision-date", ""),
            ),
        )

    @classmethod
    def _collection_url(cls, collection, url_type: str) -> str | None:
        related_urls = cls._collection_umm(collection).get("RelatedUrls", [])
        for related_url in related_urls:
            if related_url.get("Type") == url_type:
                return related_url.get("URL")
        return None

    @classmethod
    def _metadata_from_collection(
        cls,
        collection,
    ) -> dict[str, str | bool]:
        umm = cls._collection_umm(collection)
        meta = collection.get("meta", {})
        formats = {
            file_information.get("Format")
            for file_information in umm.get(
                "ArchiveAndDistributionInformation",
                {},
            ).get("FileArchiveInformation", [])
            if file_information.get("Format")
        }
        associations = meta.get("associations", {})
        return {
            "short_name": umm.get("ShortName", ""),
            "version": umm.get("Version", ""),
            "status": umm.get("CollectionProgress", ""),
            "title": umm.get("EntryTitle", ""),
            "native_format": ", ".join(sorted(formats)),
            "concept_id": meta.get("concept-id", ""),
            "dataset_url": cls._collection_url(
                collection,
                "DATA SET LANDING PAGE",
            )
            or "",
            "data_url": cls._collection_url(collection, "GET DATA") or "",
            "associated_services": bool(associations.get("services")),
        }

    @staticmethod
    def _search_datasets(**kwargs):
        return earthaccess.search_datasets(**kwargs)

    def _collection_search_kwargs(self) -> dict[str, str | int]:
        if self.short_name is None:
            raise ValueError("short_name is required to access an Earthdata product")
        kwargs: dict[str, str | int] = {
            "short_name": self.short_name,
            "count": -1,
        }
        if self.provider is not None:
            kwargs["provider"] = self.provider
        return kwargs

    def _collection_not_found_message(self) -> str:
        return f"No current NASA Earthdata collection found for {self.short_name!r}"

    def _outdated_version_message(self, current_collection) -> str:
        current_version = self._collection_umm(current_collection).get("Version")
        return (
            f"{self.short_name} version {self.requested_version} is not current. "
            f"The current NASA version is {current_version}."
        )

    def _select_collection(self, collections: list):
        current_collection = self._current_collection(collections)
        if current_collection is None:
            return None
        if self.requested_version is None:
            return current_collection

        current_version = self._collection_umm(current_collection).get("Version")
        if self.require_current_version:
            if self.requested_version != current_version:
                raise ValueError(
                    self._outdated_version_message(current_collection)
                )
            return current_collection

        requested_collections = [
            collection
            for collection in collections
            if self._collection_matches(collection)
            and self._collection_umm(collection).get("Version")
            == self.requested_version
            and self._collection_umm(collection).get("CollectionProgress")
            in self.current_collection_states
        ]
        if not requested_collections:
            return None
        return max(
            requested_collections,
            key=lambda collection: collection.get("meta", {}).get(
                "revision-date",
                "",
            ),
        )

    def resolve_product(self) -> dict[str, str | bool]:
        """Resolve and return CMR metadata for the requested product."""
        self._find_collection()
        return dict(self.product_metadata)

    def _find_collection(self) -> None:
        if self.collection is not None:
            return

        collections = self._search_datasets(**self._collection_search_kwargs())
        collection = self._select_collection(collections)
        if collection is None:
            raise FileNotFoundError(self._collection_not_found_message())

        self.collection = collection
        self.version = self._collection_umm(collection).get("Version")
        self.product_metadata = self._metadata_from_collection(collection)
        self.concept_id = self.product_metadata["concept_id"]
        if not self.concept_id:
            concept_id = collection.concept_id
            if callable(concept_id):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", FutureWarning)
                    concept_id = concept_id()
            self.concept_id = concept_id
            self.product_metadata["concept_id"] = concept_id

        self.collection_status = self.product_metadata["status"]
        self.product_title = self.product_metadata["title"]
        self.native_format = self.product_metadata["native_format"]
        self.dataset_url = self.product_metadata["dataset_url"]
        self.data_url = self.product_metadata["data_url"]

    def _harmony_download_dir(
        self,
        *,
        variables: list[str] | None,
        start: datetime,
        end: datetime,
        bbox: tuple[float, float, float, float] | None,
        dimensions: dict[str, tuple[float, float]] | None,
    ) -> Path:
        request_description = {
            "short_name": self.short_name,
            "version": self.version,
            "variables": sorted(variables) if variables is not None else None,
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
        variables: list[str] | None,
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
            variables=variables if variables is not None else ["all"],
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
                "This is not a local file-processing error.\n\n"
                "First confirm that all agreements listed here are accepted:\n"
                f"{self.eula_url}\n\n"
                "If they are already accepted, NASA's Harmony/OPeNDAP service "
                "is unavailable for this granule. Retry later, or download the "
                "complete original granule with backend='earthaccess'. The "
                "backend is not switched automatically because original "
                "granules may be large.\n\n"
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
        self._find_collection()
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
            paths = earthaccess.download(self.granules, self.download_dir)
            self.paths = [Path(path) for path in paths]
        except earthaccess.exceptions.EulaNotAccepted as error:
            raise RuntimeError(
                "NASA Earthdata login succeeded, but your account has not "
                "accepted a required data-use agreement. Open this page in a "
                "web browser, sign in, and accept the listed agreement:\n"
                f"{self.eula_url}\n"
                "Then rerun the download. atmPy cannot accept a legal "
                "agreement on your behalf."
            ) from error

    def _validate_download_request(
        self,
        start: str | date | datetime,
        end: str | date | datetime,
        *,
        variables: str | list[str] | tuple[str, ...] | None,
        bbox: tuple[float, float, float, float] | None,
        dimensions: dict[str, tuple[float, float]] | None,
        backend: DownloadBackend | None,
    ) -> tuple[
        datetime,
        datetime,
        list[str] | None,
        DownloadBackend,
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

        return (
            start_datetime,
            end_datetime,
            variable_names,
            selected_backend,
        )

    def _download_files(
        self,
        *,
        start: datetime,
        end: datetime,
        variables: list[str] | None,
        bbox: tuple[float, float, float, float] | None,
        dimensions: dict[str, tuple[float, float]] | None,
        backend: DownloadBackend,
        login: bool,
        show_progress: bool,
    ) -> list[Path]:
        self.download_dir.mkdir(parents=True, exist_ok=True)
        if login:
            self.login()

        if backend == "harmony":
            self._download_harmony(
                variables=variables,
                start=start,
                end=end,
                bbox=bbox,
                dimensions=dimensions,
                show_progress=show_progress,
            )
        else:
            self._download_granules(start=start, end=end)

        if not self.paths:
            raise RuntimeError("NASA Earthdata did not return any downloaded files")
        return list(self.paths)

    def download_files(
        self,
        start: str | date | datetime,
        end: str | date | datetime,
        *,
        variables: str | list[str] | tuple[str, ...] | None = None,
        bbox: tuple[float, float, float, float] | None = None,
        dimensions: dict[str, tuple[float, float]] | None = None,
        backend: DownloadBackend | None = None,
        login: bool = True,
        show_progress: bool = True,
    ) -> list[Path]:
        """Download Earthdata files and return their local paths.

        By default, all variables are downloaded. Harmony applies an explicit
        ``variables`` selection, ``bbox``, and ``dimensions`` before the
        download. The earthaccess backend downloads complete original granules;
        product-specific subclasses must apply requested subsets when opening
        those files.
        """
        (
            start_datetime,
            end_datetime,
            variable_names,
            selected_backend,
        ) = self._validate_download_request(
            start,
            end,
            variables=variables,
            bbox=bbox,
            dimensions=dimensions,
            backend=backend,
        )
        return self._download_files(
            start=start_datetime,
            end=end_datetime,
            variables=variable_names,
            bbox=bbox,
            dimensions=dimensions,
            backend=selected_backend,
            login=login,
            show_progress=show_progress,
        )


def create_earthdata_settings_file(path: str | Path | None = None) -> Path:
    """Create the settings file used by :class:`Earthdata`."""
    return Earthdata(settings_path=path).create_settings_file()
