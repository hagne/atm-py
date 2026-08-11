"""Download current MERRA-2-related data from NASA Earthdata.

The catalog is not hard-coded. :meth:`Merra2.available_products` queries NASA's
Common Metadata Repository (CMR) and returns the current version, title, native
format, and exact GES DISC dataset URL for every collection that NASA identifies
with a MERRA-2 project or platform tag. This includes the standard reanalysis,
climate-statistics, gridded-observation (GIO), and M2-SCREAM collections.

NASA documentation and catalogs:

* MERRA-2 project and data access:
  https://disc.gsfc.nasa.gov/information/mission-project?title=MERRA-2
* GMAO MERRA-2 documentation and file specifications:
  https://gmao.gsfc.nasa.gov/gmao-products/merra-2/documentation_merra-2/
* GES DISC collection catalog:
  https://cmr.earthdata.nasa.gov/search/site/collections/directory/GES_DISC/gov.nasa.eosdis
* Most standard collection pages use:
  https://disc.gsfc.nasa.gov/datacollection/{SHORT_NAME}_{VERSION}.html

For example, the single-level hourly collection is documented at
https://disc.gsfc.nasa.gov/datacollection/M2T1NXSLV_5.12.4.html, the current
monthly extremes collection at
https://disc.gsfc.nasa.gov/datacollection/M2SMNXEDI_2.html, a GIO collection at
https://disc.gsfc.nasa.gov/datasets/M2_MHS_METOP-B_1/summary, and M2-SCREAM at
https://disc.gsfc.nasa.gov/datasets/GMAO_M2SCREAM_INST3_CHEM_1/summary.
Use ``Merra2.available_products().dataset_url`` for the corresponding URL of
each dataset in the live catalog.
"""

import configparser
import hashlib
import json
import os
import re
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
    """Access a current MERRA-2-related collection using NASA Earthdata.

    Parameters
    ----------
    short_name
        NASA collection short name, for example ``"M2T1NXSLV"`` for hourly
        single-level diagnostics, ``"M2T1NXAER"`` for hourly aerosol
        diagnostics, or ``"M2I3NPASM"`` for three-hourly pressure-level
        assimilated meteorology.
    version
        Collection version. The default, ``None``, resolves the current version
        from NASA CMR. An explicitly requested older version is rejected so an
        outdated collection cannot be downloaded accidentally.
    backend
        ``"harmony"`` requests variable and spatial subsetting before download.
        Not every MERRA-2-related collection offers a Harmony service; use
        ``"earthaccess"`` to download its original files.

    Notes
    -----
    Dataset titles, versions, formats, and exact landing-page URLs are available
    from :meth:`available_products`. The authoritative general documentation is
    https://gmao.gsfc.nasa.gov/gmao-products/merra-2/documentation_merra-2/.
    Each product's NASA website is also exposed as ``dataset_url`` after
    :meth:`resolve_product` or :meth:`download` is called.
    """

    short_name = "M2T1NXSLV"
    version = None
    ozone_variable = "TO3"
    ozone_units = "Dobsons"
    collection_provider = "GES_DISC"
    current_collection_states = frozenset({"ACTIVE", "COMPLETE"})
    merra2_projects = frozenset({"MERRA-2", "MERRA-2 Observation"})
    merra2_platform = "MERRA-2"
    project_url = (
        "https://disc.gsfc.nasa.gov/information/mission-project?title=MERRA-2"
    )
    documentation_url = (
        "https://gmao.gsfc.nasa.gov/gmao-products/merra-2/"
        "documentation_merra-2/"
    )
    collection_catalog_url = (
        "https://cmr.earthdata.nasa.gov/search/site/collections/"
        "directory/GES_DISC/gov.nasa.eosdis"
    )
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
        version: str | None = None,
        backend: DownloadBackend = "harmony",
    ):
        if backend not in ("harmony", "earthaccess"):
            raise ValueError("backend must be 'harmony' or 'earthaccess'")
        if not isinstance(short_name, str) or not short_name.strip():
            raise ValueError("short_name must be a NASA collection short name")
        if version is not None and (
            not isinstance(version, str) or not version.strip()
        ):
            raise ValueError("version must be a non-empty string or None")

        self.short_name = short_name.strip()
        self.requested_version = version.strip() if version is not None else None
        self.version = self.requested_version
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

    @staticmethod
    def _collection_umm(collection) -> dict:
        return collection.get("umm", {})

    @classmethod
    def _is_merra2_collection(cls, collection) -> bool:
        umm = cls._collection_umm(collection)
        projects = {
            project.get("ShortName")
            for project in umm.get("Projects", [])
        }
        platforms = {
            platform.get("ShortName")
            for platform in umm.get("Platforms", [])
        }
        return bool(
            projects.intersection(cls.merra2_projects)
            or cls.merra2_platform in platforms
        )

    @classmethod
    def _version_key(cls, collection) -> tuple[int, ...]:
        version = cls._collection_umm(collection).get("Version", "")
        return tuple(int(part) for part in re.findall(r"\d+", version))

    @classmethod
    def _current_collection(cls, collections: list):
        candidates = [
            collection
            for collection in collections
            if cls._is_merra2_collection(collection)
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
        related_urls = cls._collection_umm(collection).get(
            "RelatedUrls",
            [],
        )
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
        short_name = umm.get("ShortName", "")
        version = umm.get("Version", "")
        dataset_url = cls._collection_url(
            collection,
            "DATA SET LANDING PAGE",
        )
        if dataset_url is None:
            dataset_url = (
                "https://disc.gsfc.nasa.gov/datacollection/"
                f"{short_name}_{version}.html"
            )

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
            "short_name": short_name,
            "version": version,
            "status": umm.get("CollectionProgress", ""),
            "title": umm.get("EntryTitle", ""),
            "native_format": ", ".join(sorted(formats)),
            "concept_id": meta.get("concept-id", ""),
            "dataset_url": dataset_url,
            "data_url": cls._collection_url(collection, "GET DATA") or "",
            "associated_services": bool(associations.get("services")),
        }

    @classmethod
    def available_products(cls) -> xr.Dataset:
        """Return NASA's current MERRA-2-related product catalog.

        The returned dataset has one ``short_name`` coordinate per current
        collection and includes ``version``, ``title``, ``native_format``,
        ``status``, ``concept_id``, ``dataset_url``, ``data_url``, and
        ``associated_services``. ``dataset_url`` contains the exact NASA GES
        DISC website for each product.

        CMR uses several metadata routes for the MERRA-2 product families.
        This method combines the standard MERRA-2 project, the MERRA-2 platform,
        and the MERRA-2 Observation project searches, then retains only the
        newest ``ACTIVE`` or ``COMPLETE`` version of each short name.

        See also
        --------
        NASA MERRA-2 documentation:
        https://gmao.gsfc.nasa.gov/gmao-products/merra-2/documentation_merra-2/
        NASA GES DISC catalog:
        https://cmr.earthdata.nasa.gov/search/site/collections/directory/GES_DISC/gov.nasa.eosdis
        """
        searches = (
            {"project": "MERRA-2"},
            {"platform": cls.merra2_platform},
            {"project": "MERRA-2 Observation"},
        )
        collections_by_concept_id = {}
        for search in searches:
            collections = earthaccess.search_datasets(
                provider=cls.collection_provider,
                count=-1,
                **search,
            )
            for collection in collections:
                concept_id = collection.get("meta", {}).get("concept-id")
                if concept_id:
                    collections_by_concept_id[concept_id] = collection

        collections_by_short_name: dict[str, list] = {}
        for collection in collections_by_concept_id.values():
            if not cls._is_merra2_collection(collection):
                continue
            short_name = cls._collection_umm(collection).get("ShortName")
            if short_name:
                collections_by_short_name.setdefault(short_name, []).append(
                    collection
                )

        current_collections = [
            cls._current_collection(collection_versions)
            for collection_versions in collections_by_short_name.values()
        ]
        metadata = sorted(
            (
                cls._metadata_from_collection(collection)
                for collection in current_collections
                if collection is not None
            ),
            key=lambda product: product["short_name"],
        )
        if not metadata:
            raise FileNotFoundError(
                "NASA CMR did not return any current MERRA-2 collections"
            )

        fields = (
            "version",
            "status",
            "title",
            "native_format",
            "concept_id",
            "dataset_url",
            "data_url",
            "associated_services",
        )
        return xr.Dataset(
            data_vars={
                field: (
                    "short_name",
                    [product[field] for product in metadata],
                )
                for field in fields
            },
            coords={
                "short_name": [
                    product["short_name"] for product in metadata
                ]
            },
            attrs={
                "source": "NASA Common Metadata Repository",
                "project_url": cls.project_url,
                "documentation_url": cls.documentation_url,
                "collection_catalog_url": cls.collection_catalog_url,
                "selection": (
                    "Newest ACTIVE or COMPLETE version for each short name"
                ),
            },
        )

    def resolve_product(self) -> dict[str, str | bool]:
        """Resolve and return metadata for the current requested product.

        ``dataset_url`` in the returned dictionary is the product's exact NASA
        GES DISC website. If ``version`` was omitted during construction, the
        resolved current version is stored on :attr:`version`.
        """
        self._find_collection()
        return dict(self.product_metadata)

    def _find_collection(self) -> None:
        if self.collection is not None:
            return

        collections = earthaccess.search_datasets(
            short_name=self.short_name,
            provider=self.collection_provider,
            count=-1,
        )
        current_collection = self._current_collection(collections)
        if current_collection is None:
            raise FileNotFoundError(
                f"No current MERRA-2-related collection found for "
                f"{self.short_name!r}. Use Merra2.available_products() to "
                "inspect valid short names, versions, descriptions, and NASA "
                "dataset URLs."
            )

        current_version = self._collection_umm(current_collection).get("Version")
        if (
            self.requested_version is not None
            and self.requested_version != current_version
        ):
            metadata = self._metadata_from_collection(current_collection)
            raise ValueError(
                f"{self.short_name} version {self.requested_version} is not "
                f"current. The current NASA version is {current_version}: "
                f"{metadata['dataset_url']}"
            )

        self.collection = current_collection
        self.version = current_version
        self.product_metadata = self._metadata_from_collection(self.collection)
        self.concept_id = self.product_metadata["concept_id"]
        if not self.concept_id:
            concept_id = self.collection.concept_id
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

        ``short_name`` and the current product version are configured on the
        :class:`Merra2` instance. Call :meth:`available_products` for current
        short names, product descriptions, and NASA dataset URLs.
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
        product_attributes = {
            "merra2_short_name": self.short_name,
            "merra2_version": self.version,
            "merra2_product_title": self.product_title,
            "merra2_dataset_url": self.dataset_url,
        }
        self.dataset.attrs.update(
            {
                name: value
                for name, value in product_attributes.items()
                if value is not None
            }
        )
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
