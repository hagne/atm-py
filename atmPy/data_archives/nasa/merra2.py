"""Download and interpret current MERRA-2 data from NASA Earthdata.

The product catalog is discovered from NASA's Common Metadata Repository (CMR)
rather than hard-coded. Earthdata authentication, collection resolution, and
file retrieval are implemented by :class:`~atmPy.data_archives.nasa.Earthdata`;
this module contains the MERRA-2 collection and dataset conventions.
"""

from datetime import date, datetime
from pathlib import Path

import xarray as xr

from .earthdata import (
    AuthStrategy,
    DownloadBackend,
    Earthdata,
    create_earthdata_settings_file,
)


class Merra2(Earthdata):
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
    collection_provider = "GES_DISC"
    require_current_version = True
    cache_namespace = "merra2"
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
        super().__init__(
            settings_path=settings_path,
            download_dir=download_dir,
            auth_strategy=auth_strategy,
            short_name=short_name,
            provider=self.collection_provider,
            version=version,
            backend=backend,
        )
        self.chunks = chunks
        self.dataset: xr.Dataset | None = None

    @classmethod
    def _is_merra2_collection(cls, collection) -> bool:
        umm = cls._collection_umm(collection)
        projects = {
            project.get("ShortName") for project in umm.get("Projects", [])
        }
        platforms = {
            platform.get("ShortName") for platform in umm.get("Platforms", [])
        }
        return bool(
            projects.intersection(cls.merra2_projects)
            or cls.merra2_platform in platforms
        )

    @classmethod
    def _collection_matches(cls, collection) -> bool:
        return cls._is_merra2_collection(collection)

    @classmethod
    def _metadata_from_collection(
        cls,
        collection,
    ) -> dict[str, str | bool]:
        metadata = super()._metadata_from_collection(collection)
        if not metadata["dataset_url"]:
            metadata["dataset_url"] = (
                "https://disc.gsfc.nasa.gov/datacollection/"
                f"{metadata['short_name']}_{metadata['version']}.html"
            )
        return metadata

    @classmethod
    def available_products(cls) -> xr.Dataset:
        """Return NASA's current MERRA-2-related product catalog.

        The returned dataset has one ``short_name`` coordinate per current
        collection and includes ``version``, ``title``, ``native_format``,
        ``status``, ``concept_id``, ``dataset_url``, ``data_url``, and
        ``associated_services``. CMR's MERRA-2 project, platform, and
        observation-project searches are combined before selecting the newest
        ``ACTIVE`` or ``COMPLETE`` version of each short name.
        """
        searches = (
            {"project": "MERRA-2"},
            {"platform": cls.merra2_platform},
            {"project": "MERRA-2 Observation"},
        )
        collections_by_concept_id = {}
        for search in searches:
            collections = cls._search_datasets(
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
                "short_name": [product["short_name"] for product in metadata]
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

    def _collection_not_found_message(self) -> str:
        return (
            f"No current MERRA-2-related collection found for "
            f"{self.short_name!r}. Use Merra2.available_products() to inspect "
            "valid short names, versions, descriptions, and NASA dataset URLs."
        )

    def _outdated_version_message(self, current_collection) -> str:
        metadata = self._metadata_from_collection(current_collection)
        return (
            f"{self.short_name} version {self.requested_version} is not current. "
            f"The current NASA version is {metadata['version']}: "
            f"{metadata['dataset_url']}"
        )

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
        variables: str | list[str] | tuple[str, ...] | None = None,
        bbox: tuple[float, float, float, float] | None = None,
        dimensions: dict[str, tuple[float, float]] | None = None,
        backend: DownloadBackend | None = None,
        login: bool = True,
        show_progress: bool = True,
    ) -> xr.Dataset:
        """Download variables from a MERRA-2 collection into an xarray dataset.

        By default, all variables are returned. ``bbox`` must be ordered as
        west, south, east, north. ``dimensions`` maps dimension names to
        inclusive minimum and maximum values. Harmony performs these subsets
        before download. The ``earthaccess`` backend downloads complete
        granules and applies the subsets locally.
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
        paths = self._download_files(
            start=start_datetime,
            end=end_datetime,
            variables=variable_names,
            bbox=bbox,
            dimensions=dimensions,
            backend=selected_backend,
            login=login,
            show_progress=show_progress,
        )

        ds = xr.open_mfdataset(
            paths,
            combine="by_coords",
            chunks=self.chunks,
        )
        if variable_names is not None:
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
