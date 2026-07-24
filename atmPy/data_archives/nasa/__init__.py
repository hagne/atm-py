"""Access to data products distributed by NASA."""

from .merra2 import (
    Merra2,
    create_earthdata_settings_file,
    download_total_column_ozone,
)

__all__ = [
    "Merra2",
    "create_earthdata_settings_file",
    "download_total_column_ozone",
]
