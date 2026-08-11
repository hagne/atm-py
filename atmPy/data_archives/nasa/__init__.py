"""Access to data products distributed by NASA."""

from .earthdata import Earthdata, create_earthdata_settings_file
from .merra2 import (
    Merra2,
)

__all__ = [
    "Earthdata",
    "Merra2",
    "create_earthdata_settings_file",
]
