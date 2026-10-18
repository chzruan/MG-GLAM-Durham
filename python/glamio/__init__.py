"""glamio: numpy-only readers for GLAM / MG-GLAM simulation outputs."""

from .catalogue import CATALOGUE_DTYPE, COLUMN_UNITS, RECORD_BYTES, Catalogue, read_catalogue, read_catalogue_text
from .meta import read_meta

__version__ = "0.1.0"

__all__ = ["CATALOGUE_DTYPE", "COLUMN_UNITS", "RECORD_BYTES", "Catalogue", "read_catalogue",
           "read_catalogue_text", "read_meta", "__version__"]
