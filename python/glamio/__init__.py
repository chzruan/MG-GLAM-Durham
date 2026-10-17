"""glamio: numpy-only readers for GLAM / MG-GLAM simulation outputs."""

from .catalogue import CATALOGUE_DTYPE, Catalogue, read_catalogue
from .meta import read_meta

__version__ = "0.1.0"

__all__ = ["CATALOGUE_DTYPE", "Catalogue", "read_catalogue", "read_meta", "__version__"]
