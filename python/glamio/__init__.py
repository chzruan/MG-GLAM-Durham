"""glamio: numpy-only readers for GLAM / MG-GLAM simulation outputs."""

from .catalogue import CATALOGUE_DTYPE, COLUMN_UNITS, RECORD_BYTES, Catalogue, read_catalogue, read_catalogue_text
from .frhalo import FRHALO_DTYPE, FRHalo, read_frhalo
from .lensplanes import LensPlanes, read_lensplanes
from .meta import meta_float, read_meta
from .modes import Modes, read_modes
from .powerdm import PowerDM, read_powerdm
from .product import MetaProduct
from .schedule import SCHEDULE_DTYPE, Schedule, read_schedule
from .snapshot import PMHeader, Snapshot, read_pm_header, read_snapshot
from .subsample import Subsample, read_subsample, select_ids, splitmix64

__version__ = "0.2.0"

__all__ = ["CATALOGUE_DTYPE", "COLUMN_UNITS", "RECORD_BYTES", "Catalogue", "read_catalogue",
           "read_catalogue_text", "FRHALO_DTYPE", "FRHalo", "read_frhalo", "LensPlanes", "read_lensplanes",
           "meta_float", "read_meta", "Modes", "read_modes", "PowerDM", "read_powerdm", "MetaProduct",
           "SCHEDULE_DTYPE", "Schedule", "read_schedule", "PMHeader", "Snapshot", "read_pm_header",
           "read_snapshot", "Subsample", "read_subsample", "select_ids", "splitmix64", "__version__"]
