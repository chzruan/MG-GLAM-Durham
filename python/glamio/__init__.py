"""glamio: numpy-only readers for GLAM / MG-GLAM simulation outputs."""

from .catalogue import CATALOGUE_DTYPE, COLUMN_UNITS, RECORD_BYTES, Catalogue, read_catalogue, read_catalogue_text
from .frhalo import FRHALO_DTYPE, FRHalo, read_frhalo
from .lensplanes import AXES, LensPlanes, read_lensplanes
from .meta import meta_float, read_meta
from .modes import MODE_PRODUCTS, Modes, read_modes
from .powerdm import PowerDM, read_powerdm
from .product import MetaProduct, check_size, load_array, load_meta, meta_int, resolve_pair
from .run import SIDECAR_PRODUCTS, product_path, products_at
from .schedule import SCHEDULE_DTYPE, Schedule, read_schedule
from .snapshot import MAX_PARTICLES_DEFAULT, PAGE, RECORDS_PER_FILE, RHO_CRIT_GLAM, PMHeader, Snapshot, read_pm_header, read_snapshot
from .subsample import MAX_IDS_DEFAULT, NQ, SALT_DEFAULT, Subsample, keep_mask, read_subsample, select_ids, splitmix64

__version__ = "0.2.1"

__all__ = ["CATALOGUE_DTYPE", "COLUMN_UNITS", "RECORD_BYTES", "Catalogue", "read_catalogue",
           "read_catalogue_text", "FRHALO_DTYPE", "FRHalo", "read_frhalo", "AXES", "LensPlanes", "read_lensplanes",
           "meta_float", "read_meta", "MODE_PRODUCTS", "Modes", "read_modes", "PowerDM", "read_powerdm",
           "MetaProduct", "check_size", "load_array", "load_meta", "meta_int", "resolve_pair",
           "SIDECAR_PRODUCTS", "product_path", "products_at",
           "SCHEDULE_DTYPE", "Schedule", "read_schedule",
           "MAX_PARTICLES_DEFAULT", "PAGE", "RECORDS_PER_FILE", "RHO_CRIT_GLAM", "PMHeader", "Snapshot",
           "read_pm_header", "read_snapshot",
           "MAX_IDS_DEFAULT", "NQ", "SALT_DEFAULT", "Subsample", "keep_mask", "read_subsample", "select_ids",
           "splitmix64", "__version__"]
