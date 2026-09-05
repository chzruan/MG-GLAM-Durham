"""Read a sidecar's original-row mask without changing raw files or mirrors."""
from pathlib import Path
import json

import h5py
import numpy as np

from catalogue_core import sha256


def load_cleaning(sidecar, source=None, verify_source=True):
    """Return (drop_mask, receipt); True means drop that original data row.

    Apply the full-length mask before making mass cuts. Nhalo is an output ID,
    not a zero-based row index. For HDF5 mirrors, first establish their unchanged
    alignment with the ASCII source; a matching row count alone does not prove it.
    """
    with h5py.File(sidecar,'r') as f:
        receipt=json.loads(f['receipt_json'][...].tobytes())
        drop=f['drop_mask'][...]
    if drop.shape!=(receipt['source_rows'],) or drop.dtype!=np.bool_:
        raise ValueError('Invalid source row frame')
    if int(drop.sum())!=receipt['removed_rows']:
        raise ValueError('Mask/receipt count mismatch')
    if verify_source:
        path=Path(receipt['source'] if source is None else source)
        if sha256(path)!=receipt['source_sha256']:
            raise ValueError('Source SHA256 does not match the cleaning frame')
    return drop,receipt
