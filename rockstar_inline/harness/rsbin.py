"""Reader of the in-line Rockstar catalogue CATALOGS/Rockstar.<step>.<box>.bin + .meta (test harness;
glamio has the production reader). The .meta lists the columns as space-separated ``columns``, ``dtypes``
(little-endian) and ``offsets``; the records are raw, ``record_bytes`` each.

Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later."""
import numpy as np


def read_meta(path):
    meta = {}
    with open(path) as f:
        for line in f:
            if " = " in line:
                k, v = line.rstrip("\n").split(" = ", 1)
                meta[k.strip()] = v.strip()
    return meta


def dtype_from_meta(meta):
    names, types, offsets = meta["columns"].split(), meta["dtypes"].split(), meta["offsets"].split()
    return np.dtype({"names": names, "formats": ["<" + t for t in types], "offsets": [int(o) for o in offsets],
                     "itemsize": int(meta["record_bytes"])})


def read_rsbin(base, allow_status=False):
    """base: path without .bin/.meta. Returns (records, meta). A catalogue whose status is not ``ok`` (skipped
    or failed: zero rows) raises ValueError unless allow_status."""
    meta = read_meta(base + ".meta")
    if meta.get("status", "").split()[:1] != ["ok"] and not allow_status:
        raise ValueError(f"{base}: status {meta.get('status')} ({meta.get('reason')})")
    dt = dtype_from_meta(meta)
    rows = int(meta["rows"])
    data = np.fromfile(base + ".bin", dtype=dt)
    if len(data) != rows:
        raise ValueError(f"{base}.bin holds {len(data)} records, the .meta says {rows}")
    return data, meta
