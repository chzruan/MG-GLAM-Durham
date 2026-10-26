"""Reader of the in-line Rockstar catalogue CATALOGS/Rockstar.<step>.<box>.bin + .meta (test harness;
glamio gets the production reader). The .meta lists every column as
``column_NN = name type offset unit...``; the records are raw little-endian, ``record_bytes`` each.

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
    names, formats, offsets = [], [], []
    for i in range(1, int(meta["columns"]) + 1):
        name, typ, off = meta[f"column_{i:02d}"].split()[:3]
        names.append(name)
        formats.append(typ)
        offsets.append(int(off))
    return np.dtype({"names": names, "formats": formats, "offsets": offsets,
                     "itemsize": int(meta["record_bytes"])})


def read_rsbin(base, allow_status=False):
    """base: path without .bin/.meta. Returns (records, meta). A catalogue whose status is not ``ok`` (skipped
    or failed: zero rows) raises ValueError unless allow_status."""
    meta = read_meta(base + ".meta")
    if meta.get("status") != "ok" and not allow_status:
        raise ValueError(f"{base}: status {meta.get('status')} ({meta.get('reason')})")
    dt = dtype_from_meta(meta)
    rows = int(meta["rows"])
    data = np.fromfile(base + ".bin", dtype=dt)
    if len(data) != rows:
        raise ValueError(f"{base}.bin holds {len(data)} records, the .meta says {rows}")
    return data, meta
