"""Create immutable, row-preserving BDM cleaning products from ASCII or HDF5.

micromamba run -n cosemu python3 clean_catalogue.py --help
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time

import h5py
import numpy as np

from catalogue_core import FIELDS, RULE, audit, read_ascii, sha256


def json_dataset(handle, name, obj):
    payload = json.dumps(obj, allow_nan=False).encode()
    handle.create_dataset(name, data=np.frombuffer(payload, dtype=np.uint8),
                          compression='gzip', shuffle=True)


def load_receipt(path):
    with h5py.File(path, 'r') as f:
        return json.loads(f['receipt_json'][...].tobytes())


def read_hdf5(path):
    with h5py.File(path, 'r') as f:
        data = {key: f[key][...] for key in FIELDS}
        n = len(data['x'])
        for key, arr in data.items():
            if arr.ndim != 1 or len(arr) != n or not np.isfinite(arr).all():
                raise ValueError(f'Invalid row field {key}')
        if np.any(data['Mtot'] <= 0) or np.any(data['Mbound'] <= 0):
            raise ValueError('Nonpositive mass')
        ids = data['Nhalo']
        if np.any(ids != np.rint(ids)) or np.any(np.diff(ids) <= 0):
            raise ValueError('Nhalo must be strictly increasing integer IDs')
    return data


def write_ascii(source, output, drop, receipt, expected_hash):
    """Copy the retained lines verbatim, then independently read back and hash."""
    raw_hash, keep_hash = hashlib.sha256(), hashlib.sha256()
    count = 0
    with open(source, 'rb') as src, open(output, 'xb') as dst:
        for _ in range(8):
            line = src.readline()
            raw_hash.update(line)
            dst.write(line)
        dst.write(('# BDM-cleaning strict-v1; original row order and Nhalo retained; '
                   'receipt=' + str(receipt.resolve()) + '\n').encode())
        for line in src:
            raw_hash.update(line)
            if not line.strip() or line.lstrip().startswith(b'#'):
                raise ValueError('Unexpected non-data line in source row frame')
            if count >= len(drop):
                raise ValueError('Source has more rows than parsed')
            if not drop[count]:
                keep_hash.update(line)
                dst.write(line)
            count += 1
    if count != len(drop) or raw_hash.hexdigest() != expected_hash:
        raise ValueError('Source changed between parsing and copying')
    copied_hash = hashlib.sha256()
    copied_rows = 0
    with open(output, 'rb') as f:
        for _ in range(9):
            f.readline()
        for line in f:
            copied_hash.update(line)
            copied_rows += 1
    if copied_rows != int((~drop).sum()) or keep_hash.digest() != copied_hash.digest():
        raise ValueError('Retained rows are not byte-identical')
    return dict(retained_source_rows_sha256=keep_hash.hexdigest(),
                output_data_rows_sha256=copied_hash.hexdigest(),
                all_surviving_columns_byte_identical=True,
                output_header_lines=9, source_header_lines=8,
                output_rows=copied_rows)


def write_hdf5(source, output, drop, receipt):
    """Preserve every dataset dtype/value and attr; mask aligned row datasets."""
    n = len(drop)
    checked = []
    with h5py.File(source, 'r') as src, h5py.File(output, 'x') as dst:
        dst.attrs.update(src.attrs)
        dst.attrs['catalogue_cleaning'] = 'strict-v1'
        dst.attrs['cleaning_receipt'] = str(receipt.resolve())
        dst.attrs['original_rows'] = n
        for key, dataset in src.items():
            if not isinstance(dataset, h5py.Dataset):
                raise ValueError('Only flat catalogue HDF5 datasets are supported')
            if dataset.ndim >= 1 and dataset.shape[0] == n:
                options = {}
                if dataset.compression is not None:
                    options.update(compression=dataset.compression,
                                   compression_opts=dataset.compression_opts,
                                   shuffle=dataset.shuffle, fletcher32=dataset.fletcher32)
                out = dst.create_dataset(key, shape=(int((~drop).sum()),)+dataset.shape[1:],
                                         dtype=dataset.dtype, **options)
                out.attrs.update(dataset.attrs)
                offset = 0
                for start in range(0, n, 131072):
                    values = dataset[start:start+131072][~drop[start:start+131072]]
                    out[offset:offset+len(values)] = values
                    offset += len(values)
                checked.append(key)
            else:
                src.copy(key, dst)
    with h5py.File(source, 'r') as src, h5py.File(output, 'r') as dst:
        for key in checked:
            offset = 0
            for start in range(0, n, 131072):
                values = src[key][start:start+131072][~drop[start:start+131072]]
                if values.tobytes() != dst[key][offset:offset+len(values)].tobytes():
                    raise ValueError(f'Changed retained values in {key}')
                offset += len(values)
    return dict(all_surviving_columns_byte_identical=True,
                checked_datasets=checked, output_rows=int((~drop).sum()))


def clean(source, output, sidecar, identity, box=1024.0, velocities=False):
    t = time.monotonic()
    source, output, sidecar = [Path(p).resolve() for p in [source, output, sidecar]]
    forbidden = [Path('/cosma8/data/dp203/dc-ruan1/DESI_MGx100'),
                 Path('/cosma8/data/dp203/dc-ruan1/degrace_pilot'),
                 Path('/cosma8/data/dp203/dc-ruan1/mg_glam')]
    for p in [output, sidecar]:
        if p == source or any(p.is_relative_to(root) for root in forbidden):
            raise ValueError(f'Output path would modify a production tree: {p}')
    if output == sidecar:
        raise ValueError('Catalogue and sidecar must have distinct paths')
    if output.exists() or sidecar.exists():
        raise FileExistsError('Cleaning outputs must be new; use verified receipts to resume')
    original_stat = source.stat()
    source_hash = sha256(source)
    is_hdf = source.suffix.lower() in ['.hdf5', '.h5']
    if is_hdf:
        data = read_hdf5(source)
    else:
        data, header = read_ascii(source)
        # Actual a differs slightly from the nominal redshift label.
        identity = dict(identity, source_scale_factor=float(header[1].split()[2]))
    validation, parent, drop, edges, exact_drop = audit(data, box, velocities)
    if validation['remaining_strict_pairs'] != 0:
        raise ValueError('Strict duplicate graph still has surviving edges')
    git_root = Path(__file__).resolve().parents[4]
    git_hash = subprocess.check_output(['git', '-C', str(git_root), 'rev-parse', 'HEAD'], text=True).strip()
    tool_hashes = {p.name:sha256(p) for p in [Path(__file__), Path(__file__).with_name('catalogue_core.py')]}
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    receipt_path = sidecar
    temp = output.with_name(output.name + f'.partial-{os.getpid()}')
    if is_hdf:
        bitwise = write_hdf5(source, temp, drop, receipt_path)
        if sha256(source) != source_hash:
            raise ValueError('Source changed during HDF5 copy')
    else:
        bitwise = write_ascii(source, temp, drop, receipt_path, source_hash)
    final_stat = source.stat()
    if (original_stat.st_size, original_stat.st_mtime_ns) != (final_stat.st_size, final_stat.st_mtime_ns):
        raise ValueError('Source metadata changed during cleaning')
    n = len(drop)
    grouped_rows = np.unique(edges)
    groups = {}
    for row in grouped_rows:
        groups.setdefault(int(parent[row]), []).append(int(row))
    group_records = [dict(kept_row_index=root, kept_Nhalo=int(data['Nhalo'][root]),
                          member_row_indices=rows,
                          member_Nhalo=[int(data['Nhalo'][i]) for i in rows])
                     for root, rows in sorted(groups.items())]
    # Publish the validated data without clobbering a concurrently created path.
    os.link(temp, output)
    temp.unlink()
    receipt = dict(schema_version=1, identity=identity, rule=RULE, box_size_mpc_h=box,
        source=str(source), source_sha256=source_hash, source_rows=n,
        source_bytes=original_stat.st_size, source_mtime_ns=original_stat.st_mtime_ns,
        output=str(output), output_sha256=sha256(output),
        sidecar=str(sidecar), row_drop_mask_dataset='drop_mask',
        groups_dataset='groups_json', validation_dataset='validation_json',
        tool_git_commit=git_hash, tool_sha256=tool_hashes,
        original_row_frame='zero-based data rows; header excluded; no mass selection',
        removed_rows=int(drop.sum()), removed_fraction=float(drop.mean()),
        selected_rows=validation['selected_rows'],
        selected_removed_rows=validation['selected_removed_global_mask'],
        selected_only_strict_removed=validation['selected_strict_removed'],
        exact_removed_rows=int(exact_drop.sum()), groups_count=len(groups),
        remaining_same_bound_pairs_below_0p2=validation['remaining_same_bound_pairs_below_0p2'],
        particle_membership_validated_fraction=None,
        bitwise_validation=bitwise, elapsed_seconds=time.monotonic()-t,
        maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        job_id=os.environ.get('SLURM_JOB_ID'), array_job_id=os.environ.get('SLURM_ARRAY_JOB_ID'),
        array_task_id=os.environ.get('SLURM_ARRAY_TASK_ID'))
    # One sidecar per catalogue: masks, groups, validation and receipt share a
    # container to avoid consuming thousands of extra filesystem inodes.
    side_temp = sidecar.with_name(sidecar.name + f'.partial-{os.getpid()}')
    with h5py.File(side_temp, 'x') as f:
        for key, values in dict(drop_mask=drop, exact_drop_mask=exact_drop,
                drop_row_indices=np.flatnonzero(drop),
                dropped_Nhalo=data['Nhalo'][drop].astype(np.int64),
                dropped_representative_row_index=parent[drop]).items():
            f.create_dataset(key, data=values, compression='gzip', shuffle=True)
        f.attrs.update(source=str(source), source_nrows=n, source_sha256=source_hash,
                       tool_git_commit=git_hash, rule_json=json.dumps(RULE),
                       row_index_base=0)
        json_dataset(f, 'groups_json', group_records)
        json_dataset(f, 'validation_json', validation)
        json_dataset(f, 'receipt_json', receipt)
    os.link(side_temp, sidecar)  # Published last; presence marks a validated product.
    side_temp.unlink()
    print(json.dumps(receipt), flush=True)
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('sidecar', type=Path, help='new .cleaning.hdf5 file with masks and receipt')
    p.add_argument('--gravity', required=True)
    p.add_argument('--imodel', type=int, required=True)
    p.add_argument('--ibox', type=int, required=True)
    p.add_argument('--redshift', type=float, required=True)
    p.add_argument('--box', type=float, default=1024.0)
    p.add_argument('--velocities', action='store_true')
    args = p.parse_args()
    clean(args.source, args.output, args.sidecar,
          {k:getattr(args,k) for k in ['gravity','imodel','ibox','redshift']},
          args.box, args.velocities)


if __name__ == '__main__':
    main()
