"""Bind matched IC controls to checked input files and an optional master receipt."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path
import struct


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def f32(value: float | str) -> float:
    return struct.unpack('f', struct.pack('f', float(value)))[0]


def read_receipt(path: Path) -> dict:
    return dict(tuple(part.strip() for part in line.split('=',1)) for line in path.read_text().splitlines())


def inputs(run: Path, realization: int) -> dict:
    if realization < 1:
        raise ValueError('Realization must be positive')
    lines = (run.parent/'Setup.dat').read_text().splitlines()
    scalar = lambda row: f32(lines[row].split()[0])
    seed = int((run.parent/'TableSeeds.dat').read_text().splitlines()[realization].split()[0])
    return dict(realization=realization, seed=seed, a_position=scalar(2), astep=scalar(3),
                setup_amplt=scalar(4), box_mpc_h=scalar(5), sigma8=scalar(6), hubble=scalar(7),
                omega_matter=scalar(8), omega_lambda=scalar(9), omega_baryon=scalar(10),
                nrow=int(scalar(11)), ngrid=int(scalar(12)),
                input_sha256={name:sha(run.parent/name) for name in ('Setup.dat','PkTable.dat','TableSeeds.dat')})


def configure(run: Path, master: Path | None, realization: int = 1,
              master_nrow: int = 1024, origin_ngrid: int = 2048, normalize_only: bool = False) -> dict:
    current = inputs(run,realization)
    if current['nrow'] > master_nrow or master_nrow % current['nrow']:
        raise ValueError('NROW must divide the chosen master NROW')
    alpha = -1.
    master_evidence = None
    if master is None:
        if current['nrow'] != master_nrow:
            raise ValueError('A lower-resolution IC requires the completed master receipt')
    else:
        source = inputs(master, realization)
        recorded = json.loads((master/'matched_ic_inputs.json').read_text())
        if recorded['inputs'] != source:
            raise ValueError('Master input files changed after its controls were written')
        if recorded['controls_sha256'] != sha(master/'matched_ic.nml'):
            raise ValueError('Master controls changed after they were recorded')
        actual = read_receipt(master/'matched_ic_receipt.txt')
        for key in ('realization','seed','nrow','ngrid'):
            if int(actual[key]) != source[key]:
                raise ValueError(f'Master {key} differs from configured inputs')
        for key in ('a_position','astep','box_mpc_h','setup_amplt'):
            if float(actual[key]) != source[key]:
                raise ValueError(f'Master {key} differs from configured inputs')
        if int(actual['master_nrow']) != master_nrow or int(actual['origin_ngrid']) != origin_ngrid:
            raise ValueError('Master Fourier stride or physical origin differs')
        if source['nrow'] != master_nrow:
            raise ValueError('Master receipt is from a lower-resolution IC')
        for key in ('realization','seed','a_position','box_mpc_h','sigma8','hubble',
                    'omega_matter','omega_lambda','omega_baryon'):
            if source[key] != current[key]:
                raise ValueError(f'Matched inputs disagree: {key}')
        for name in ('PkTable.dat','TableSeeds.dat'):
            if source['input_sha256'][name] != current['input_sha256'][name]:
                raise ValueError(f'Matched input bytes disagree: {name}')
        alpha = float(actual['alpha'])
        if not (alpha > 0.) or f32(alpha) != alpha:
            raise ValueError('Master alpha must round-trip exactly through float32')
        master_evidence = dict(run_dir=str(master.resolve()), receipt_sha256=sha(master/'matched_ic_receipt.txt'),
                               controls_sha256=sha(master/'matched_ic.nml'),
                               inputs_sha256=sha(master/'matched_ic_inputs.json'))
    controls = ('&matched_ic\n'
                f' ic_master_nrow={master_nrow}, ic_origin_ngrid={origin_ngrid},\n'
                f' ic_alpha={alpha:.17g}, ic_normalize_only={".true." if normalize_only else ".false."}\n/\n')
    result = dict(format='native-luxury-master-stride-v1', inputs=current,
                  master_nrow=master_nrow, origin_ngrid=origin_ngrid,
                  origin_mpc_h=current['box_mpc_h']/(2*origin_ngrid), alpha_requested=alpha,
                  normalize_only=normalize_only, master=master_evidence,
                  controls_sha256=hashlib.sha256(controls.encode()).hexdigest())
    run.mkdir(parents=True, exist_ok=True)
    outputs = {run/'matched_ic.nml': controls.encode(),
               run/'matched_ic_inputs.json': (json.dumps(result,indent=2)+'\n').encode()}
    # Check both destinations before publishing either; never replace provenance.
    for path, content in outputs.items():
        if path.exists() and path.read_bytes() != content:
            raise FileExistsError(f'Existing matched IC controls differ: {path}')
    for path, content in outputs.items():
        if path.exists():
            continue
        with tempfile.NamedTemporaryFile(prefix='.matched-ic-control-',dir=run,delete=False) as out:
            temporary = Path(out.name)
            out.write(content)
            out.flush()
            os.fsync(out.fileno())
        try:
            try:
                os.link(temporary,path)  # Atomic publication without clobbering a concurrent writer.
            except FileExistsError:
                if path.read_bytes() != content:
                    raise FileExistsError(f'Concurrent matched IC controls differ: {path}')
        finally:
            temporary.unlink()
    return result


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--master-run-dir',type=Path)
    parser.add_argument('--realization',type=int,default=1)
    parser.add_argument('--master-nrow',type=int,default=1024)
    parser.add_argument('--origin-ngrid',type=int,default=2048)
    parser.add_argument('--normalize-only',action='store_true')
    args=parser.parse_args()
    print(json.dumps(configure(args.run_dir,args.master_run_dir,args.realization,
                               args.master_nrow,args.origin_ngrid,args.normalize_only),indent=2))


if __name__ == '__main__':
    main()
