"""Small shared helpers; use micromamba run -n cosemu python3 -B."""
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(os.environ.get('BDM_CONVERGENCE_ROOT', Path(__file__).resolve().parent)).resolve()
REPO = ROOT.parents[2]
WORK = ROOT / 'work'
BIN = WORK / 'bin'
BASE = 'e289c3f7ea28d38390a1663e50bb2784b85abc84'
FINDER_SHA = '39f7e514d1db9fb3b9c7545fc7f3269b4e20b5ee2d1366945f8916b997e4cd5d'
MATRIX = {'A': (256, 2048), 'B': (512, 1024), 'C': (512, 2048),
          'D': (512, 4096), 'E': (1024, 2048), 'F': (1024, 4096),
          'T': (1024, 4096)}
CONFIG = ('! Common requested configuration for legacy and v3\n'
          'iVirial = 1 ! virial overdensity\n'
          'MassMin = 2.5e12 ! bound mass in Msun/h\n'
          'Rext = 0.15 ! retained aperture correction\n')


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while chunk := stream.read(8 * 1024**2):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def native_env(threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE',
               OMP_PROC_BIND='close', OMP_PLACES='cores', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', LD_LIBRARY_PATH=os.environ['BDM_AUDIT_NATIVE_LIBS'])
    for key in ['LIBRARY_PATH', 'MAKEFLAGS', 'MFLAGS', 'MAKEOVERRIDES']:
        env.pop(key, None)
    return env


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def checker():
    module = load_module('bdm_prior_checker', REPO / 'BDM-refine/validation/20260906-n1024/run_validation.py')
    module.ROOT = REPO
    return module


def git(*args):
    return subprocess.check_output(['git', *args], cwd=REPO, text=True).strip()


def file_manifest(paths):
    return {str(Path(p).resolve()): {'sha256': sha(p), 'bytes': Path(p).stat().st_size} for p in paths}


def verify_manifest(manifest):
    for name, expected in manifest.items():
        path = Path(name)
        if not path.is_file() or path.stat().st_size != expected['bytes'] or sha(path) != expected['sha256']:
            raise ValueError(f'Changed or missing frozen file: {name}')
