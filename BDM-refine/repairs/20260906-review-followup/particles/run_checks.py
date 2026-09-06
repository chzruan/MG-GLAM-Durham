"""Bounded N3/N4 controls; always run via micromamba run -n cosemu python3 -B."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
BASE = "54f53aae2aac6700d1ddb8490bebb0719f71165d"
for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"


def digest(value):
    return hashlib.sha256(value).hexdigest()


def extract(source, name, kind="subroutine"):
    pattern = rf"^\s*(?:pure\s+)?(?:real\*8\s+)?{kind}\s+{name}\b.*?^\s*end\s+{kind}\s+{name}\b[^\n]*"
    found = re.search(pattern, source, re.I | re.M | re.S)
    if not found:
        raise ValueError(name)
    return found.group()


def build(work, source, compiler, checked):
    structures = re.search(r"^module\s+Structures\b.*?^end module\s+Structures", source,
                           re.I | re.M | re.S).group()
    stub = """module Tools
real :: Box=32.,AEXPN=.8
integer :: NGRID=128,NROW=32
integer*8 :: Nparticles=0,memoryWords=0,peakWords=0
real,allocatable :: Xpar(:),Ypar(:),Zpar(:),VX(:),VY(:),VZ(:)
contains
real function seconds()
use omp_lib, only: omp_get_wtime
seconds=real(omp_get_wtime())
end function
real function Memory(n)
integer*8 :: n
memoryWords=memoryWords+n
peakWords=max(peakWords,memoryWords)
Memory=real(dble(memoryWords)*4.d0/1024.d0**3)
end function
end module
module LinkerList
use Structures
use Tools
contains
"""
    names = ["List", "Limits", "RescaleCoords", "AddBuffer", "RemoveBuffer", "PrepareParticleSearch",
             "SizeList", "BdmParticlePosition"]
    generated = structures + "\n" + stub + "\n".join(extract(source, name) for name in names)
    generated += "\n" + extract(source, "BdmParticleCoordinate", "function")
    generated += "\nend module\n" + (HERE / "cases.f90").read_text()
    (work / "cases.f90").write_text(generated)
    if compiler == "gfortran":
        flags = ["-O0", "-g", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", "-ftrapv"] if checked else ["-O3"]
        command = [compiler, *flags, "-fopenmp", "-ffree-line-length-none", "cases.f90", "-o", "cases"]
    else:
        flags = ["-O0", "-g", "-check", "bounds", "-fpe0"] if checked else ["-O3"]
        command = [compiler, *flags, "-fp-model", "precise", "-qopenmp", "-extend-source", "cases.f90", "-o", "cases"]
    env = dict(os.environ)
    if compiler == "ifx":
        env.pop("LIBRARY_PATH", None)
        env["LD_LIBRARY_PATH"] = os.environ["BDM_AUDIT_NATIVE_LIBS"]
    result = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True, timeout=60)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    return {"command": command, "stderr": result.stderr, "generated_sha256": digest(generated.encode())}


def run(work, case, n, threads, compiler, repetitions=1, expected_error=None):
    env = {**os.environ, "OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE",
           "OMP_PROC_BIND": "false", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if compiler == "ifx":
        env["LD_LIBRARY_PATH"] = os.environ["BDM_AUDIT_NATIVE_LIBS"]
    started = time.monotonic()
    result = subprocess.run([str(work / "cases"), case, str(n), str(repetitions)], cwd=work,
                            env=env, text=True, capture_output=True, timeout=40)
    record = {"case": case, "n": n, "threads": threads, "returncode": result.returncode,
              "elapsed_seconds": time.monotonic()-started}
    if expected_error:
        errors = [expected_error] if isinstance(expected_error, str) else expected_error
        output = (result.stdout+result.stderr).lower()
        if not result.returncode or "missing" in output or not any(error.lower() in output for error in errors):
            raise RuntimeError(json.dumps(record)+"\n"+result.stdout+result.stderr)
        record["expected_error"] = expected_error
    else:
        if result.returncode or "FOLLOWUP_PASS" not in result.stdout:
            raise RuntimeError(json.dumps(record)+"\n"+result.stdout+result.stderr)
        record["output_sha256"] = digest((work / "result.bin").read_bytes())
        match = re.search(r"FOLLOWUP_METRICS\s+(\w+)\s+([\d.Ee+-]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", result.stdout)
        if match:
            record.update(kernel=match[1], seconds_per_call=float(match[2]), repetitions=int(match[3]),
                          original_particles=int(match[4]), buffered_particles=int(match[5]), measured_memory_words=int(match[6]))
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler", choices=["gfortran", "ifx"], default="gfortran")
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--particles", type=int, default=262144)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--repetitions", type=int, default=6)
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = (REPO / "PMP2linker.f90").read_text()
    baseline = subprocess.check_output(["git", "show", f"{BASE}:PMP2linker.f90"], cwd=REPO, text=True)
    compilation, records, hashes = {}, [], {}
    cases = [("buffer", n) for n in [0, 1, 16383, 16384, 16385, 32769]]
    cases += [(case, 32769) for case in ["buffer_fractional", "buffer_halfbox", "list_periodic",
        "list_i16_range", "list_i16_low", "list_i16_high", "list_i32_low", "list_i32_high", "list_empty_slabs"]]
    extended = [("list_i32_min", 10001), ("list_i32_max", 10001)]
    errors = {"invalid_negative": "[0,Box)", "invalid_upper": "[0,Box)",
              "invalid_nan": ["original analysis coordinates", "SIGFPE", "floating invalid"],
              "invalid_count": "count mismatch", "invalid_active": "already active", "invalid_capacity": "safe int64",
              "invalid_negative_count": "safe int64", "invalid_cell": "Invalid BDM particle list geometry",
              "invalid_bounds": "Invalid BDM particle list bounds", "invalid_list_nan": "Nonfinite BDM particle cell coordinate",
              "invalid_cache_memory": "BDM z-cell cache exceeds"}
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="bdm-review-particles-") as temporary:
        scratch = Path(temporary)
        for checked in ([False] if args.benchmark else [True, False]):
            mode = "checked" if checked else "optimized"
            builds = {}
            for variant, text in [("baseline", baseline), ("updated", source)]:
                work = scratch / (mode+"-"+variant); work.mkdir(); builds[variant] = work
                compilation[mode+"-"+variant] = build(work, text, args.compiler, checked)
            for threads in args.threads:
                if args.benchmark:
                    for case in ["list_benchmark", "buffer_benchmark"]:
                        for sample in range(args.samples):
                            for variant in (["baseline", "updated"] if sample % 2 == 0 else ["updated", "baseline"]):
                                record = run(builds[variant], case, args.particles, threads, args.compiler, args.repetitions)
                                record.update(mode=mode, variant=variant, sample=sample); records.append(record)
                                key = (case, args.particles)
                                if key in hashes:
                                    assert record["output_sha256"] == hashes[key], record
                                hashes[key] = record["output_sha256"]
                else:
                    for case, n in cases:
                        for variant in ["baseline", "updated"]:
                            record = run(builds[variant], case, n, threads, args.compiler)
                            record.update(mode=mode, variant=variant); records.append(record)
                            key = (case, n)
                            if key in hashes:
                                assert record["output_sha256"] == hashes[key], record
                            hashes[key] = record["output_sha256"]
                    for case, n in extended:
                        record = run(builds["updated"], case, n, threads, args.compiler)
                        record.update(mode=mode, variant="updated", baseline_skipped="Original default-integer loop/conversion is undefined at int32 extrema")
                        records.append(record)
                    for case, error in errors.items():
                        if case == "invalid_cache_memory" and threads == 1:
                            continue
                        record = run(builds["updated"], case, 10001, threads, args.compiler, expected_error=error)
                        record.update(mode=mode, variant="updated"); records.append(record)
    medians = {}
    if args.benchmark:
        for kernel in ["list", "buffer"]:
            medians[kernel] = {}
            for threads in args.threads:
                values = {variant: statistics.median(r["seconds_per_call"] for r in records
                    if r["kernel"] == kernel and r["threads"] == threads and r["variant"] == variant)
                    for variant in ["baseline", "updated"]}
                values["baseline_over_updated"] = values["baseline"] / values["updated"]
                medians[kernel][str(threads)] = values
    report = {"baseline_commit": BASE, "compiler": args.compiler,
        "compiler_version": subprocess.check_output([args.compiler, "--version"], text=True).splitlines()[0],
        "source_sha256": digest(source.encode()),
        "routine_sha256": {name: {"baseline": digest(extract(baseline,name).encode()), "updated": digest(extract(source,name).encode())}
                           for name in ["List", "AddBuffer"]},
        "test_sha256": {p.name:digest(p.read_bytes()) for p in [Path(__file__),HERE/"cases.f90"]},
        "compilation": compilation, "experiments": records, "passed": len(records), "median_seconds": medians,
        "elapsed_seconds": time.monotonic()-started,
        "scope": "Small extracted kernels, independent cell/image oracles and bitwise cross-variant/thread comparison; no production finder-speedup claim"}
    args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"passed":len(records),"elapsed_seconds":report["elapsed_seconds"],"median_seconds":medians},indent=2))


if __name__ == "__main__":
    main()
