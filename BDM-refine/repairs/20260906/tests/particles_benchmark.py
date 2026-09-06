"""Compare particle slab prefiltering with the audited List routine.

Run using micromamba run -n cosemu python3 -B particles_benchmark.py.
This small kernel experiment is not a production-size scaling result.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

from particles_tests import HERE, REPO, build, extract


def main():
    source = (REPO/'PMP2linker.f90').read_text()
    audited = subprocess.check_output(['git', 'show', 'a8c7715:PMP2linker.f90'], cwd=REPO, text=True)
    old_list = extract(audited, 'List')
    source_old = source.replace(extract(source, 'List'), old_list)
    samples = []
    compilation = {}
    with tempfile.TemporaryDirectory(prefix='bdm-particle-list-benchmark-') as scratch:
        work = Path(scratch)
        for label, text in [('baseline', source_old), ('prefilter', source)]:
            destination = work/label
            destination.mkdir()
            compilation[label] = build(destination, text, 'optimized')
        for threads in [1, 2, 4, 8]:
            for repeat in range(3):
                for label in (['baseline', 'prefilter'] if repeat % 2 == 0 else ['prefilter', 'baseline']):
                    env = {**os.environ, 'OMP_NUM_THREADS': str(threads), 'OMP_DYNAMIC': 'FALSE',
                           'OMP_PROC_BIND': 'false', 'OPENBLAS_NUM_THREADS': '1'}
                    result = subprocess.run([str(work/label/'optimized'), 'list_benchmark', '128'],
                                            cwd=work, env=env, text=True, capture_output=True, timeout=30)
                    if result.returncode or 'PARTICLES_TEST_PASS' not in result.stdout:
                        raise RuntimeError(result.stdout+'\n'+result.stderr)
                    seconds = float(re.search(r'LIST_SECONDS_COUNT\s+([\d.Ee+-]+)', result.stdout).group(1))
                    samples.append(dict(variant=label, threads=threads, repeat=repeat,
                                        seconds_per_build=seconds, stdout=result.stdout, stderr=result.stderr))
    medians = {}
    for threads in [1, 2, 4, 8]:
        times = {label: statistics.median(s['seconds_per_build'] for s in samples
                                         if s['threads'] == threads and s['variant'] == label)
                 for label in ['baseline', 'prefilter']}
        times['baseline_over_prefilter_speedup'] = times['baseline']/times['prefilter']
        medians[str(threads)] = times
    result = dict(particles=128**3, rebuilds_per_sample=10, samples_per_variant=3,
                  source_sha256=hashlib.sha256(source.encode()).hexdigest(),
                  old_list_sha256=hashlib.sha256(old_list.encode()).hexdigest(),
                  benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  compilation=compilation, median_seconds=medians, samples=samples,
                  correctness='Every link and cell head equals an independent serial row-order oracle; every row reachable exactly once',
                  limitations='Small kernel timing on the current host; shared resource, allocation and cache noise; not a production finder benchmark')
    (HERE/'particles_benchmark_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(medians, indent=2))


if __name__ == '__main__':
    main()
