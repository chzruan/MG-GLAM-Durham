"""Check unbinding work counters against physical fixtures and unmodified code.

Use micromamba run -n cosemu python3 -B. Temporary native builds are removed;
the original independent review and its results remain unchanged.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
TESTS = REPO/'BDM-refine/repairs/20260906/tests'
sys.path.insert(0, str(TESTS))
import numpy as np
import halo_review as review
import halo_gaps as gaps


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def program(text, measured):
    # Equal-mass fixture consistent with the chosen box/grid mean density.
    text = text.replace('Box=32.;NGRID=128;Np=n;', 'Box=32.;NGRID=128;NROW=16;Np=n;')
    if measured:
        needle = "  write(*,*) 'VALUES'"
        assert text.count(needle) == 1
        text = text.replace(needle, "  write(*,*) 'UNBINDING',HaloUnbindingPasses(1),HaloUnbindingWork(1)\n"+needle)
        needle = "  write(*,*) 'EMPTY_RECALL_PASS'"
        assert text.count(needle) == 1
        text = text.replace(needle, "  if(HaloUnbindingPasses(1)/=0.or.HaloUnbindingWork(1)/=0_8) &\n"
            "    error stop 'stale unbinding diagnostics after empty call'\n"+needle)
    return text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'diagnostics-results.json')
    # The uninstrumented SO-v3 commit is the like-for-like arithmetic control.
    # The first retained receipt explicitly records its earlier v2 baseline.
    parser.add_argument('--baseline', default='e2a327a5fdb5de5c8dc7e0fe30c700fb147f2df4')
    args = parser.parse_args()
    source = (REPO/'PMP2linker.f90').read_text()
    before = subprocess.check_output(['git', 'show', args.baseline+':PMP2linker.f90'], cwd=REPO, text=True)
    mass = float(np.float32(float(np.float32(.3))*2.774e11*(32/16)**3))
    fixtures = [
        ('empty', np.empty((0,6)), 0, 0),
        ('one_pass', np.c_[review.sphere(128,.1)+5., np.zeros((128,3))], 1, 128),
        ('all_unbound', np.c_[review.sphere(64,.1)+5., review.sphere(64,1.e6)], 1, 64),
        ('singular_centre', np.c_[np.full((20,3),5.), np.zeros((20,3))], 1, 20),
    ]
    for pairs, removed in [(100,90),(200,190)]:
        phase = gaps.energy_ladder(pairs, removed, mass=mass)
        passes, survivors = gaps.reference_unbinding(phase, mass)
        assert passes == removed+1 and len(survivors) == 20
        fixtures.append((f'ladder_{2*pairs}',phase,passes,sum(range(20,2*pairs+1,2))))
    records, builds = [], {}
    with tempfile.TemporaryDirectory(prefix='bdm-unbinding-diagnostics-') as directory:
        root = Path(directory)
        for variant, text in [('before',before),('measured',source)]:
            work = root/variant
            work.mkdir()
            review.SOURCE = text
            builds[variant] = review.build(work, lambda p: program(p,variant=='measured'))
        for name, phase, passes, evaluations in fixtures:
            for mode in ['checked','optimized']:
                old = review.run_halo(root/'before',name,phase,mode,mass=mass,identity=True)
                new = review.run_halo(root/'measured',name,phase,mode,mass=mass,identity=True)
                assert old['ids'] == new['ids'] and old['values'] == new['values'], (name,mode)
                found = re.search(r'UNBINDING\s+(\d+)\s+(\d+)',new['stdout'])
                assert found and tuple(map(int,found.groups())) == (passes,evaluations), (name,new['stdout'])
                records.append(dict(case=name,mode=mode,rows=len(phase),passes=passes,
                    active_particle_rows=evaluations,survivors=len(new['ids']),
                    all_properties_and_membership_unchanged=True,empty_recall_clears_counters=True))
    report = dict(completed=True,created_at_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(source),baseline=args.baseline,baseline_source_sha256=sha(before),
        test_sha256=sha(Path(__file__).read_text()),mass_one=mass,fixture_nrow=16,fixture_box=32,
        comparison_count=len(records),native_executions=2*len(records),compilation=builds,records=records)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(f'{len(records)} unbinding-counter comparisons passed; all {2*len(records)} native executions succeeded')


if __name__ == '__main__':
    main()
