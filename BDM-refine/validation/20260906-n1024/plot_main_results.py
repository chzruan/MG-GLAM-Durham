"""Render property comparisons only from completed, hash-verified main replays.

Run with micromamba run -n cosemu python3 -B under the sized Slurm allocation.
The fixed-density conclusion is added to the slides after its separate job.
"""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import traceback

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np

ROOT = Path(os.environ.get('BDM_VALIDATION_ROOT', Path(__file__).resolve().parent)).resolve()
SCRIPT_DIR = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block := stream.read(8 * 1024**2):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Production catalogue matching requires its measured allocation'
    simulation_path = ROOT / 'main-simulation.json'
    validation_path = ROOT / 'main-validation.json'
    preparation_path = ROOT / 'main-preparation.json'
    simulation = json.loads(simulation_path.read_text())
    validation = json.loads(validation_path.read_text())
    preparation = json.loads(preparation_path.read_text())
    assert simulation['completed'] and validation['completed']
    assert validation['physics_identity_checks_complete']
    assert validation['simulation_receipt_sha256'] == sha(simulation_path)
    assert simulation['spec'] == validation['spec']
    assert simulation['spec']['sources_sha256'] == preparation['sources_sha256']
    comparison_source = ROOT / 'main-comparison-catalogues.npz'
    assert validation['comparison_sha256'] == sha(comparison_source)
    assert {s['z'] for s in simulation['snapshots']} == {0, 1, 2}
    # The normal density accumulation can change between calls. Link every
    # plotted refined table to the exact membership replay that was checked,
    # rather than transferring a raw-membership claim to an inline catalogue.
    catalogues = ROOT / 'main-verified-comparison-catalogues.npz'
    assert not catalogues.exists(), 'Preserve an existing verified comparison'
    arrays, membership_provenance = {}, {}
    with np.load(comparison_source, allow_pickle=False) as archive:
        for z in [0, 1, 2]:
            members = [r for r in validation['records'] if r['z'] == z and r['variant'] == 'members']
            old = [r for r in validation['records'] if r['z'] == z and r['variant'] == 'old']
            assert len(members) == len(old) == 1
            members, old = members[0], old[0]
            assert members['scientific_verification_completed'] and old['scientific_verification_completed']
            membership = members['membership']
            assert all(membership[key] == 0 for key in ['exact_duplicate_member_sets',
                'repeated_original_ids', 'mass_count_mismatches', 'host_exclusion_violations'])
            raw = ROOT / membership['retained_raw']
            assert sha(raw) == membership['raw_sha256']
            assert sha(raw.with_suffix('.index.npz')) == membership['index_sha256']
            for record in [old, members]:
                assert sha(ROOT / record['catalogue']['path']) == record['catalogue']['sha256']
            arrays[f'old_z{z}'] = archive[f'old_z{z}']
            assert np.array_equal(arrays[f'old_z{z}'], np.loadtxt(
                ROOT / old['catalogue']['path'], skiprows=8, ndmin=2), equal_nan=True)
            arrays[f'new_z{z}'] = np.loadtxt(ROOT / members['catalogue']['path'], skiprows=8, ndmin=2)
            assert len(arrays[f'new_z{z}']) == membership['selected']
            membership_provenance[f'z{z}'] = dict(catalogue=members['catalogue'],
                membership=membership, normal_density_comparison_to_inline=members['normal_density_comparison'],
                byte_identical_to_inline=members['byte_identical_to_inline'])
    np.savez_compressed(catalogues, **arrays)
    del arrays
    metadata = dict(simulation['spec'],
                    snapshot_provenance=simulation['snapshots'],
                    simulation_receipt_sha256=sha(simulation_path),
                    validation_receipt_sha256=sha(validation_path),
                    preparation_receipt_sha256=sha(preparation_path),
                    source_comparison_npz_sha256=sha(comparison_source),
                    plotted_refined_catalogues=membership_provenance,
                    plotted_refined_scope='Standalone refined finder with a read-only membership '
                        'dump after publication. Raw identity, mass/count and host checks apply '
                        'directly to these plotted refined tables. Inline catalogues are retained '
                        'separately; ordinary density accumulation can produce differences.',
                    sample_label='Same-snapshot finder comparison',
                    reference_logmass=12.5,
                    reference_source='https://arxiv.org/pdf/2110.00328, Appendix A',
                    reference_scope='Published z=0 abundance-resolution guidance only; '
                        'the line is a common visual reference at all three epochs, '
                        'not new convergence evidence for this finder or z1/z2',
                    old_baseline_scope='Independent standalone processes retain the '
                        'historical first-call configuration-order defect',
                    old_finder_binary_sha256=preparation['binaries_sha256']['PMP2BDM.old.exe'],
                    new_finder_binary_sha256=preparation['binaries_sha256']['PMP2BDM.members.exe'],
                    new_inline_binary_sha256=preparation['binaries_sha256']['PMP2main.exe'])
    metadata_path = ROOT / 'comparison-metadata.json'
    output = ROOT / 'slides/figs'
    assert not output.exists(), 'Preserve existing figure evidence; inspect before rerunning'
    write_json(metadata_path, metadata)
    report = dict(completed=False, started_at_utc=datetime.now(timezone.utc).isoformat(),
                  job_id=os.environ['SLURM_JOB_ID'], runner_sha256=sha(__file__),
                  metadata_sha256=sha(metadata_path), catalogue_sha256=sha(catalogues))
    log = ROOT / 'work/main-property-plots.log'
    timing = ROOT / 'work/main-property-plots.time'
    start = time.monotonic()
    try:
        with tempfile.TemporaryDirectory(prefix='main-plot-cache-') as cache:
            env = dict(os.environ, MPLCONFIGDIR=cache, OPENBLAS_NUM_THREADS='1',
                       MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
            command = ['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(timing),
                       'micromamba', 'run', '-n', 'cosemu', 'python3', '-B',
                       str(SCRIPT_DIR / 'compare_properties.py'), '--catalogues', str(catalogues),
                       '--metadata', str(metadata_path), '--outdir', str(output), '--expect-ordered-new']
            report['command'] = command
            with log.open('x') as stream:
                process = subprocess.run(command, cwd=ROOT, env=env, stdout=stream,
                                         stderr=subprocess.STDOUT, timeout=1800)
            process.check_returncode()
        elapsed, user, system, rss = timing.read_text().strip().splitlines()[-1].split()
        receipt = json.loads((output / 'comparison-figure-receipt.json').read_text())
        assert len(receipt['figures_sha256']) == 7
        for name, expected in receipt['figures_sha256'].items():
            assert sha(output / name) == expected
        assert receipt['script_sha256'] == sha(SCRIPT_DIR / 'compare_properties.py')
        report.update(completed=True, elapsed_seconds=time.monotonic()-start,
                      timed_elapsed_seconds=float(elapsed), cpu_user_seconds=float(user),
                      cpu_system_seconds=float(system), maxrss_kib=int(rss),
                      log_sha256=sha(log), figure_receipt_sha256=sha(output/'comparison-figure-receipt.json'),
                      finished_at_utc=datetime.now(timezone.utc).isoformat())
    except BaseException:
        report['failure_traceback'] = traceback.format_exc()
        raise
    finally:
        write_json(ROOT / 'main-property-plots.json', report)
    print('Seven main comparison PDFs completed and hash-verified', flush=True)


if __name__ == '__main__':
    main()
