"""Freeze the measured one-core comparison job after the main validation."""
import json
import shutil
import subprocess
import run_validation as v


def main():
    jobs = json.loads((v.ROOT / 'jobs.json').read_text())
    assert not any(j.get('stage') == 'main-plots' for j in jobs), 'Already submitted'
    parent = next(j for j in jobs if j.get('stage') == 'main-validation')
    pilot = v.ROOT / 'plotter-resource-pilot.json'
    plan = dict(partition='cosma8-serial', account='dp004', cpus=1, mem='2G',
                time_limit_seconds=600, expected_seconds_range=[60, 300],
                expected_core_hours_range=[1/60, 1/12], time_limit_core_hours=1/6,
                pilot_receipt_sha256=v.sha(pilot), pilot_elapsed_seconds=118.69,
                pilot_cpu_seconds=57.98, pilot_maxrss_kib=194484,
                memory_headroom_over_pilot=2*1024**2/194484,
                memory_reason='One-epoch N512 pilot used 189.93 MiB. 2 GiB provides 10.8x '
                    'headroom for three epochs and their larger selected populations. '
                    'Completed N1024 inline outputs contain 52,370 rows at z2 and '
                    '119,688 at z1, versus 71,587 pilot refined rows. Large raw tapes '
                    'are hashed in 8 MiB blocks and never loaded whole.',
                cache='Node-local temporary Matplotlib/LaTeX cache, removed on exit',
                full_node_exclusivity_required=False,
                dependencies=[f'afterok:{parent["job_id"]}'])
    v.write_json(v.ROOT / 'main-plots-resource-plan.json', plan)
    frozen = v.WORK / 'launch-main-plots'
    frozen.mkdir(exist_ok=False)
    for name in ['plot_main_results.py', 'compare_properties.py', 'house_style.py',
                 'chz-paper.mplstyle', 'plots.sbatch']:
        shutil.copy2(v.ROOT / name, frozen / name)
    command = ['sbatch', '--parsable', f'--dependency=afterok:{parent["job_id"]}',
               f'--export=ALL,BDM_VALIDATION_ROOT={v.ROOT}',
               str(frozen / 'plots.sbatch'), str(frozen / 'plot_main_results.py')]
    record = dict(stage='main-plots', purpose='Seven property PDFs from membership-verified '
                  'refined and original standalone catalogues on identical saved snapshots',
                  command=command, dependencies=plan['dependencies'],
                  script_sha256={p.name: v.sha(p) for p in frozen.iterdir()},
                  submitter_sha256=v.sha(__file__), resource_plan_sha256=v.sha(v.ROOT/'main-plots-resource-plan.json'),
                  commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=v.REPO, text=True).strip())
    v.write_json(frozen / 'submission.json', record)
    process = subprocess.run(command, cwd=v.REPO, text=True, capture_output=True)
    record.update(submitted_at_utc=v.now(), submission_returncode=process.returncode,
                  stdout=process.stdout, stderr=process.stderr)
    if process.returncode == 0:
        record['job_id'] = process.stdout.strip().split(';')[0]
        assert record['job_id'].isdigit()
        jobs.append(record)
        v.write_json(v.ROOT / 'jobs.json', jobs)
        record['scontrol'] = subprocess.check_output(['scontrol', 'show', 'job', record['job_id']], text=True)
    v.write_json(frozen / 'submission.json', record)
    process.check_returncode()
    print(record['job_id'])


if __name__ == '__main__':
    main()
