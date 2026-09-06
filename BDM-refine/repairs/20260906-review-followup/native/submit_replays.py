"""Seal source/build/input provenance and submit two measured shared jobs.

Execute with micromamba run -n cosemu python3 -B. Requires a completed native
build and its small preflight. Does not evolve a simulation or alter original
snapshots, density tapes, validation receipts, or catalogues.
"""
import hashlib
import inspect
import json
import os
from pathlib import Path
import shutil
import subprocess

import run_replays as r


def main():
    here=r.HERE
    assert not (here/'plan.json').exists(),'Refuse to overwrite a frozen plan'
    build=json.loads((here/'build.json').read_text())
    assert build['completed']
    for variant in ['reference','optimized-v2','v3']:
        assert r.sha(build['variants'][variant]['binary_path'])==build['variants'][variant]['binary_sha256']
    simulation_path=r.VALIDATION/'main-simulation.json'
    thread_path=r.VALIDATION/'main-thread-control/results.json'
    simulation=json.loads(simulation_path.read_text())
    threading=json.loads(thread_path.read_text())
    assert simulation['completed'] and threading['completed'] and threading['fixed_density_controls_passed']
    spec={key:simulation['spec'][key] for key in
          ['nrow','ngrid','box_mpc_h','cosmology','halo_config','particle_limit_exclusive']}
    assert (spec['nrow'],spec['ngrid'],spec['box_mpc_h'])==(1024,2048,512.)
    assert spec['nrow']**3<1200**3
    snapshots=[dict(z=s['z'],header=s['header'],files_sha256=s['files_sha256'],
                    directory=str(r.VALIDATION/'work/main/Run1')) for s in simulation['snapshots']]
    for snapshot in snapshots:
        for name in snapshot['files_sha256']:
            assert (Path(snapshot['directory'])/name).is_file()
    density=r.VALIDATION/'main-thread-control/work/run/d64.density.bin'
    density_receipt=threading['density_tapes'][density.name]
    assert density.stat().st_size==density_receipt['bytes']
    groups={}
    for threads in [32,64]:
        variants=['optimized-v2','reference','v3'] if threads==32 else ['reference','optimized-v2','v3']
        groups[str(threads)]=[dict(label='z0-'+variant,z=0,variant=variant,passes=2,
                                  density_mode='read',density_path=str(density),
                                  density_sha256=density_receipt['sha256']) for variant in variants]
    for z in [2,1]:
        field=r.ROOT/f'work/replays-t64/z{z}.density.bin'
        groups['64'].extend([
            dict(label=f'z{z}-reference',z=z,variant='reference',passes=1,density_mode='write',density_path=str(field)),
            dict(label=f'z{z}-v3',z=z,variant='v3',passes=1,density_mode='read',density_path=str(field))])
    command=['sacct','-j','11948467,11948491','--format=JobID,State,ElapsedRaw,TotalCPU,ReqTRES%100,AllocTRES%100,MaxRSS,ExitCode','--parsable2']
    accounting=subprocess.check_output(command,text=True)
    assert 'COMPLETED' in accounting and '198231384K' in accounting
    resources=dict(partition='cosma8-serial',account='dp004',nodes=1,ntasks=1,exclusive=False,
                   mem='256G',time_limit_seconds=3600,
                   pilot_jobs=['11948467','11948491'],pilot_accounting_command=command,pilot_accounting=accounting,
                   measured_fixed_density_batch_maxrss_kib=198231384,
                   measured_fixed_density_batch_gib=198231384/1024**2,
                   measured_fixed_density_native_gib=152.,
                   memory_reason='The same N1024 saved-PM-bit and immutable-FI wrapper previously used '
                     '189.05 GiB batch RSS. 256 GiB gives 35% headroom for variant/tape I/O. '
                     'The optimized buffer count reduces scratch by about 7 GiB; this saving is not needed to fit the request.',
                   packing='Independent 32-core and 64-core jobs use their actual thread counts on the shared '
                     'queue. Concurrent requests total 96 cores and 512 GiB, below one COSMA8 node; '
                     'the scheduler may place them separately. No exclusive-node request.',
                   groups={'32':dict(cpus=32,expected_seconds_range=[1200,2100],expected_core_hours_range=[32*1200/3600,32*2100/3600],time_limit_core_hours=32),
                           '64':dict(cpus=64,expected_seconds_range=[1200,2100],expected_core_hours_range=[64*1200/3600,64*2100/3600],time_limit_core_hours=64)})
    v=r.checker()
    plan=dict(created_at_utc=r.now(),completed_build=True,build_sha256=r.sha(here/'build.json'),
              source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=r.REPO,text=True).strip(),
              driver_sha256=r.sha(here/'run_replays.py'),submitter_sha256=r.sha(__file__),
              batch_sha256=r.sha(here/'replays.sbatch'),spec=spec,snapshots=snapshots,groups=groups,
              historical_z0_d64_catalogue_sha256=threading['fixed_density_comparisons'][0]['first_sha256'],
              resources=resources,checker_function_sha256=hashlib.sha256(inspect.getsource(v.check_memberships).encode()).hexdigest(),
              evidence_sha256={str(path):r.sha(path) for path in
                  [simulation_path,thread_path,r.VALIDATION/'run_validation.py']},
              interpretation='Reference is pre-follow-up v2 plus output-invariant telemetry. Optimized-v2 '
                'uses final computational code with only the diagnostic SO coefficient/header reverted. '
                'Final v3 changes the SO normalization. Each physical pair consumes exactly the same FI '
                'and original snapshot. Fixed-FI reproducibility does not imply ordinary-density peak stability.',
              performance_scope='Two stage timings per variant/thread setting at z0; report individual values '
                'and medians as a controlled replay measurement, not a general scaling law. Last-pass full '
                'finder time includes membership/diagnostic I/O and must not be used as production wall time.')
    launch=r.ROOT/'work/launch-native'
    launch.mkdir(exist_ok=False)
    for name in ['run_replays.py','replays.sbatch']:
        shutil.copy2(here/name,launch/name)
        assert r.sha(launch/name)==r.sha(here/name)
    r.write_json(here/'plan.json',plan)
    submissions=[]
    for threads in [64,32]:
        command=['sbatch','--parsable','--partition=cosma8-serial','--account=dp004',
                 '--nodes=1','--ntasks=1',f'--cpus-per-task={threads}','--mem=256G','--time=01:00:00',
                 f'--job-name=bdm-v3-t{threads}',f'--chdir={r.REPO}',
                 f'--output={r.ROOT}/work/native-t{threads}-%j.log',
                 f'--error={r.ROOT}/work/native-t{threads}-%j.log',
                 '--export=ALL,BDM_REVIEW_ROOT='+str(r.ROOT),str(launch/'replays.sbatch'),
                 str(launch/'run_replays.py'),r.sha(here/'plan.json')]
        record=dict(command=command,dependencies=[],threads=threads,prepared_at_utc=r.now(),
                    plan_sha256=r.sha(here/'plan.json'),script_sha256={p.name:r.sha(p) for p in launch.iterdir()})
        process=subprocess.run(command,capture_output=True,text=True)
        record.update(submitted_at_utc=r.now(),returncode=process.returncode,stdout=process.stdout,stderr=process.stderr)
        if process.returncode==0:
            record['job_id']=process.stdout.strip().split(';')[0]
            state=subprocess.run(['scontrol','show','job',record['job_id'],'-o'],capture_output=True,text=True)
            record['scontrol']=state.stdout
        submissions.append(record)
        r.write_json(here/'submission.json',submissions)
        assert process.returncode==0,process.stderr
        print('Submitted',record['job_id'],'with',threads,'shared cores',flush=True)


if __name__=='__main__':
    main()
