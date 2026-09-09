"""Consolidate only this campaign's finished build and fixture trees.

Executables, active simulations, replay evidence, active launch bundles, reference
inputs and timestep tables remain in place. Every archived byte and symlink
is verified before removal; interrupted cleanup can resume from its receipt.
"""
import argparse
import getpass
import hashlib
import io
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import tarfile

from common import ROOT, WORK, now, sha, verify_manifest, write_json

PREFIXES=['work/native-build','work/ic-build','work/replay-build',
          'work/schedule-preflight','work/pilots']


def digest(stream):
    result=hashlib.sha256()
    while block:=stream.read(8*1024**2):result.update(block)
    return result.hexdigest()


def inspect(path):
    info=path.lstat()
    result=dict(bytes=info.st_size,mode=stat.S_IMODE(info.st_mode),mtime_ns=info.st_mtime_ns,
                device=info.st_dev,inode=info.st_ino)
    if stat.S_ISLNK(info.st_mode):result.update(kind='symlink',target=os.readlink(path))
    else:
        assert stat.S_ISREG(info.st_mode),str(path)
        result.update(kind='file',sha256=sha(path))
    return result


def verify_original(path,item):
    assert inspect(path)==item,f'Original changed: {path}'


def verify_archive(archive,manifest):
    with tarfile.open(archive,'r:gz') as tar:
        members=tar.getmembers()
        assert len(members)==len(manifest)+1
        assert {p.name for p in members}==set(manifest)|{'MANIFEST.json'}
        assert json.loads(tar.extractfile('MANIFEST.json').read())==manifest
        for name,item in manifest.items():
            member=tar.getmember(name);assert member.mode==item['mode']
            if item['kind']=='symlink':assert member.issym() and member.linkname==item['target']
            else:
                assert member.isfile() and member.size==item['bytes']
                with tar.extractfile(member) as stream:assert digest(stream)==item['sha256']


def active_job_ids():
    result=subprocess.run(['squeue','--noheader','--user',getpass.getuser(),
                           '--format=%A'],check=True,capture_output=True,text=True)
    return set(result.stdout.split())


def completed_launch_plan(jobs,states,active,archived):
    """Select exact registered launch paths; unknown/live jobs and files stay put."""
    manifest={};selected={}
    terminal={'COMPLETED','FAILED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL',
              'PREEMPTED','BOOT_FAIL','DEADLINE','CANCELLED'}
    for job in jobs:
        identifier=job['job_id'];state=states.get(identifier,'UNKNOWN')
        if state.split()[0] not in terminal or identifier in active:continue
        script=Path(job['script'])
        assert script.parent==WORK/'slurm' and script.parent.resolve()==WORK/'slurm',str(script)
        paths=[(script,job['script_sha256']),
               (script.with_suffix('.pyz'),job['bundle_sha256']),
               (script.parent/(script.stem+'-'+identifier+'.log'),None)]
        for path,expected in paths:
            name=str(path.relative_to(ROOT))
            if name in archived:
                if expected:assert archived[name]['sha256']==expected,name
                continue
            if not path.exists():
                assert expected is None,f'Missing unarchived frozen launch file: {path}'
                continue  # A job cancelled before starting need not have a log.
            assert path.is_file() and not path.is_symlink(),str(path)
            item=inspect(path)
            if expected:assert item['sha256']==expected,f'Changed frozen launch file: {path}'
            manifest[name]=item;selected[identifier]=state
    return dict(sorted(manifest.items())),selected


def archive_completed_launches(batch,require_finished=False):
    assert re.fullmatch(r'[a-z0-9][a-z0-9-]*',batch),'Use a simple archive batch label'
    if require_finished:
        analysis=json.loads((ROOT/'convergence.json').read_text())
        render=json.loads((ROOT/'render-validation.json').read_text())
        presentation=json.loads((ROOT/'presentation-validation.json').read_text())
        assert analysis['completed'] and len(analysis['inputs'])==21
        assert render['completed'] and render['complete_campaign']
        assert presentation['completed'] and presentation['complete_campaign']
        assert render['input_sha256']==presentation['input_sha256']==sha(ROOT/'convergence.json')
        assert render['presentation_manifest_sha256']==sha(ROOT/'presentation-validation.json')
        assert render['plot_manifest_sha256']==sha(ROOT/'plot-manifest-n300.json')
        assert sha(presentation['pdf'])==presentation['pdf_sha256']
        assert sha(render['figure_pdf'])==presentation['plot_pdf_sha256']
    archive=ROOT/f'launches-{batch}.tar.gz'
    receipt=ROOT/f'launches-cleanup-{batch}.json'
    if receipt.exists():
        record=json.loads(receipt.read_text())
        assert sha(archive)==record['archive_sha256']
        manifest=record['manifest'];verify_archive(archive,manifest)
        if record['completed']:
            print('Verified completed launch cleanup:',record['files_removed']);return
    else:
        import accounting
        accounting.main()
        jobs=json.loads((ROOT/'jobs.json').read_text())
        accounting_record=json.loads((ROOT/'accounting.json').read_text())
        states={job['job_id']:job['state'] for job in accounting_record['jobs']}
        archived={};prior_archives={}
        for previous in sorted(ROOT.glob('launches-cleanup-*.json')):
            prior=json.loads(previous.read_text())
            assert prior['completed'],f'Resume interrupted cleanup first: {previous}'
            assert sha(prior['archive'])==prior['archive_sha256']
            verify_archive(Path(prior['archive']),prior['manifest'])
            archived.update(prior['manifest'])
            prior_archives[str(previous)]=sha(previous)
        manifest,selected=completed_launch_plan(jobs,states,active_job_ids(),archived)
        if not manifest:
            print('No additional completed launch files to consolidate');return
        staged=archive.with_name(archive.name+'.tmp')
        assert not archive.exists() and not staged.exists(),'Preserve prior archive'
        assert sum(item['bytes'] for item in manifest.values())<=64*1024**2,'Resize the cleanup allocation for larger logs'
        payload=(json.dumps(manifest,indent=2)+'\n').encode()
        with staged.open('xb') as stream:
            with tarfile.open(fileobj=stream,mode='w:gz',compresslevel=4,dereference=False) as tar:
                member=tarfile.TarInfo('MANIFEST.json');member.size=len(payload)
                tar.addfile(member,io.BytesIO(payload))
                for name,item in manifest.items():
                    verify_original(ROOT/name,item);tar.add(ROOT/name,arcname=name,recursive=False)
            stream.flush();os.fsync(stream.fileno())
        verify_archive(staged,manifest)
        for name,item in manifest.items():verify_original(ROOT/name,item)
        assert not set(selected)&active_job_ids(),'A selected job became active; preserve all originals'
        staged.replace(archive)
        record=dict(completed=False,started_at_utc=now(),archive=str(archive),
                    archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,
                    archive_verified_member_by_member=True,manifest=manifest,
                    selected_terminal_jobs=selected,prior_archive_receipts=prior_archives,
                    accounting_sha256=sha(ROOT/'accounting.json'),files_removed=0,
                    require_finished=require_finished)
        write_json(receipt,record)
    assert not set(record['selected_terminal_jobs'])&active_job_ids(),'A selected job became active'
    present=[name for name in manifest if (ROOT/name).exists() or (ROOT/name).is_symlink()]
    for name in present:verify_original(ROOT/name,manifest[name])
    for name in present:
        verify_original(ROOT/name,manifest[name]);(ROOT/name).unlink()
    record.update(completed=True,completed_at_utc=now(),files_removed=len(manifest),
                  net_files_reduced=len(manifest)-2,
                  retained='All scientific inputs/outputs, executables, configurations, and active/unknown launch files',
                  restore=f'tar -xzf {archive.name} -C . (run in the campaign directory; restore only if needed)')
    write_json(receipt,record)
    print(f'Archived, verified and removed {len(manifest)} completed launch files')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--controls',action='store_true')
    parser.add_argument('--launches');parser.add_argument('--require-finished',action='store_true')
    args=parser.parse_args()
    if args.launches:
        assert not args.controls
        archive_completed_launches(args.launches,args.require_finished);return
    assert not args.require_finished
    prefixes=['work/driver-controls','work/driver-resume-controls-20260908'] if args.controls else PREFIXES
    archive=ROOT/('work-controls.tar.gz' if args.controls else 'work-artifacts.tar.gz')
    receipt=ROOT/('controls-cleanup.json' if args.controls else 'scratch-cleanup.json')
    if receipt.exists():
        record=json.loads(receipt.read_text());assert sha(archive)==record['archive_sha256']
        manifest=record['manifest']
        if record['completed']:
            print('Verified completed scratch cleanup:',record['files_removed']);return
        verify_archive(archive,manifest)
    else:
        for name in ['schedule-preflight.json','driver-controls.json','ic-production-validation.json',
                     'replay/preflight-ifx2024.json','ic/checks-ifx2024.json']:
            evidence=json.loads((ROOT/name).read_text())
            assert evidence.get('completed',evidence.get('all_passed',False)),name
        frozen=json.loads((ROOT/'executables.json').read_text())
        verify_manifest(frozen['binaries']);verify_manifest(frozen['build_receipts'])
        if args.controls:
            assert json.loads((ROOT/'replay-resume-controls.json').read_text())['completed']
        else:
            for pilot in (WORK/'pilots').iterdir():
                assert json.loads((pilot/'Run1/evolve.json').read_text())['completed'],str(pilot)
        manifest={}
        for prefix in prefixes:
            base=ROOT/prefix;assert base.is_dir() and not base.is_symlink(),str(base)
            for directory,dirs,files in os.walk(base,followlinks=False):
                links=[name for name in dirs if (Path(directory)/name).is_symlink()]
                dirs[:]=[name for name in dirs if name not in links]
                for name in sorted(files+links):
                    path=Path(directory)/name;assert path.parent.resolve().is_relative_to(base.resolve())
                    manifest[str(path.relative_to(ROOT))]=inspect(path)
        manifest=dict(sorted(manifest.items()))
        staged=archive.with_name(archive.name+'.tmp')
        assert not archive.exists() and not staged.exists(),'Preserve prior archive'
        payload=(json.dumps(manifest,indent=2)+'\n').encode()
        with staged.open('xb') as stream:
            with tarfile.open(fileobj=stream,mode='w:gz',compresslevel=4,dereference=False) as tar:
                member=tarfile.TarInfo('MANIFEST.json');member.size=len(payload)
                tar.addfile(member,io.BytesIO(payload))
                for name,item in manifest.items():
                    verify_original(ROOT/name,item);tar.add(ROOT/name,arcname=name,recursive=False)
            stream.flush();os.fsync(stream.fileno())
        verify_archive(staged,manifest)
        # No original is removed until every member and original is verified.
        for name,item in manifest.items():verify_original(ROOT/name,item)
        staged.replace(archive)
        record=dict(completed=False,started_at_utc=now(),archive=str(archive),archive_sha256=sha(archive),
                    archive_bytes=archive.stat().st_size,archive_verified_member_by_member=True,
                    prefixes=prefixes,manifest=manifest,files_removed=0)
        write_json(receipt,record)
    present=[name for name in manifest if (ROOT/name).exists() or (ROOT/name).is_symlink()]
    for name in present:verify_original(ROOT/name,manifest[name])
    for name in present:
        verify_original(ROOT/name,manifest[name]);(ROOT/name).unlink()
    for prefix in prefixes:
        for directory,dirs,files in os.walk(ROOT/prefix,topdown=False):
            path=Path(directory)
            if not any(path.iterdir()):path.rmdir()
    record.update(completed=True,completed_at_utc=now(),files_removed=len(manifest),
                  net_files_reduced=len(manifest)-2,
                  retained='work/bin, scientific IC/snapshots/replays, bundles/logs, reference inputs and timestep tables',
                  restore=f'tar -xzf {archive.name} -C . (run in the campaign directory; restore only if needed)')
    write_json(receipt,record)
    print(f'Archived, verified and removed {len(manifest)} finished scratch files/links')


if __name__=='__main__':main()
