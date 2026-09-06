"""Consolidate only completed follow-up scratch after verified archive readback.

Use micromamba run -n cosemu python3 -B under a one-core shared allocation.
Retain production membership/density tapes and indexes. Never follow snapshot
symlinks into the earlier validation campaign or modify its archives.
"""
from datetime import datetime, timezone
import importlib.util
import io
import json
import os
from pathlib import Path
import tarfile

ROOT=Path(os.environ.get('BDM_REVIEW_ROOT',Path(__file__).resolve().parent)).resolve()
REPO=ROOT.parents[2]


def main():
    assert os.environ.get('SLURM_JOB_ID'),'Archive compression requires its shared allocation'
    definition=importlib.util.spec_from_file_location('archive_helpers',
        REPO/'BDM-refine/validation/20260906-n1024/archive_work.py')
    utility=importlib.util.module_from_spec(definition)
    definition.loader.exec_module(utility)
    summary=ROOT/'native/comparison.json'
    accounting=ROOT/'native/accounting.json'
    assert json.loads(summary.read_text())['completed']
    assert json.loads(accounting.read_text())['all_completed']
    archive=ROOT/'work-artifacts.tar.gz'
    staged=archive.with_name(archive.name+'.tmp')
    receipt=ROOT/'work-archive.json'
    assert not archive.exists() and not staged.exists() and not receipt.exists()
    live_log=ROOT/f'work/archive-{os.environ["SLURM_JOB_ID"]}.log'
    files,retained=utility.scan(ROOT,['work'],live_log)
    # The shared helper conservatively retains every PM header/member tape.
    # Here the native preflight fixtures are disposable small test data; only
    # actual N1024 replay data and the running job log need remain unpacked.
    for name in list(retained):
        if not name.startswith('work/replays-t') and ROOT/name!=live_log:
            item=retained.pop(name)
            assert item['bytes']<256*1024**2,'Review an unexpected large test fixture'
            item.pop('reason',None)
            item['sha256']=utility.sha(ROOT/name)
            files[name]=item
    files=dict(sorted(files.items()))
    manifest=dict(schema_version=1,original_root=str(ROOT),prefixes=['work'],archived=files,retained=retained,
                  comparison_sha256=utility.sha(summary),accounting_sha256=utility.sha(accounting))
    manifest_bytes=(json.dumps(manifest,indent=2)+'\n').encode()
    report=dict(completed=False,job_id=os.environ['SLURM_JOB_ID'],
                started_at_utc=datetime.now(timezone.utc).isoformat(),archiver_sha256=utility.sha(__file__),
                helper_sha256=utility.sha(definition.origin),manifest=manifest,removed_files=[],removed_directories=[])
    try:
        with tarfile.open(staged,'w:gz',compresslevel=6,dereference=False) as tar:
            for name,item in files.items():
                utility.verify_original(ROOT/name,item)
                tar.add(ROOT/name,arcname=name,recursive=False)
            member=tarfile.TarInfo('MANIFEST.json')
            member.size=len(manifest_bytes);member.mode=0o644
            tar.addfile(member,io.BytesIO(manifest_bytes))
        with staged.open('rb') as stream:os.fsync(stream.fileno())
        with tarfile.open(staged,'r:gz') as tar:
            members={m.name:m for m in tar.getmembers()}
            assert set(members)==set(files)|{'MANIFEST.json'}
            assert tar.extractfile(members['MANIFEST.json']).read()==manifest_bytes
            for name,item in files.items():
                member=members[name]
                assert member.mode==item['mode']
                if item['kind']=='symlink':
                    assert member.issym() and member.linkname==item['target']
                else:
                    assert member.isfile() and member.size==item['bytes']
                    assert utility.digest(tar.extractfile(member))==item['sha256']
        staged.replace(archive)
        report.update(archive_sha256=utility.sha(archive),archive_bytes=archive.stat().st_size,
                      archive_verified_before_removal=True)
        utility.write_json(receipt,report)
        for name,item in files.items():utility.verify_original(ROOT/name,item)
        for name,item in files.items():
            utility.verify_original(ROOT/name,item)
            (ROOT/name).unlink()
            report['removed_files'].append(name)
        for directory,dirs,names in os.walk(ROOT/'work',topdown=False):
            path=Path(directory)
            if not any(path.iterdir()):
                path.rmdir();report['removed_directories'].append(str(path.relative_to(ROOT)))
        for name,item in retained.items():
            path=ROOT/name
            if path==live_log:continue
            stat=path.lstat()
            assert (stat.st_dev,stat.st_ino,stat.st_mtime_ns,stat.st_size)==(
                item['device'],item['inode'],item['mtime_ns'],item['bytes'])
        report.update(completed=True,finished_at_utc=datetime.now(timezone.utc).isoformat(),
                      archived_file_count=len(files),retained_file_count=len(retained),
                      removed_directory_count=len(report['removed_directories']))
    finally:
        utility.write_json(receipt,report)
    print('Verified archive; removed',len(files),'finished scratch files and',len(report['removed_directories']),'directories')


if __name__=='__main__':
    main()
