"""Consolidate only this campaign's finished build and fixture trees.

Executables, active simulations, replay evidence, launch bundles, reference
inputs and timestep tables remain in place. Every archived byte and symlink
is verified before removal; interrupted cleanup can resume from its receipt.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import stat
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


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--controls',action='store_true');args=parser.parse_args()
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
