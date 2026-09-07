"""Archive verified generated root compiler files, leaving binaries and science data."""
import hashlib
import io
import json
import os
import stat
import subprocess
import tarfile

from common import REPO, ROOT, now, sha, write_json


def main():
    receipt=ROOT/'root-cleanup.json';archive=ROOT/'root-compiler-artifacts.tar.gz'
    if receipt.exists():
        previous=json.loads(receipt.read_text())
        assert sha(archive)==previous['archive_sha256']
        print('Verified prior root compiler cleanup:',previous['files_removed']);return
    tracked=set(subprocess.check_output(['git','ls-files','-z'],cwd=REPO).decode().split('\0'))
    candidates=[]
    for pattern in ['*.o','*.mod','*.opt.yaml','*.optrpt']:
        for path in sorted(REPO.glob(pattern)):
            if path.name in tracked or not stat.S_ISREG(path.lstat().st_mode):continue
            candidates.append(path)
    candidates=sorted(set(candidates))
    manifest={p.name:dict(bytes=p.stat().st_size,sha256=sha(p),mode=stat.S_IMODE(p.stat().st_mode),
                          mtime_ns=p.stat().st_mtime_ns,inode=p.stat().st_ino) for p in candidates}
    temporary=archive.with_name(archive.name+'.tmp')
    with temporary.open('xb') as stream:
        with tarfile.open(fileobj=stream,mode='w:gz') as tar:
            payload=(json.dumps(manifest,indent=2)+'\n').encode()
            info=tarfile.TarInfo('MANIFEST.json');info.size=len(payload)
            tar.addfile(info,io.BytesIO(payload))
            for path in candidates:tar.add(path,arcname=path.name,recursive=False)
        stream.flush();os.fsync(stream.fileno())
    with tarfile.open(temporary,'r:gz') as tar:
        assert sorted(tar.getnames())==sorted([*manifest,'MANIFEST.json'])
        for name,item in manifest.items():
            member=tar.getmember(name);assert member.isfile() and member.size==item['bytes']
            digest=hashlib.sha256()
            with tar.extractfile(member) as stream:
                while chunk:=stream.read(8*1024**2):digest.update(chunk)
            assert digest.hexdigest()==item['sha256']
    # Validate every original before deleting any; no unrelated directory walk.
    for path in candidates:
        item=manifest[path.name];info=path.stat()
        assert info.st_size==item['bytes'] and info.st_mtime_ns==item['mtime_ns'] and info.st_ino==item['inode']
        assert sha(path)==item['sha256']
    os.replace(temporary,archive)
    record=dict(started_at_utc=now(),archive=str(archive),archive_sha256=sha(archive),
                archive_verified_member_by_member=True,manifest=manifest,files_removed=0)
    write_json(receipt,record)
    for path in candidates:path.unlink()
    record.update(completed_at_utc=now(),files_removed=len(candidates),net_files_reduced=len(candidates)-2,
                  preserved='Executables, source, tracked files, simulation data and all unrelated directories')
    write_json(receipt,record)
    print(f'Archived and verified {len(candidates)} generated compiler files; net reduction {len(candidates)-2}')


if __name__=='__main__':main()
