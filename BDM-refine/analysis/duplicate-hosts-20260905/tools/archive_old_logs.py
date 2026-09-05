"""Archive reviewed historical logs, verify every byte, then remove loose copies.

Only entries in the saved review manifest are eligible. Scientific .DAT files,
tracked files, symlinks, and files changed since review are rejected.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import time

from catalogue_core import sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('manifest',type=Path)
    p.add_argument('archive',type=Path)
    p.add_argument('receipt',type=Path)
    mode=p.add_mutually_exclusive_group()
    mode.add_argument('--prepare-only',action='store_true',help='Verify archive but keep all loose files')
    mode.add_argument('--remove-verified',action='store_true',help='Reverify a prepared archive before unlinking')
    a=p.parse_args();root=a.root.resolve()
    tracked=set(subprocess.check_output(['git','-C',str(root),'ls-files','-z']).decode().split('\0'))
    if a.remove_verified:
        with tarfile.open(a.archive,'r:gz') as archive:
            entries=json.load(archive.extractfile('CLEANUP-MANIFEST.json'))
    else:
        entries=json.loads(a.manifest.read_text())
    for entry in entries:
        rel=entry['path'];path=root/rel
        if rel in tracked or path.is_symlink() or not path.resolve().is_relative_to(root):
            raise ValueError(f'Unsafe archive target {rel}')
        if path.suffix.lower()=='.dat':
            raise ValueError(f'Scientific data is not an old log: {rel}')
        stat=path.stat()
        if (stat.st_size,stat.st_mtime_ns)!=(entry['bytes'],entry['mtime_ns']):
            raise ValueError(f'Changed since cleanup review: {rel}')
        digest=sha256(path)
        if a.remove_verified and digest!=entry['sha256']:
            raise ValueError(f'Changed since archive preparation: {rel}')
        entry['sha256']=digest
    a.archive.parent.mkdir(parents=True,exist_ok=True)
    if not a.remove_verified:
        with tarfile.open(a.archive,'x:gz',compresslevel=6) as archive:
            for entry in entries:
                archive.add(root/entry['path'],arcname=entry['path'],recursive=False)
            payload=json.dumps(entries,indent=2).encode()
            info=tarfile.TarInfo('CLEANUP-MANIFEST.json');info.size=len(payload);info.mtime=int(time.time())
            archive.addfile(info,io.BytesIO(payload))
    with tarfile.open(a.archive,'r:gz') as archive:
        for entry in entries:
            source=archive.extractfile(entry['path'])
            h=hashlib.sha256()
            for block in iter(lambda:source.read(1<<20),b''):
                h.update(block)
            if h.hexdigest()!=entry['sha256']:
                raise ValueError(f'Archive verification failed: {entry["path"]}')
    # Recheck every source before the first deletion, then each one at unlink.
    for entry in entries:
        if sha256(root/entry['path'])!=entry['sha256']:
            raise ValueError(f'Source changed while archiving: {entry["path"]}')
    if not a.prepare_only:
        for entry in entries:
            path=root/entry['path'];stat=path.stat()
            if (stat.st_size,stat.st_mtime_ns)!=(entry['bytes'],entry['mtime_ns']):
                raise ValueError(f'Source changed before unlink: {entry["path"]}')
            path.unlink()
    receipt=dict(archived_files=len(entries),removed_loose_files=0 if a.prepare_only else len(entries),
                 prepared_only=a.prepare_only,
                 original_bytes=sum(e['bytes'] for e in entries),
                 archive=str(a.archive.resolve()),archive_bytes=a.archive.stat().st_size,
                 archive_sha256=sha256(a.archive),all_members_sha256_verified=True,
                 manifest='CLEANUP-MANIFEST.json inside archive',
                 restore_command=f'tar -xzf {a.archive.resolve()} -C {root}',
                 scope='Reviewed historical text logs only; no science catalogues or tracked files')
    with open(a.receipt,'x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    main()
