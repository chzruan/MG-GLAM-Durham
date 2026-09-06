"""Verify an archive before removing only this experiment's finished scratch.

Run with micromamba run -n cosemu python3 -B under Slurm after presentation
verification. Keep large particle/density/membership tapes and their PM headers
in place; consolidate the numerous small build, replay and launch files.
"""
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import stat
import tarfile

ROOT = Path(os.environ.get('BDM_VALIDATION_ROOT', Path(__file__).resolve().parent)).resolve()
PREFIXES = ['work', 'numerical-threading/work', 'main-thread-control/work']


def digest(stream):
    value = hashlib.sha256()
    while block := stream.read(8*1024**2):
        value.update(block)
    return value.hexdigest()


def sha(path):
    with Path(path).open('rb') as stream:
        return digest(stream)


def write_json(path, value):
    staged = path.with_name(path.name+'.tmp')
    with staged.open('w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    staged.replace(path)


def scan(root, prefixes, live_log):
    archived, retained = {}, {}
    for prefix in prefixes:
        base = root / prefix
        assert base.is_dir() and not base.is_symlink(), f'Expected scratch root has not been consolidated: {base}'
        for directory, dirs, files in os.walk(base, followlinks=False):
            links = [name for name in dirs if (Path(directory)/name).is_symlink()]
            dirs[:] = [name for name in dirs if name not in links]
            for name in sorted(files+links):
                path = Path(directory)/name
                assert path.parent.resolve().is_relative_to(root.resolve())
                relative = str(path.relative_to(root))
                info = path.lstat()
                item = dict(bytes=info.st_size, mode=stat.S_IMODE(info.st_mode),
                            mtime_ns=info.st_mtime_ns, device=info.st_dev, inode=info.st_ino)
                if path.is_symlink():
                    item.update(kind='symlink', target=os.readlink(path))
                    archived[relative] = item
                else:
                    assert path.is_file(), relative
                    item['kind'] = 'file'
                    if (info.st_size >= 256*1024**2 or
                        (name.startswith('PMcr') and name.endswith('.DAT')) or
                        name == 'repair-members.bin' or
                        name.endswith('.index.npz') or path == live_log):
                        item['reason'] = 'Large retained data, associated header/index, or live archive-job log'
                        retained[relative] = item
                    else:
                        item['sha256'] = sha(path)
                        archived[relative] = item
    return dict(sorted(archived.items())), dict(sorted(retained.items()))


def verify_original(path, record):
    info = path.lstat()
    assert stat.S_IMODE(info.st_mode) == record['mode'], str(path)
    assert (info.st_dev, info.st_ino, info.st_mtime_ns, info.st_size) == (
        record['device'], record['inode'], record['mtime_ns'], record['bytes']), str(path)
    if record['kind'] == 'symlink':
        assert stat.S_ISLNK(info.st_mode) and os.readlink(path) == record['target']
    else:
        assert stat.S_ISREG(info.st_mode) and sha(path) == record['sha256']


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Archive compression requires its sized allocation'
    summary = ROOT/'validation-summary.json'
    presentation = ROOT/'presentation-validation.json'
    assert json.loads(summary.read_text())['completed']
    assert json.loads(presentation.read_text())['completed']
    job = os.environ['SLURM_JOB_ID']
    live_log = ROOT/f'work/slurm-{job}.log'
    archive = ROOT/'work-artifacts.tar.gz'
    staged = archive.with_name(archive.name+'.tmp')
    receipt = ROOT/'work-archive.json'
    assert not archive.exists() and not staged.exists() and not receipt.exists(), 'Preserve existing archive evidence'
    files, retained = scan(ROOT, PREFIXES, live_log)
    manifest = dict(schema_version=1, original_root=str(ROOT), prefixes=PREFIXES,
        archived=files, retained=retained,
        validation_summary_sha256=sha(summary), presentation_validation_sha256=sha(presentation))
    manifest_bytes = (json.dumps(manifest, indent=2, allow_nan=False)+'\n').encode()
    report = dict(completed=False, job_id=job, started_at_utc=datetime.now(timezone.utc).isoformat(),
                  archiver_sha256=sha(__file__), manifest=manifest, removed_files=[], removed_directories=[])
    try:
        with tarfile.open(staged, 'w:gz', compresslevel=6, dereference=False) as tar:
            for name, record in files.items():
                verify_original(ROOT/name, record)
                tar.add(ROOT/name, arcname=name, recursive=False)
            member = tarfile.TarInfo('MANIFEST.json')
            member.size, member.mode = len(manifest_bytes), 0o644
            tar.addfile(member, io.BytesIO(manifest_bytes))
        with staged.open('rb') as stream:
            os.fsync(stream.fileno())
        with tarfile.open(staged, 'r:gz') as tar:
            members = {member.name: member for member in tar.getmembers()}
            assert set(members) == set(files) | {'MANIFEST.json'}
            assert tar.extractfile(members['MANIFEST.json']).read() == manifest_bytes
            for name, record in files.items():
                member = members[name]
                assert member.mode == record['mode']
                if record['kind'] == 'symlink':
                    assert member.issym() and member.linkname == record['target']
                else:
                    assert member.isfile() and member.size == record['bytes']
                    assert digest(tar.extractfile(member)) == record['sha256']
        staged.replace(archive)
        report.update(archive_sha256=sha(archive), archive_bytes=archive.stat().st_size,
                      archive_verified_before_removal=True)
        write_json(receipt, report)
        # All originals are checked again before the first deletion. Only
        # archive-listed files/links are unlinked; no recursive deletion.
        for name, record in files.items():
            verify_original(ROOT/name, record)
        for name in files:
            verify_original(ROOT/name, files[name])
            (ROOT/name).unlink()
            report['removed_files'].append(name)
        for prefix in PREFIXES:
            for directory, dirs, names in os.walk(ROOT/prefix, topdown=False):
                path = Path(directory)
                if not any(path.iterdir()):
                    path.rmdir()
                    report['removed_directories'].append(str(path.relative_to(ROOT)))
        for name, record in retained.items():
            path = ROOT/name
            if path == live_log:
                continue
            info = path.lstat()
            assert (info.st_dev, info.st_ino, info.st_mtime_ns, info.st_size) == (
                record['device'], record['inode'], record['mtime_ns'], record['bytes'])
        report.update(completed=True, finished_at_utc=datetime.now(timezone.utc).isoformat(),
            archived_file_count=len(files), retained_file_count=len(retained),
            removed_directory_count=len(report['removed_directories']))
    finally:
        write_json(receipt, report)
    print(f'Archive verified: {len(files)} small files consolidated; {len(retained)} data/log files retained')


if __name__ == '__main__':
    main()
