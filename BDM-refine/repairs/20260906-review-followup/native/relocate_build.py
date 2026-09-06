"""Move a completed build atomically, preserving commands and proving file identity.

Run with micromamba run -n cosemu python3 -B. The destination must not exist.
This handles small build/preflight artifacts; it never copies a build tree.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import stat


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while data := stream.read(8 * 1024**2):
            result.update(data)
    return result.hexdigest()


def manifest(root):
    result = {}
    for path in sorted(root.rglob('*')):
        metadata = path.lstat()
        if stat.S_ISDIR(metadata.st_mode):
            continue
        result[str(path.relative_to(root))] = dict(device=metadata.st_dev, inode=metadata.st_ino,
            bytes=metadata.st_size, mtime_ns=metadata.st_mtime_ns,
            type='symlink' if path.is_symlink() else 'file')
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    receipt = args.receipt.resolve()
    record = json.loads(receipt.read_text())
    assert record['completed'] and 'relocation' not in record
    source = Path(record['frozen_build_path'])
    destination = args.destination.resolve()
    assert source.is_dir() and not destination.exists()
    assert source.stat().st_dev == destination.parent.stat().st_dev, 'Atomic same-filesystem move required'
    previous_receipt_sha = sha(receipt)
    native = Path(record['frozen_source_path'])
    for name, expected in record['source_sha256'].items():
        assert sha(native/name) == expected
    for name, expected in record['common_objects_sha256'].items():
        assert sha(native/name) == expected
    for details in record['variants'].values():
        assert sha(details['binary_path']) == details['binary_sha256']
        assert sha(details['finder_source_path']) == details['source_sha256']
        assert sha(details['uninstrumented_source_path']) == details['uninstrumented_source_sha256']
    # Internal absolute fixture links would break on relocation. Making only
    # those links relative preserves the exact target file identity and bytes.
    links = []
    for path in sorted(source.rglob('*')):
        if not path.is_symlink():
            continue
        target = Path(os.readlink(path))
        if not target.is_absolute() or not target.is_relative_to(source):
            continue
        identity = target.stat()
        relative = os.path.relpath(target, path.parent)
        temporary = path.with_name(path.name + '.relative-link')
        assert not temporary.exists()
        temporary.symlink_to(relative)
        os.replace(temporary, path)
        assert path.stat().st_ino == identity.st_ino and path.stat().st_dev == identity.st_dev
        links.append(dict(path=str(path.relative_to(source)), original_target=str(target), relative_target=relative,
                          target_device=identity.st_dev, target_inode=identity.st_ino))
    before = manifest(source)
    os.rename(source, destination)
    after = manifest(destination)
    assert before == after, 'Atomic relocation changed per-file identity or metadata'
    for link in links:
        target = (destination/link['path']).stat()
        assert (target.st_dev, target.st_ino) == (link['target_device'], link['target_inode'])

    def moved(path):
        return str(destination/Path(path).relative_to(source))

    record['frozen_build_path'] = str(destination)
    record['frozen_source_path'] = moved(record['frozen_source_path'])
    record['preflight_long_density_path'] = moved(record['preflight_long_density_path'])
    record['production_binary']['binary_path'] = moved(record['production_binary']['binary_path'])
    for details in record['variants'].values():
        for key in ['binary_path', 'finder_source_path', 'uninstrumented_source_path']:
            details[key] = moved(details[key])
        details['source_path'] = details['finder_source_path']
        assert sha(details['binary_path']) == details['binary_sha256']
        assert sha(details['source_path']) == details['source_sha256']
    native = Path(record['frozen_source_path'])
    for name, expected in record['source_sha256'].items():
        assert sha(native/name) == expected
    for name, expected in record['common_objects_sha256'].items():
        assert sha(native/name) == expected
    relocation = dict(completed=True, completed_at_utc=datetime.now(timezone.utc).isoformat(),
        method='os.rename; same-filesystem atomic directory move; no build files copied',
        source=str(source), destination=str(destination), previous_build_receipt_sha256=previous_receipt_sha,
        relocator_sha256=sha(__file__), original_commands_and_working_directories_preserved=True,
        internal_fixture_links_made_relative_before_move=links,
        per_file_metadata_before=before, per_file_metadata_after=after,
        source_object_and_executable_hashes_reverified=True)
    relocation_path = receipt.parent/'relocation.json'
    assert not relocation_path.exists()
    relocation_path.write_text(json.dumps(relocation, indent=2) + '\n')
    record['relocation'] = dict(receipt_path=str(destination.parents[1]/'native/relocation.json'),
        receipt_sha256=sha(relocation_path), original_frozen_build_path=str(source),
        original_build_receipt_sha256=previous_receipt_sha)
    temporary = receipt.with_suffix('.json.part')
    temporary.write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    os.replace(temporary, receipt)
    print('Atomically relocated', len(before), 'unchanged build entries to', destination)


if __name__ == '__main__':
    main()
