"""Tiny synthetic adversarial controls; no native particle or mesh allocations."""
import copy
import struct
import tempfile
from pathlib import Path

import run_control as c


def rejected(action):
    try:
        action()
    except (AssertionError, KeyError):
        return
    raise AssertionError('Invalid experiment accepted')


def main():
    passed = []
    with tempfile.TemporaryDirectory(prefix='bdm-main-control-tests-') as directory:
        path = Path(directory) / 'header.DAT'
        raw = bytearray(537)
        raw[:4] = raw[-4:] = struct.pack('>i', 529)
        struct.pack_into('>f', raw, 49, 1.)
        struct.pack_into('>i', raw, 65, 173)
        struct.pack_into('>ii', raw, 97, 1024, 2048)
        struct.pack_into('>q', raw, 125, 1024**3)
        path.write_bytes(raw)
        header = c.read_header(path)
        assert header == dict(scale_factor=1., step=173, nrow=1024, ngrid=2048, particles=1024**3)
        passed.append('big-endian header and actual non-157 step')
        struct.pack_into('>q', raw, 125, 2**33+1)
        path.write_bytes(raw)
        assert c.read_header(path)['particles'] == 2**33+1
        passed.append('int64 particle header above int32 range')
        path.write_bytes(raw[:-1])
        rejected(lambda:c.read_header(path))
        passed.append('truncated header rejected')
    spec = dict(nrow=1024, ngrid=2048, box_mpc_h=512., epochs=[2, 1, 0], sources_sha256={'finder':'abc'},
                binaries_sha256={'PMP2main.exe':'def'}, halo_config='config', realization=1, gravity='GR')
    plan = dict(parent_job_id='11948372', source_sha256={'finder':'abc'},
                production_binaries_sha256={'PMP2main.exe':'def'}, config='config')
    snapshot = dict(z=0, header=header, files_sha256={name:'sha' for name in ['PMcrd.0173.DAT', 'PMcrs0.0173.DAT', 'PMcrs1.0173.DAT']})
    simulation = dict(completed=True, job_id='11948372', spec=spec, snapshots=[snapshot],
                      records=[dict(binary='PMP2main.exe', completed=True, returncode=0, binary_sha256='def', threads=64)])
    assert c.select_snapshot(simulation, plan) == snapshot
    passed.append('actual final step selected without hardcoded 157')
    for name, change in [('incomplete simulation', lambda x:x.update(completed=False)),
                         ('wrong source', lambda x:x['spec'].update(sources_sha256={'finder':'changed'})),
                         ('missing z0', lambda x:x.update(snapshots=[])),
                         ('duplicate z0', lambda x:x['snapshots'].append(copy.deepcopy(snapshot))),
                         ('wrong snapshot filename step', lambda x:x['snapshots'][0]['header'].update(step=157)),
                         ('wrong particle count', lambda x:x['snapshots'][0]['header'].update(particles=512**3)),
                         ('failed main stage', lambda x:x['records'][0].update(returncode=1))]:
        invalid = copy.deepcopy(simulation)
        change(invalid)
        rejected(lambda:c.select_snapshot(invalid, plan))
        passed.append(name+' rejected')
    text = ''.join(f'PROBE FIXED DENSITY VERIFIED {tag} bit_differences=0\n'
                   f'PROBE FINDER COMPLETE {tag} seconds=1.0\n'
                   f'PROBE PARTICLE RESTORATION VERIFIED {tag}\n' for tag in c.TAGS)
    text += 'BDM FIXED DENSITY THREAD CONTROL COMPLETE\n'
    assert c.restoration_markers(text)[0]['particle_restoration_verified'] == 6
    passed.append('six ordered density and particle restoration markers')
    rejected(lambda:c.restoration_markers(text.replace('bit_differences=0', 'bit_differences=1', 1)))
    passed.append('nonzero restored density difference rejected')
    rejected(lambda:c.restoration_markers(text.replace('PROBE PARTICLE RESTORATION VERIFIED d64-t32-p1', 'omitted')))
    passed.append('missing particle restoration rejected')
    rejected(lambda:c.restoration_markers(text.replace('BDM FIXED DENSITY THREAD CONTROL COMPLETE', '')))
    passed.append('exit-zero without native completion rejected')
    c.require_fixed_groups([dict(byte_identical=True)]*4)
    rejected(lambda:c.require_fixed_groups([dict(byte_identical=True)]*3+[dict(byte_identical=False)]))
    passed.append('one unequal fixed-FI catalogue makes experiment fail')
    c.write_json(c.ROOT / 'driver-controls.json', dict(completed=True, tests=passed,
                 driver_sha256=c.sha(c.ROOT / 'run_control.py'), test_sha256=c.sha(__file__)))
    print(len(passed), 'main threading controls passed')


if __name__ == '__main__':
    main()
