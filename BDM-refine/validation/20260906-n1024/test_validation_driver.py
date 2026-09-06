"""Adversarial controls for false-success and stale-output validation hazards."""
import argparse
import json
from pathlib import Path
import struct
import tempfile

import numpy as np
import run_validation as v


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--receipt', type=Path, default=v.ROOT/'driver-controls.json')
    args = parser.parse_args()
    passed = []

    def rejected(name, action):
        try:
            action()
        except (AssertionError, RuntimeError, ValueError):
            passed.append(name)
        else:
            raise AssertionError(f'Bad input accepted: {name}')

    with tempfile.TemporaryDirectory(dir=v.WORK, prefix='driver-controls-') as temporary:
        root = Path(temporary)
        with np.load(v.REPAIRS/'membership-n128.npz') as members, \
                np.load(v.REPAIRS/'native-n128-catalogues.npz') as tables:
            reference = json.loads((v.REPAIRS/'native-n128.json').read_text())
            meta = next(row['membership'] for row in reference['records'] if 'membership' in row)
            table = tables['baseline_0_t1']
            raw = bytearray(struct.pack('>qqf', meta['selected'], meta['candidates'], meta['mass_one']))
            starts = []
            for i, candidate in enumerate(members['candidates']):
                starts.append(len(raw))
                a, b = members['offsets'][i:i+2]
                ids = members['ids'][a:b]
                raw.extend(struct.pack('>qq', int(candidate), len(ids)))
                raw.extend(members['properties'][i].astype('>f4').tobytes())
                raw.extend(ids.astype('>i8').tobytes())
            spec = v.case_spec('pilot')
            spec.update(nrow=128, ngrid=256, box_mpc_h=128)
            path = root/'members.bin'
            path.write_bytes(raw)
            result = v.check_memberships(path, spec, table)
            assert result['selected'] == 1215
            passed.append('verified native128 tape agrees with23 printed fields')

            bad = bytearray(raw)
            struct.pack_into('>q', bad, starts[0], 0)
            path.write_bytes(bad)
            rejected('zero candidate ID', lambda: v.check_memberships(path, spec, table))

            bad = bytearray(raw)
            bad[starts[1]:starts[1]+8] = bad[starts[0]:starts[0]+8]
            path.write_bytes(bad)
            rejected('repeated candidate ID', lambda: v.check_memberships(path, spec, table))

            bad = bytearray(raw)
            struct.pack_into('>f', bad, 16, meta['mass_one']*1.01)
            path.write_bytes(bad)
            rejected('wrong snapshot particle mass', lambda: v.check_memberships(path, spec, table))

            path.write_bytes(raw)
            wrong_table = table.copy()
            wrong_table[0, 0] += .1
            rejected('same-size catalogue from another halo set', lambda: v.check_memberships(path, spec, wrong_table))

        log = root/'old.log'
        log.write_text('iVirial = 1\nRext = .15\nMassMin = 2.5e12\nSTOP Too many buffer particles\n')
        rejected('exit-zero STOP with header-only old catalogue',
                 lambda: v.verify_old_completion(log, np.empty((0, 24))))
        rejected('partial nonempty old catalogue before WriteFiles completion',
                 lambda: v.verify_old_completion(log, np.zeros((1, 24))))
        log.write_text('time for WriteFiles  =  1.00 2.00\n')
        v.verify_old_completion(log, np.zeros((1, 24)))
        passed.append('actual legacy completion-marker spacing')

        # Exercise the reuse path without executing any native process.
        stage = root/'reuse'
        stage.mkdir()
        (stage/'BDM.config').write_text(v.CONFIG)
        (stage/'PMcrd.0007.DAT').write_bytes(b'controlled input')
        (stage/'baseline.log').write_text('controlled completed process')
        product = stage/'catalogue.DAT'
        product.write_bytes(b'controlled completed output')
        receipt = dict(binary='PMP2BDM.exe', binary_sha256=v.sha(v.BIN/'PMP2BDM.exe'),
                       threads=1, stdin='7\n', completed=True,
                       log_path=str(stage/'baseline.log'), log_sha256=v.sha(stage/'baseline.log'),
                       inputs_sha256={str(p): v.sha(p) for p in [stage/'BDM.config', stage/'PMcrd.0007.DAT']},
                       outputs_sha256={str(product): v.sha(product)})
        v.write_json(stage/'baseline.json', receipt)
        action = lambda: v.run_native('PMP2BDM.exe', stage, '7\n', 1, 'baseline', [])
        assert action()['completed']
        passed.append('unchanged completed stage reused')
        product.write_bytes(b'changed same-length output')
        rejected('changed completed-stage output', action)
        product.write_bytes(b'controlled completed output')
        (stage/'PMcrd.0007.DAT').write_bytes(b'changed controlled input')
        rejected('changed completed-stage input', action)
        (stage/'PMcrd.0007.DAT').write_bytes(b'controlled input')
        aggregate = root/'aggregate.json'
        v.write_json(aggregate, dict(spec=spec, completed=True, records=[receipt]))
        before = v.sha(aggregate)
        assert v.verify_report(aggregate, spec)
        assert v.sha(aggregate) == before
        passed.append('completed aggregate preserved')
        product.write_bytes(b'changed same-length output')
        rejected('changed output cannot be blessed in a new aggregate', lambda: v.verify_report(aggregate, spec))

        reference = np.ones((3, 24))
        replay = reference.copy()
        replay[[0, 2], 0] += .001
        comparison = v.compare_normal_density_catalogues(reference, replay)
        assert comparison['changed_rows'] == 2 and not comparison['identical']
        assert comparison['bound_mass_counts_velocities_identical']
        assert [row['row_one_based'] for row in comparison['first_changed_rows']] == [1, 3]
        passed.append('density-dependent position changes remain explicit')
        replay[1, 6] += 1
        comparison = v.compare_normal_density_catalogues(reference, replay)
        assert not comparison['bound_mass_counts_velocities_identical']
        passed.append('changed bound mass is not classified as identical')
        comparison = v.compare_normal_density_catalogues(reference, replay[:2])
        assert not comparison['same_shape'] and 'changed_rows' not in comparison
        passed.append('changed selection prevents rowwise residual comparison')

    result = dict(driver_sha256=v.sha(v.__file__), test_sha256=v.sha(__file__),
                  passed=passed, completed=True,
                  native_fixture='Archived repaired128^3 integration data; no large simulation or new native execution')
    v.write_json(args.receipt, result)
    print(len(passed), 'adversarial driver controls passed')


if __name__ == '__main__':
    main()
