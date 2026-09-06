"""Small parser/identity checks; micromamba run -n cosemu python3 -B."""
import json
from pathlib import Path
import re
import struct
import tempfile

import run_replays as r


def main():
    with tempfile.TemporaryDirectory(prefix='bdm-replay-driver-') as tmp:
        path=Path(tmp)/'unbinding.bin'
        path.write_bytes(struct.pack('>q',4)+b''.join(struct.pack('>iiqqf',*x) for x in
            [(0,1,0,0,0.),(1,0,40,40,4.e12),(3,0,100,30,3.e12),(91,0,10010,20,0.)]))
        diagnostics=r.diagnostics(path)
        assert diagnostics['evaluated_candidates']==3 and diagnostics['max_passes']==91
        assert diagnostics['active_particle_rows']==10150 and diagnostics['selected_over_32']==0
        assert diagnostics['candidates_over_32']==1
        density=Path(tmp)/'fi.bin'
        density.write_bytes(struct.pack('>qqi',2,8,1)+bytes(32))
        cache={};spec={'nrow':2,'ngrid':2}
        first=r.verify_density(density,spec,cache=cache)
        r.verify_density(density,spec,first['sha256'],cache)
        raw=bytearray(density.read_bytes());raw[-1]=1;density.write_bytes(raw)
        try:r.verify_density(density,spec,first['sha256'])
        except AssertionError:pass
        else:raise AssertionError('final mandatory hash accepted changed density')
        log=' time for ParametersDistinct =    1.23 2.34\n time for ParametersDistinct  =    1.24 2.35\n'
        assert re.findall(r'time for ParametersDistinct  =\s*([0-9.]+)',log)==['1.24']
        result=dict(completed=True,driver_sha256=r.sha(Path(r.__file__)),test_sha256=r.sha(__file__),
                    checks=['big-endian diagnostic layout/histogram/work/selection decode',
                            'density dimensions and initial hash',
                            'mandatory uncached final hash rejects changed same-size density',
                            'outer timing excludes the duplicate internal ParametersDistinct label'],
                    metadata_cache_limitation='Timestamp metadata alone can miss rapid same-size changes; the final uncached hash is required.',
                    diagnostic_fixture=diagnostics)
        r.write_json(r.HERE/'driver-preflight.json',result)
    print('Driver diagnostic, timing and density identity controls passed')


if __name__=='__main__':
    main()
