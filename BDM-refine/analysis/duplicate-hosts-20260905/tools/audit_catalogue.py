"""Read-only pilot: micromamba run -n cosemu python3 audit_catalogue.py ..."""
import argparse
import json
import os
from pathlib import Path
import resource
import time

from catalogue_core import audit, read_ascii, sha256, RULE

p = argparse.ArgumentParser()
p.add_argument('source', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--velocities', action='store_true')
a = p.parse_args()
t = time.monotonic()
data, header = read_ascii(a.source)
result, *_ = audit(data, velocities=a.velocities)
result.update(source=str(a.source.resolve()), source_sha256=sha256(a.source),
              rule=RULE, elapsed_seconds=time.monotonic()-t,
              maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              job_id=os.environ.get('SLURM_JOB_ID'))
a.output.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in
                 ['remaining_ambiguous_pairs','mass_bins']}, indent=2), flush=True)
