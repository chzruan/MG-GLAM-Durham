"""Capped shared-queue workers; one immutable catalogue and sidecar per input."""
import argparse
import json
import os
from pathlib import Path
import traceback

from catalogue_core import sha256
from clean_catalogue import clean, load_receipt


def output_paths(entry, output_root):
    relative = Path(entry['source']).relative_to('/cosma8/data/dp203/dc-ruan1')
    output = output_root/relative
    return output, output.with_suffix('.cleaning.hdf5')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest',type=Path)
    p.add_argument('output_root',type=Path)
    p.add_argument('--workers',type=int,default=8)
    p.add_argument('--worker',type=int,default=int(os.environ.get('SLURM_ARRAY_TASK_ID','0')))
    a=p.parse_args()
    if a.workers <= 0 or not 0 <= a.worker < a.workers:
        p.error('Invalid worker assignment')
    entries=json.loads(a.manifest.read_text())
    failed=[]
    for index in range(a.worker,len(entries),a.workers):
        entry=entries[index]
        output,sidecar=output_paths(entry,a.output_root)
        try:
            if sidecar.exists():
                receipt=load_receipt(sidecar)
                if (receipt['source']!=entry['source'] or not output.exists() or
                    sha256(output)!=receipt['output_sha256'] or
                    sha256(entry['source'])!=receipt['source_sha256']):
                    raise ValueError('Existing receipt/output/source mismatch')
                print(json.dumps(dict(index=index,status='verified-existing',source=entry['source'])),flush=True)
                continue
            clean(entry['source'],output,sidecar,
                  {key:entry[key] for key in ['gravity','imodel','ibox','redshift','snapnum']},
                  velocities=(entry['imodel']==0 and entry['ibox']==1))
        except Exception:
            traceback.print_exc()
            failed.append(dict(index=index,source=entry['source']))
    if failed:
        print(json.dumps(dict(failed=failed)),flush=True)
        raise SystemExit(1)
    print(json.dumps(dict(worker=a.worker,status='complete')),flush=True)


if __name__=='__main__':
    main()
