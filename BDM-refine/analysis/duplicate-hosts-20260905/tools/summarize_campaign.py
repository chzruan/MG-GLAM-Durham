"""Collect validated per-catalogue receipts without re-reading science arrays."""
import argparse
import csv
import json
from pathlib import Path

import h5py

from catalogue_core import sha256
from clean_catalogue import load_receipt
from run_campaign import output_paths


def summarize(manifest, output_root, report_root, require_complete=False):
    entries=json.loads(manifest.read_text())
    rows=[];mass_rows=[];missing=[];invalid=[]
    for entry in entries:
        output,sidecar=output_paths(entry,output_root)
        if not sidecar.exists():
            missing.append(entry)
            continue
        try:
            r=load_receipt(sidecar)
            if not output.exists() or r['source']!=entry['source']:
                raise ValueError('Missing output or source path mismatch')
            if not r['bitwise_validation']['all_surviving_columns_byte_identical']:
                raise ValueError('Failed bitwise check')
            with h5py.File(sidecar,'r') as f:
                if f['drop_mask'].shape!=(r['source_rows'],):
                    raise ValueError('Mask not aligned to the complete source')
                v=json.loads(f['validation_json'][...].tobytes())
            if v['remaining_strict_pairs']!=0:
                raise ValueError('Remaining strict duplicate edges')
            row={key:entry[key] for key in ['gravity','imodel','ibox','redshift','snapnum']}
            row.update(source_rows=r['source_rows'],removed_rows=r['removed_rows'],
                removed_fraction=r['removed_fraction'],selected_rows=r['selected_rows'],
                selected_removed_rows=r['selected_removed_rows'],
                selected_only_strict_removed=r['selected_only_strict_removed'],
                selected_fraction=r['selected_removed_rows']/r['selected_rows'] if r['selected_rows'] else None,
                exact_removed_rows=r['exact_removed_rows'],
                remaining_same_bound_pairs_below_0p2=r['remaining_same_bound_pairs_below_0p2'],
                remaining_strict_pairs=v['remaining_strict_pairs'],
                source_sha256=r['source_sha256'],output_sha256=r['output_sha256'],
                sidecar_sha256=sha256(sidecar),tool_git_commit=r['tool_git_commit'],
                source=r['source'],output=str(output.resolve()),sidecar=str(sidecar.resolve()),
                elapsed_seconds=r['elapsed_seconds'],job_id=r['job_id'])
            rows.append(row)
            for m in v['mass_bins']:
                mass_rows.append({**{key:entry[key] for key in ['gravity','imodel','ibox','redshift']},**m})
        except Exception as error:
            invalid.append(dict(source=entry['source'],error=str(error)))
    report_root.mkdir(parents=True,exist_ok=True)
    for name,records in [('catalogue_validation.csv',rows),('mass_bin_validation.csv',mass_rows)]:
        if records:
            with open(report_root/name,'w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n')
                w.writeheader();w.writerows(records)
    totals={key:sum(r[key] for r in rows) for key in ['source_rows','removed_rows',
             'selected_rows','selected_removed_rows','remaining_same_bound_pairs_below_0p2']}
    groups=[]
    for gravity,model in sorted(set((r['gravity'],r['imodel']) for r in rows)):
        rs=[r for r in rows if (r['gravity'],r['imodel'])==(gravity,model)]
        count=sum(r['selected_rows'] for r in rs)
        removed=sum(r['selected_removed_rows'] for r in rs)
        groups.append(dict(gravity=gravity,imodel=model,catalogues=len(rs),
                           selected_rows=count,selected_removed_rows=removed,
                           selected_fraction=removed/count if count else None))
    result=dict(expected_catalogues=len(entries),complete_catalogues=len(rows),
                totals=totals,per_model=groups,missing=missing,invalid=invalid)
    (report_root/'campaign_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['missing','per_model']},indent=2))
    if require_complete and (missing or invalid):
        raise RuntimeError('Campaign incomplete; see campaign_summary.json')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest',type=Path)
    p.add_argument('output_root',type=Path)
    p.add_argument('report_root',type=Path)
    p.add_argument('--require-complete',action='store_true')
    a=p.parse_args()
    summarize(a.manifest,a.output_root,a.report_root,a.require_complete)
