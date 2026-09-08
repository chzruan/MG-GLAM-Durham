"""Convergence from retained catalogues and original-particle membership.

Expensive membership matching is separate from plotting. Use a small shared
Slurm allocation for the full campaign; memory is bounded by halo counts and
the largest two member sets, not by the full particle cube.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

from common import MATRIX, REPO, ROOT, WORK, now, sha, verify_manifest, write_json

PAIRS=[('A','C','particle'),('C','E','particle'),('A','E','particle'),
       ('B','C','force'),('C','D','force'),('E','F','force'),('F','T','time')]
EDGES=np.arange(12.25,15.76,.25)


def finite_json(value):
    if isinstance(value,np.ndarray):return finite_json(value.tolist())
    if isinstance(value,dict):return {k:finite_json(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [finite_json(v) for v in value]
    if isinstance(value,(np.integer,)):return int(value)
    if isinstance(value,(float,np.floating)):return float(value) if np.isfinite(value) else None
    if isinstance(value,np.bool_):return bool(value)
    return value


def load(name,z,verify=True):
    path=ROOT/f'v3-validation-{name}-z{z}-ng2048.json'
    if not path.exists():return None
    report=json.loads(path.read_text());assert report['completed'] and 'v3' in report['variants']
    if verify:verify_manifest(report['science_files'])
    evidence=report['variants']['v3']['membership']
    raw=REPO/evidence['retained_raw'];index=raw.with_suffix('.index.npz')
    assert sha(index)==evidence['index_sha256']
    with np.load(index) as data:contents={k:data[k] for k in data.files}
    return dict(name=name,z=z,report=report,raw=raw,receipt=path,receipt_sha256=sha(path),**contents)


def project_ids(ids,nfine,ncoarse):
    """Restrict fine IDs to the shared initial lattice, preserving original IDs."""
    if nfine==ncoarse:return ids
    assert nfine%ncoarse==0
    ratio=nfine//ncoarse;zero=ids-1
    x=zero%nfine;y=(zero//nfine)%nfine;z=zero//nfine**2
    keep=(x%ratio==0)&(y%ratio==0)&(z%ratio==0)
    return 1+x[keep]//ratio+ncoarse*(y[keep]//ratio)+ncoarse**2*(z[keep]//ratio)


def read_ids(stream,sample,index):
    stream.seek(int(sample['offsets'][index]))
    return np.frombuffer(stream.read(int(sample['counts'][index])*8),dtype='>i8').astype(np.int64)


def match(left,right):
    nleft=left['report']['spec']['nrow'];nright=right['report']['spec']['nrow']
    assert nleft<=nright and nright%nleft==0
    a=left['properties'];b=right['properties'];box=left['report']['spec']['box_mpc_h']
    assert box==right['report']['spec']['box_mpc_h']
    tree=cKDTree(b[:,:3]%box,boxsize=box)
    radius=np.maximum(a[:,8],np.max(b[:,8],initial=0.))+box/2048
    neighbours=tree.query_ball_point(a[:,:3]%box,radius,workers=1)
    best_left={};best_right={};examined=0;spatial_candidates=0
    with left['raw'].open('rb') as lhs,right['raw'].open('rb') as rhs:
        for i,candidates in enumerate(neighbours):
            ids_left=None
            for j in candidates:
                delta=a[i,:3]-b[j,:3];delta-=box*np.rint(delta/box)
                distance=float(np.linalg.norm(delta))
                if distance>max(a[i,8],b[j,8])+box/2048:continue
                spatial_candidates+=1
                if ids_left is None:ids_left=read_ids(lhs,left,i)
                ids_right=project_ids(read_ids(rhs,right,j),nright,nleft)
                examined+=1
                if not len(ids_right):continue
                overlap=len(np.intersect1d(ids_left,ids_right,assume_unique=True))
                fa=overlap/len(ids_left);fb=overlap/len(ids_right)
                if min(fa,fb)<.5:continue
                score=fa*fb
                record=(score,-distance,-j,j,overlap,fa,fb,distance)
                if i not in best_left or record[:3]>best_left[i][:3]:best_left[i]=record
                reverse=(score,-distance,-i,i)
                if j not in best_right or reverse[:3]>best_right[j][:3]:best_right[j]=reverse
    rows=[]
    for i,item in best_left.items():
        j=item[3]
        if best_right[j][3]==i:rows.append((i,j,item[0],item[4],item[5],item[6],item[7]))
    pairs=np.asarray(rows,dtype=np.float64).reshape(-1,7)
    return pairs,dict(method='mutual best shared-lattice membership overlap',
                      minimum_fraction_each=.5,spatial_radius='max(R_left,R_right)+one common analysis cell',
                      physical_cell_mpc_h=box/2048,spatial_candidates=spatial_candidates,
                      member_sets_examined=examined,matches=len(pairs),
                      coarse_rows=len(a),reference_rows=len(b),
                      identity='Original PM IDs mapped to nested initial lattice; no present-day nearest-particle IDs')


def abundance(sample):
    p=sample['properties'];box=sample['report']['spec']['box_mpc_h']
    octant=np.floor((p[:,:3]%box)/(box/2)).astype(int)
    octant=octant[:,0]+2*octant[:,1]+4*octant[:,2]
    counts=np.histogram(np.log10(p[:,6]),EDGES)[0]
    octants=np.array([np.histogram(np.log10(p[octant==j,6]),EDGES)[0] for j in range(8)])
    return dict(counts=counts,octant_counts=octants,dn_dlog10m=counts/(box**3*np.diff(EDGES)),
                mass_one=sample['report']['variants']['v3']['membership']['mass_one'],rows=len(p))


def quantiles(x,y,mask):
    result=[]
    for lo,hi in zip(EDGES[:-1],EDGES[1:]):
        selected=mask&(x>=lo)&(x<hi)&np.isfinite(y)
        vals=y[selected]
        result.append(dict(count=len(vals),q16_median_q84=np.quantile(vals,[.16,.5,.84]) if len(vals) else [None]*3))
    return result


def compare(left,right,pairs,description,floor=300):
    aa=abundance(left);bb=abundance(right)
    ratio=np.divide(aa['counts'],bb['counts'],out=np.full(len(EDGES)-1,np.nan),where=bb['counts']>0)
    leave_a=aa['counts']-aa['octant_counts'];leave_b=bb['counts']-bb['octant_counts']
    ratios=np.divide(leave_a,leave_b,out=np.full_like(leave_a,np.nan,dtype=float),where=leave_b>0)
    jk=np.sqrt(7/8*np.sum((ratios-np.mean(ratios,axis=0))**2,axis=0))
    ii=pairs[:,0].astype(int);jj=pairs[:,1].astype(int)
    a=left['properties'][ii];b=right['properties'][jj]
    resolved=(left['counts'][ii]>=floor)&(right['counts'][jj]>=floor)
    x=np.log10(b[:,6])
    stats={}
    for label,column in [('bound_mass',6),('so_mass',7),('radius',8),('vmax',11),('axis_ba',16),('axis_ca',17)]:
        value=100*np.divide(a[:,column]-b[:,column],b[:,column],out=np.full(len(a),np.nan),where=b[:,column]>0)
        stats[label]=quantiles(x,value,resolved)
    stats['bulk_velocity_km_s']=quantiles(x,np.linalg.norm(a[:,3:6]-b[:,3:6],axis=1),resolved)
    stats['centre_distance_mpc_h']=quantiles(x,pairs[:,6],resolved)
    mass_floor=max(2.5e12,floor*max(aa['mass_one'],bb['mass_one']))
    eligible=right['properties'][:,6]>=mass_floor
    reference_counts=np.histogram(np.log10(right['properties'][eligible,6]),EDGES)[0]
    selected_counts=np.histogram(x[resolved&eligible[jj]],EDGES)[0]
    completeness=np.divide(selected_counts,reference_counts,out=np.full(len(EDGES)-1,np.nan),where=reference_counts>0)
    left_eligible=left['counts']>=floor
    left_counts=np.histogram(np.log10(left['properties'][left_eligible,6]),EDGES)[0]
    matched_left=np.histogram(np.log10(a[resolved,6]),EDGES)[0]
    purity=np.divide(matched_left,left_counts,out=np.full(len(EDGES)-1,np.nan),where=left_counts>0)
    return dict(kind=description,coarse=left['name'],reference=right['name'],redshift=left['z'],
                particle_floor=floor,mass_floor_msun_h=mass_floor,
                left_counts=aa['counts'],right_counts=bb['counts'],abundance_ratio=ratio,
                abundance_ratio_jackknife8_sigma=jk,matched_statistics=stats,
                reference_eligible_counts=reference_counts,reference_completeness=completeness,
                left_eligible_counts=left_counts,left_matched_fraction=purity,
                uncertainty='Eight spatial delete-one octants, paired ratio; correlated and limited-volume estimate',
                shifts='100*(coarse/reference-1); velocity and centre distances are absolute',
                catalogue_edge='Abundance uses raw audited membership masses, avoiding ASCII rounding at bin edges')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--allow-partial',action='store_true')
    args=parser.parse_args()
    destination=WORK/'analysis';destination.mkdir(parents=True,exist_ok=True)
    samples={};missing=[]
    for z in [2,1,0]:
        for name in MATRIX:
            sample=load(name,z)
            if sample is None:missing.append(f'{name}/z{z}')
            else:samples[name,z]=sample
    if missing and not args.allow_partial:raise RuntimeError('Missing verified catalogues: '+', '.join(missing))
    result=dict(completed=not missing,created_at_utc=now(),missing=missing,log10_mass_edges=EDGES,
                particle_floors=[100,300,1000],abundances={},comparisons=[],
                caveat='Finest is a comparison reference, not ground truth; convergence range is empirical.')
    for (name,z),sample in samples.items():result['abundances'][f'{name}/z{z}']=abundance(sample)
    for z in [2,1,0]:
        for lhs,rhs,kind in PAIRS:
            if (lhs,z) not in samples or (rhs,z) not in samples:continue
            left=samples[lhs,z];right=samples[rhs,z]
            path=destination/f'matches-{lhs}-{rhs}-z{z}.npz';receipt=path.with_suffix('.json')
            identity=dict(left_membership=left['report']['variants']['v3']['membership']['raw_sha256'],
                          right_membership=right['report']['variants']['v3']['membership']['raw_sha256'],
                          left_nrow=left['report']['spec']['nrow'],right_nrow=right['report']['spec']['nrow'],
                          analysis_sha256=sha(sys.argv[0]))
            if receipt.exists():
                saved=json.loads(receipt.read_text());assert saved['identity']==identity
                assert sha(path)==saved['sha256']
                with np.load(path) as arrays:pairs=arrays['pairs']
                description=saved['matching']
            else:
                pairs,description=match(left,right)
                np.savez_compressed(path,pairs=pairs)
                write_json(receipt,dict(identity=identity,matching=description,sha256=sha(path)))
            for floor in [100,300,1000]:
                row=compare(left,right,pairs,kind,floor);row['matching']=description
                result['comparisons'].append(row)
            print(lhs,rhs,'z=',z,description['matches'],'mutual membership matches',flush=True)
    write_json(ROOT/'convergence.json',finite_json(result))
    print('Convergence measurements written; completed=',not missing)


if __name__=='__main__':main()
