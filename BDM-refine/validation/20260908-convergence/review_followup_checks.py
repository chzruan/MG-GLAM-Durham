"""Small controls for the review response; run through micromamba cosemu."""
import ast
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import types
import zipfile
import zipimport

os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy.spatial import cKDTree

import assess
import ic_validate
from common import REPO, ROOT, WORK, now, sha, write_json


def main():
    data=json.loads((ROOT/'convergence.json').read_text())
    old=json.loads(subprocess.check_output(['git','show',
        'cb87021:BDM-refine/validation/20260908-convergence/convergence-assessment.json'],cwd=REPO))
    current=json.loads((ROOT/'convergence-assessment.json').read_text())
    assert current['comparisons']==old['comparisons'] and current['criteria']==old['criteria']
    assert len(current['comparisons'])==63
    assert current['assess_source_sha256']==sha(ROOT/'assess.py')
    assert sha(ROOT/'convergence.json')=='31d96e95f8a6df13ecfdefd11df88b3c811cc73eb15cb21bc4677f7a5221fa3d'
    with tempfile.TemporaryDirectory(prefix='review-followup-',dir=WORK) as temp:
        frozen=Path(temp)/'sources.pyz'
        with zipfile.ZipFile(frozen,'w') as archive:
            for name in ['common.py','assess.py','ic/configure.py']:
                archive.write(ROOT/name,name)
            archive.write(ROOT/'ic_validate.py','__main__.py')
        assert ic_validate.source_identity(frozen)==ic_validate.source_identity()
        loader=zipimport.zipimporter(str(frozen));module=types.ModuleType('frozen_assess')
        module.__loader__=loader;module.__file__=str(frozen/'assess.py')
        exec(loader.get_code('assess'),module.__dict__)
        assert module.source_sha256()==assess.source_sha256()
        review=REPO/'BDM-refine/analysis/review-20260909-convergence/review_checks.py'
        spec=importlib.util.spec_from_file_location('independent_review',review)
        reviewer=importlib.util.module_from_spec(spec);spec.loader.exec_module(reviewer)
        reviewer.RESULTS=Path(temp)/'independent.json'
        reviewer.stage_assess()
        recomputed=json.loads(reviewer.RESULTS.read_text())['independent_assessment']
        assert recomputed['reproduced']==63 and not recomputed['mismatches']
        independent=dict(reproduced=recomputed['reproduced'],mismatches=recomputed['mismatches'],
                         reviewer_source_sha256=sha(review))

    # Execute the frozen checker's actual host block, not a rewritten predicate.
    checker=REPO/'BDM-refine/validation/20260906-n1024/run_validation.py'
    tree=ast.parse(checker.read_text());block=None
    for node in tree.body:
        if not isinstance(node,ast.FunctionDef):continue
        start=next((i for i,s in enumerate(node.body) if isinstance(s,ast.Assign)
                    and any(isinstance(t,ast.Name) and t.id=='box' for t in s.targets)),None)
        stop=next((i for i,s in enumerate(node.body) if isinstance(s,ast.Assert)
                   and isinstance(s.test,ast.UnaryOp) and isinstance(s.test.operand,ast.Name)
                   and s.test.operand.id=='violations'),None)
        if start is not None and stop is not None:
            block=ast.Module(body=copy.deepcopy(node.body[start:stop+1]),type_ignores=[])
            break
    assert block is not None
    program=compile(ast.fix_missing_locations(block),str(checker),'exec')
    cases=[('separated',[0.,2.],[2.,1.],[.2,.2],False),
           ('lower_priority_inside_higher_aperture',[0.,.5],[2.,1.],[1.,.1],True),
           ('allowed_reverse_orientation',[0.,.5],[2.,1.],[.1,1.],False),
           ('stable_index_breaks_mass_tie',[0.,.5],[2.,2.],[1.,.1],True),
           ('periodic_violation',[.05,9.95],[2.,1.],[.2,.1],True),
           ('exact_aperture_boundary',[0.,1.],[2.,1.],[1.,.1],False)]
    host=[]
    for name,x,mass,radius,reject in cases:
        props=np.zeros((2,21));props[:,0]=x;props[:,6]=mass;props[:,8]=radius
        scope=dict(np=np,cKDTree=cKDTree,spec={'box_mpc_h':10.},selected=2,
                   props=props,index=np.array([1,2]))
        failed=False
        try:exec(program,scope)
        except AssertionError:failed=True
        assert failed==reject,name
        host.append(dict(case=name,rejected=failed,examined=scope['examined'],violations=scope['violations']))

    octants=np.array([10,9,7,5,2,6,4,4]);counts=int(octants.sum())
    ratios=(counts-octants)/(counts-octants)
    paired_sigma=float(np.sqrt(7/8*np.sum((ratios-ratios.mean())**2)))
    assert paired_sigma==0 and counts==47
    count_scale=math.sqrt(2/counts)
    diagnostics=current['supplemental_diagnostics']
    for row in diagnostics:
        if (row['coarse'],row['reference']) in [('E','F'),('F','T')]:
            assert row['minimum_particles_at_common_mass_floor']=={'coarse':1868,'reference':1868}
            assert row['minimum_particles_at_first_whole_bin']=={'coarse':2362,'reference':2362}
    review=json.loads((REPO/'BDM-refine/analysis/review-20260909-convergence/review_results.json').read_text())
    for row in diagnostics:
        key=f"{row['coarse']}/{row['reference']} z{row['redshift']}"
        if row['nominal_particle_floor']!=300 or key not in review['unscreened_properties_inside_passing_intervals']:continue
        for axis in ['axis_ba','axis_ca']:
            value=row['unscreened_shapes_inside_passing_intervals'][axis]['maximum_absolute_bin_median_shift_percent']
            assert round(value,2)==review['unscreened_properties_inside_passing_intervals'][key][f'max_abs_median_{axis}_pct']
    legacy_with_membership=[]
    for path in ROOT.glob('replay-*-ng2048.json'):
        record=json.loads(path.read_text())
        if 'membership' in record['variants']['legacy']:legacy_with_membership.append(path.name)
    assert not legacy_with_membership
    write_json(ROOT/'review-followup-controls.json',dict(completed=True,completed_at_utc=now(),
        source_sha256=sha(__file__),assessment_source_sha256=sha(ROOT/'assess.py'),
        numerical_input_unchanged=True,all_63_original_decisions_and_criteria_unchanged=True,
        independent_reassessment=independent,plain_and_zipapp_source_hashes_match=True,
        host_checker_sha256=sha(checker),host_block_controls=host,
        paired_identical_counts_example=dict(count_each=counts,paired_octant_sigma=paired_sigma,
            hypothetical_independent_count_relative_scale=count_scale,
            meaning='Independent-count scale omits paired covariance; not the uncertainty of this paired ratio.'),
        effective_particle_cuts_checked=True,shape_summaries_match_independent_review=True,
        legacy_membership_receipts_available=legacy_with_membership))
    print('Review follow-up controls passed; numerical measurements and acceptance decisions unchanged.')


if __name__=='__main__':main()
