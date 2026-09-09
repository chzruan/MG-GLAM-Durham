"""Describe measured convergence, separating count limits from chosen tolerances."""
import json
from pathlib import Path

import numpy as np

from common import ROOT, now, read_json_snapshot, sha, verify_json_snapshot, write_json

CRITERIA=dict(minimum_bin_count=30,abundance_difference_percent=5.,
              maximum_abundance_jackknife_sigma_percent=5.,median_bound_mass_difference_percent=5.,
              median_resolved_vmax_difference_percent=2.,reference_completeness=.9)


def finite(value):return value is not None and np.isfinite(value)


def intervals(mask,edges):
    result=[];start=None
    for i,selected in enumerate([*mask,False]):
        if selected and start is None:start=i
        if not selected and start is not None:
            result.append([float(edges[start]),float(edges[i])]);start=None
    return result


def evaluate(row,edges):
    result=[]
    for i,valid in enumerate(row['valid_mass_bins']):
        mass=row['matched_statistics']['bound_mass'][i]
        vmax=row['matched_statistics']['vmax'][i]
        counts=[row['left_counts'][i],row['right_counts'][i],row['reference_eligible_counts'][i],
                mass['count'],vmax['count']]
        enough=bool(valid and min(counts)>=CRITERIA['minimum_bin_count'])
        ratio=row['abundance_ratio'][i];sigma=row['abundance_ratio_jackknife8_sigma'][i]
        median_mass=mass['q16_median_q84'][1];median_vmax=vmax['q16_median_q84'][1]
        completeness=row['reference_completeness'][i]
        enough=enough and all(map(finite,[ratio,sigma,median_mass,median_vmax,completeness]))
        checks={} if not enough else dict(
            abundance=abs(100*(ratio-1))<=CRITERIA['abundance_difference_percent'],
            abundance_uncertainty=100*sigma<=CRITERIA['maximum_abundance_jackknife_sigma_percent'],
            bound_mass=abs(median_mass)<=CRITERIA['median_bound_mass_difference_percent'],
            resolved_vmax=abs(median_vmax)<=CRITERIA['median_resolved_vmax_difference_percent'],
            completeness=completeness>=CRITERIA['reference_completeness'])
        result.append(dict(log10_mass_interval=[edges[i],edges[i+1]],eligible=bool(enough),checks=checks,
                           meets_working_criteria=bool(enough and all(checks.values()))))
    selected=[r['meets_working_criteria'] for r in result]
    return dict(coarse=row['coarse'],reference=row['reference'],kind=row['kind'],redshift=row['redshift'],
                particle_floor=row['particle_floor'],mass_floor_msun_h=row['mass_floor_msun_h'],
                matching=row['matching'],eligible_bins=sum(r['eligible'] for r in result),
                meets_working_criteria_log10_mass_intervals=intervals(selected,edges),bins=result)


def assess(path):
    data,input_sha=read_json_snapshot(path);edges=data['log10_mass_edges']
    comparisons=[evaluate(row,edges) for row in data['comparisons']]
    overview=[]
    for left,right in [('C','E'),('E','F'),('F','T')]:
        rows=[row for row in data['comparisons'] if (row['coarse'],row['reference'])==(left,right)
              and row['redshift']==0 and row['particle_floor']==300]
        if not rows:continue
        row=rows[0];verdict=evaluate(row,edges)
        eligible=[i for i,b in enumerate(verdict['bins']) if b['eligible']]
        if not eligible:continue
        first=eligible[0]
        overview.append(dict(pair=left+'/'+right,particle_floor=300,
            log10_mass_intervals=verdict['meets_working_criteria_log10_mass_intervals'],
            maximum_absolute_eligible_bin_shifts_percent=dict(
                abundance=max(abs(100*(row['abundance_ratio'][i]-1)) for i in eligible),
                median_bound_mass=max(abs(row['matched_statistics']['bound_mass'][i]['q16_median_q84'][1]) for i in eligible),
                median_vmax=max(abs(row['matched_statistics']['vmax'][i]['q16_median_q84'][1]) for i in eligible)),
            lowest_eligible_bin=dict(log10_mass_interval=[edges[first],edges[first+1]],
                abundance_shift_percent=100*(row['abundance_ratio'][first]-1),
                median_bound_mass_shift_percent=row['matched_statistics']['bound_mass'][first]['q16_median_q84'][1],
                median_vmax_shift_percent=row['matched_statistics']['vmax'][first]['q16_median_q84'][1])))
    result=dict(completed=bool(data['completed']),created_at_utc=now(),input=str(path),input_sha256=input_sha,
                criteria=CRITERIA,comparisons=comparisons,z0_overview=overview,
                initial_conditions=dict(generator='PMP2start.matched.exe, campaign variant of native GLAM',
                    lpt_order=1,z_init=100,
                    scope='Existing first-order campaign; not a convergence test of the default 2LPTIC workflow.'),
                interpretation='Descriptive working tolerances, not a formal confidence interval or proof of absolute accuracy. '
                    'The highest resolution is a reference. One realization and eight spatial octants cannot establish volume convergence.')
    verify_json_snapshot(path,input_sha)
    write_json(ROOT/'convergence-assessment.json',result)
    lines=['# BDM convergence measurements','',
        '**Complete seven-run measurement.**' if data['completed'] else '**Partial campaign: conclusions await the missing runs.**','',
        'All primary catalogues use the same 2048^3 analysis mesh. The simulation box is 256 Mpc/h; outputs are z=2,1,0. '
        'Particle resolution, evolved force resolution and timestep size are compared separately.','',
        '**IC scope:** all seven runs use native GLAM first-order (Zel\'dovich) initial conditions at z_init=100. '
        'This completed suite is an existing-run exception to the default 2LPTIC workflow for new simulations; '
        'it does not establish convergence with 2LPTIC initial conditions.','',
        'Bound mass is the original-member count times the stored particle mass, evaluated in float64. '
        'Reported aperture mass and radius include the empirical Rext expansion. Positive Vmax values are compared only when both haloes resolve them; '
        'unresolved values are counted separately. Halo matches require mutual best shared-lattice membership overlap, with at least 50% in each object.','',
        'The table applies a descriptive working screen: at least 30 objects in each required bin; abundance within 5% with '
        'jackknife sigma at most 5%; median bound mass within 5%; median resolved Vmax within 2%; reference completeness at least 90%. '
        'These choices are adjustable. They do not certify those absolute accuracies. '
        'Scatter, individual criteria and 100/300/1000-particle floors are retained in the JSON results.','',
        '| Comparison | z | Particle floor | Eligible bins | log10 mass intervals meeting working criteria |',
        '|---|---:|---:|---:|---|']
    for row in comparisons:
        if row['particle_floor']!=300:continue
        ranges='; '.join(f'{lo:.2f}–{hi:.2f}' for lo,hi in row['meets_working_criteria_log10_mass_intervals']) or 'None'
        lines.append(f"| {row['coarse']}/{row['reference']} ({row['kind']}) | {row['redshift']} | 300 | {row['eligible_bins']} | {ranges} |")
    lines+=['','## Membership checks','',
        '| Catalogue | Published haloes | Exact duplicate sets | Repeated particle IDs | Mass/count mismatches | Host violations | Higher-priority neighbours examined |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for key,source in data['inputs'].items():
        m=source['membership_checks']
        lines.append('| '+key+' | '+' | '.join(str(m[field]) for field in ['selected','exact_duplicate_member_sets',
            'repeated_original_ids','mass_count_mismatches','host_exclusion_violations','higher_priority_neighbours_examined'])+' |')
    lines+=['','A zero host-violation count is not a production positive control when no eligible higher-priority neighbours were examined. '
        'The independent small fixtures exercise the failure branch.','',
        '## Limits of the measurement','',
        '- The finest simulation is a comparison reference, not an independent physical truth.',
        '- One matched realization isolates numerical changes but does not measure box-size or cosmology dependence.',
        '- At z=0, particle refinement C/E and timestep refinement F/T meet the working screen over broader mass ranges than force refinement E/F. '
          'The force mesh is the limiting tested setting for lower-mass z=0 haloes in this suite.',
        '- The z=0 timestep result does not extend to all higher-redshift masses; inspect the separate z=1 and z=2 ranges.',
        '- The initial E/F positions differ by at most 6.103515625e-5 Mpc/h; physical velocities match exactly. '
          'The native periodic edge guard is retained. F/T initial positions match exactly.',
        '- Native output velocities are staggered by half a timestep. The F/T velocity difference includes this output-time effect.',
        '- The normal schedule has 158 steps and T has 316. All normal endpoints and output epochs are retained.',
        '- Shape and velocity differences are reported separately; the displayed mass/Vmax screen does not certify every halo property.','',
        'Source and input receipt, membership, index and density hashes are recorded in convergence.json. '
        'Completed job allocations and billing are recorded in accounting.json.','']
    (ROOT/'CONVERGENCE.md').write_text('\n'.join(lines))
    return result
