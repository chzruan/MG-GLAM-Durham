"""Describe measured convergence, separating count limits from chosen tolerances."""
import hashlib
import math
from pathlib import Path

import numpy as np

from common import ROOT, now, read_json_snapshot, sha, verify_json_snapshot, write_json

CRITERIA=dict(minimum_bin_count=30,abundance_difference_percent=5.,
              maximum_abundance_jackknife_sigma_percent=5.,median_bound_mass_difference_percent=5.,
              median_resolved_vmax_difference_percent=2.,reference_completeness=.9)

PROPERTY_DEFINITIONS=dict(
    axis_ba='Reported intermediate/major ratio: reduced-inertia eigenvalue ratio after the legacy '
            'empirical exponent correction depending on RadRms / extended aperture, then transverse-axis reordering.',
    axis_ca='Reported minor/major ratio: reduced-inertia eigenvalue ratio after the legacy '
            'empirical exponent correction depending on RadRms / extended aperture, then transverse-axis reordering.')


def source_sha256():
    # The loader reads the actual module bytes in both ordinary and zipapp runs.
    return hashlib.sha256(__loader__.get_data(__file__)).hexdigest()


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


def diagnostics(data,comparisons):
    """Supplement the frozen acceptance decisions without changing their scope."""
    result=[]
    for row,verdict in zip(data['comparisons'],comparisons):
        passed=[i for i,b in enumerate(verdict['bins']) if b['meets_working_criteria']]
        shape={}
        for name in ['axis_ba','axis_ca']:
            usable=[row['matched_statistics'][name][i] for i in passed
                    if row['matched_statistics'][name][i]['count']>=CRITERIA['minimum_bin_count']
                    and finite(row['matched_statistics'][name][i]['q16_median_q84'][1])]
            shape[name]=dict(bins_with_enough_objects=len(usable),
                maximum_absolute_bin_median_shift_percent=max(
                    (abs(s['q16_median_q84'][1]) for s in usable),default=None))
        counts=[]
        for i,b in enumerate(verdict['bins']):
            left,right=row['left_counts'][i],row['right_counts'][i]
            scale=math.sqrt(1/left+1/right) if left>0 and right>0 else None
            ratio=row['abundance_ratio'][i]
            counts.append(dict(log10_mass_interval=b['log10_mass_interval'],eligible=b['eligible'],
                passes_full_screen=b['meets_working_criteria'],left_count=left,right_count=right,
                paired_octant_sigma_percent=100*row['abundance_ratio_jackknife8_sigma'][i]
                    if finite(row['abundance_ratio_jackknife8_sigma'][i]) else None,
                independent_count_relative_scale_percent=100*scale if scale is not None else None,
                independent_count_ratio_scale_percent=100*abs(ratio)*scale
                    if scale is not None and finite(ratio) else None))
        masses={side:data['abundances'][row[side]+'/z'+str(row['redshift'])]['mass_one']
                for side in ['coarse','reference']}
        whole=next((edges for edges,valid in zip(
            zip(data['log10_mass_edges'][:-1],data['log10_mass_edges'][1:]),row['valid_mass_bins']) if valid),None)
        result.append(dict(coarse=row['coarse'],reference=row['reference'],redshift=row['redshift'],
            nominal_particle_floor=row['particle_floor'],mass_floor_msun_h=row['mass_floor_msun_h'],
            minimum_particles_at_common_mass_floor={side:math.ceil(row['mass_floor_msun_h']/mass)
                for side,mass in masses.items()},
            first_whole_bin_log10_mass_interval=whole,
            minimum_particles_at_first_whole_bin={side:math.ceil(10**whole[0]/mass)
                for side,mass in masses.items()} if whole else None,
            unscreened_shapes_inside_passing_intervals=shape,abundance_count_diagnostics=counts))
    return result


def assess(path):
    producer_sha=source_sha256()
    data,input_sha=read_json_snapshot(path);edges=data['log10_mass_edges']
    comparisons=[evaluate(row,edges) for row in data['comparisons']]
    extra=diagnostics(data,comparisons)
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
                assess_source_sha256=producer_sha,property_definitions=PROPERTY_DEFINITIONS,
                criteria=CRITERIA,comparisons=comparisons,z0_overview=overview,
                supplemental_diagnostics=extra,
                abundance_uncertainty_interpretation=dict(
                    paired_octant_sigma='Delete-one-octant paired-ratio scatter in this realization. '
                        'Can be zero for identical octant counts; not a confidence interval or proof of zero ensemble uncertainty.',
                    independent_count_scale='Hypothetical independent Poisson counts: relative scale sqrt(1/Nleft + 1/Nright); '
                        'ratio scale additionally multiplies by |Nleft/Nright|. Actual paired catalogues are correlated. '
                        'This diagnostic omits their covariance and is neither a replacement uncertainty nor an acceptance criterion.'),
                initial_conditions=dict(generator='PMP2start.matched.exe, campaign variant of native GLAM',
                    lpt_order=1,z_init=100,
                    scope='Existing first-order campaign; not a convergence test of the default 2LPTIC workflow.'),
                interpretation='Descriptive working tolerances, not a formal confidence interval or proof of absolute accuracy. '
                    'The highest resolution is a reference. One realization and eight spatial octants cannot establish volume convergence.')
    verify_json_snapshot(path,input_sha)
    assert source_sha256()==producer_sha, 'Assessment source changed during reporting'
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
        'The table applies a descriptive working screen **only to abundance, median bound mass, median resolved Vmax and reference completeness**: '
        'at least 30 objects in each required bin; abundance within 5% with paired octant sigma at most 5%; '
        'median bound mass within 5%; median resolved Vmax within 2%; reference completeness at least 90%. '
        'These choices are adjustable. They do not certify those absolute accuracies. '
        'Scatter, individual criteria and 100/300/1000-particle floors are retained in the JSON results.','',
        '| Comparison | z | Nominal particle floor | Eligible bins | log10 mass intervals meeting the stated screen |',
        '|---|---:|---:|---:|---|']
    for row in comparisons:
        if row['particle_floor']!=300:continue
        ranges='; '.join(f'{lo:.2f}–{hi:.2f}' for lo,hi in row['meets_working_criteria_log10_mass_intervals']) or 'None'
        lines.append(f"| {row['coarse']}/{row['reference']} ({row['kind']}) | {row['redshift']} | 300 | {row['eligible_bins']} | {ranges} |")
    lines+=['','## Shape shifts inside the passing intervals','',
        'The screen does not test shape. The following are maximum absolute bin-median percentage shifts **inside its passing intervals**, '
        'with at least 30 objects for the shape statistic. They are not halo-to-halo scatter or shape acceptance limits. '
        'Both reported axis ratios include the inherited concentration-dependent empirical correction and transverse-axis reordering '
        '(PMP2linker.f90:1105–1117); they are not the uncorrected tensor ratios.','',
        '| Pair | z | max abs median delta(b/a), % | max abs median delta(c/a), % |',
        '|---|---:|---:|---:|']
    for d in extra:
        if d['nominal_particle_floor']!=300:continue
        values=[d['unscreened_shapes_inside_passing_intervals'][p]['maximum_absolute_bin_median_shift_percent']
                for p in ['axis_ba','axis_ca']]
        lines.append(f"| {d['coarse']}/{d['reference']} | {d['redshift']} | "+
                     ' | '.join('Unmeasured' if v is None else f'{v:.2f}' for v in values)+' |')
    lines+=['',('In the completed suite, E/F at z=0 reaches 5.86% in median c/a inside its mass/Vmax passing interval. '
        'F/T reaches 4.70–5.42% across the three outputs; the coarser A/C particle comparison reaches 6.24% at z=0. '
        if data['completed'] else '')+'Shape, velocity, tails and scatter require their own criteria.','',
        '## Effective particle and publication cuts','',
        'The common mass cut is max(2.5e12 Msun/h, nominal_floor × max(m_particle)). Only whole bins above that cut are used. '
        'The table below uses the nominal 300-particle floor; the values do not depend on redshift in this campaign.','',
        '| Pair | Common mass cut, Msun/h | Minimum bound particles, left / reference | First whole-bin lower log10 mass |',
        '|---|---:|---:|---:|']
    for d in extra:
        if d['nominal_particle_floor']!=300 or d['redshift']!=0:continue
        n=d['minimum_particles_at_common_mass_floor'];edge=d['first_whole_bin_log10_mass_interval']
        lines.append(f"| {d['coarse']}/{d['reference']} | {d['mass_floor_msun_h']:.7g} | {n['coarse']} / {n['reference']} | "+
                     ('None' if edge is None else f'{edge[0]:.2f}')+' |')
    lines+=['','For E/F and F/T, the publication mass cut corresponds to 1867.22 particle masses, hence at least **1868 bound particles**. '
        'Their first whole bin begins at log10 mass 12.50, requiring at least **2362 particles**. '
        'Nominal floors of 100, 300 and 1000 therefore give the same results for these pairs; they do not constitute three independent resolution tests.','',
        '## What the abundance diagnostic measures','',
        'The paired delete-one-octant jackknife measures variation of the abundance ratio across spatial omissions in this realization. '
        'It can legitimately be zero when both catalogues have identical octant counts. This does not establish zero ensemble uncertainty. '
        'The screen uses measured shifts and this diagnostic; it is not a confidence statement that the true abundance shift is below 5%.','',
        'For comparison, supplemental_diagnostics in convergence-assessment.json reports sqrt(1/Nleft + 1/Nright) as a '
        '**hypothetical independent-count relative scale**, and also that scale multiplied by the absolute ratio for comparison with delta n. '
        'The actual catalogues are paired and correlated, so their covariance must be included for a sampling-error interpretation. '
        'The independent-count scale is not used to replace the jackknife, bound the actual uncertainty, or change the acceptance decisions.','',
        'C/E at z=1, log10 mass 14.00–14.25 has identical octant counts [10,9,7,5,2,6,4,4]: '
        '47 objects in each catalogue, paired sigma 0, and independent-count relative scale 20.63%. '
        'This is one unique bin repeated at three nominal floors. At floor 300, 21 bins meet the abundance-difference cut while paired sigma exceeds 2.5%; '
        'that count includes bins failing other parts of the full screen. High-mass abundance agreement has limited statistical discrimination.','']
    lines+=['','## Membership checks','',
        'The frozen checker applies the production **priority-ordered extended-aperture rule**: a lower-priority centre must lie outside '
        'the higher-priority halo\'s reported aperture. Priority is bound mass, then stable candidate index. '
        'This is the rule in PMP2linker.f90:696–750; it is not symmetric centre exclusion or a test using unextended SO radii. '
        'The zero fields are assertions that passed; the examined count counts neighbours already inside the queried aperture and is not a broad coverage measure.','',
        '| Catalogue | Published haloes | Exact duplicate sets | Repeated IDs within a halo | Mass/count mismatches | Priority-aperture violations | In-aperture priority pairs examined |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for key,source in data['inputs'].items():
        m=source['membership_checks']
        lines.append('| '+key+' | '+' | '.join(str(m[field]) for field in ['selected','exact_duplicate_member_sets',
            'repeated_original_ids','mass_count_mismatches','host_exclusion_violations','higher_priority_neighbours_examined'])+' |')
    lines+=['','The independent review found four reverse-orientation centre-in-aperture pairs: the higher-priority centre lies inside '
        'the lower-priority aperture. This orientation is allowed by the production rule. All four separations exceed both unextended SO radii. '
        'They are not host-rule failures. The production data provide no positive control of the violation branch; fixtures remain necessary.','',
        'Exact member-set uniqueness does not require disjoint memberships. The review measured excess memberships '
        '(occurrences after the first appearance of a particle ID) at 0.06–0.13% of total memberships. '
        'This global rate does not bound the fractional mass error of an individual halo. See '
        '[the independent review](../../analysis/review-20260909-convergence/REVIEW.md) for the full scan.','',
        '## Limits of the measurement','',
        '- The finest simulation is a comparison reference, not an independent physical truth.',
        '- One matched realization isolates numerical changes but does not measure box-size or cosmology dependence.',
        ('- At z=0, particle refinement C/E and timestep refinement F/T meet the working screen over broader mass ranges than force refinement E/F. '
          'The force mesh is the limiting tested setting for lower-mass z=0 haloes in this suite.' if data['completed'] else
          '- The missing comparisons prevent a conclusion about which numerical setting limits the complete suite.'),
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
