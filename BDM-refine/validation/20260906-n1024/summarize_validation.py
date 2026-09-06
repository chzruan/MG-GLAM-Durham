"""Generate scientific summary and Beamer numbers from completed evidence.

Use micromamba run -n cosemu python3 -B. This does not infer raw membership
identity for inline catalogues, or convergence from one simulation.
"""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np

ROOT = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block := stream.read(8*1024**2):
            h.update(block)
    return h.hexdigest()


def number(value):
    return f'{int(value):,}'.replace(',', '{,}')


def atomic_text(path, text):
    staged = path.with_name(path.name+'.tmp')
    with staged.open('w') as stream:
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    staged.replace(path)


def main():
    inputs = {name: ROOT/path for name, path in dict(
        simulation='main-simulation.json', validation='main-validation.json',
        threading='main-thread-control/results.json', accounting='accounting.json',
        plots='main-property-plots.json', comparison='slides/figs/comparison-summary.json',
        figure_receipt='slides/figs/comparison-figure-receipt.json',
        metadata='comparison-metadata.json').items()}
    reports = {key: json.loads(path.read_text()) for key, path in inputs.items()}
    simulation, validation, control = (reports[k] for k in ['simulation', 'validation', 'threading'])
    comparison, accounting = reports['comparison'], reports['accounting']
    assert all(reports[k]['completed'] for k in ['simulation', 'validation', 'threading', 'plots'])
    assert validation['simulation_receipt_sha256'] == sha(inputs['simulation'])
    assert control['simulation_receipt_sha256'] == sha(inputs['simulation'])
    assert validation['spec']['sources_sha256'] == simulation['spec']['sources_sha256']
    assert control['source_sha256'] == simulation['spec']['sources_sha256']
    assert control['fixed_density_controls_passed']
    assert all(pair['byte_identical'] for pair in control['fixed_density_comparisons'])
    assert control['restoration']['fixed_density_verified'] == 6
    assert control['restoration']['particle_restoration_verified'] == 6
    assert control['restoration']['exact_original_particle_bits']
    assert comparison['input_npz_sha256'] == sha(ROOT/'main-verified-comparison-catalogues.npz')
    assert comparison['metadata_sha256'] == sha(inputs['metadata'])
    assert reports['metadata']['validation_receipt_sha256'] == sha(inputs['validation'])
    assert reports['metadata']['simulation_receipt_sha256'] == sha(inputs['simulation'])
    assert reports['plots']['figure_receipt_sha256'] == sha(inputs['figure_receipt'])
    assert reports['figure_receipt']['summary_sha256'] == sha(inputs['comparison'])
    for name, expected in reports['figure_receipt']['figures_sha256'].items():
        assert sha(ROOT/'slides/figs'/name) == expected
    assert accounting['all_finished']
    active_jobs = [str(reports[key]['job_id']) for key in ['simulation', 'validation', 'threading', 'plots']]
    assert all(accounting['jobs'][job]['state'] == 'COMPLETED' for job in active_jobs)

    physics = {}
    normal = []
    for record in validation['records']:
        assert record['scientific_verification_completed']
        if record['variant'] != 'old':
            normal.append(dict(z=record['z'], variant=record['variant'], threads=record['threads'],
                               byte_identical_to_inline=record['byte_identical_to_inline'],
                               **record['normal_density_comparison']))
        if record['variant'] != 'members':
            continue
        z = record['z']
        membership = record['membership']
        for key in ['exact_duplicate_member_sets', 'repeated_original_ids',
                    'mass_count_mismatches', 'host_exclusion_violations']:
            assert membership[key] == 0
        epoch = comparison['epochs'][f'z{z}']
        provenance = reports['metadata']['plotted_refined_catalogues'][f'z{z}']
        assert provenance['catalogue'] == record['catalogue']
        assert provenance['membership'] == membership
        assert epoch['catalogues']['new']['rows'] == membership['selected']
        assert epoch['quality']['new']['axis_order']['invalid_count'] == 0
        assert epoch['catalogues']['new']['nonfinite_any_rows'] == 0
        physics[f'z{z}'] = dict(old_rows=epoch['catalogues']['old']['rows'],
            refined_rows=membership['selected'], matched_pairs=epoch['matching']['pairs'],
            count_change_percent=100*(membership['selected']/epoch['catalogues']['old']['rows']-1),
            exact_duplicate_sets=0, host_violations=0, membership=membership,
            candidate_diagnostics=record['diagnostics'],
            catalogue_sha256=record['catalogue']['sha256'],
            scope='Exact plotted standalone refined membership replay; not an inferred inline membership')
    assert set(physics) == {'z0', 'z1', 'z2'}
    state = next(r for r in validation['records'] if r['variant'] == 'state')
    assert state['two_calls_bitwise_identical']

    prepared = ROOT/'slides/figs/comparison-plot-ready.npz'
    assert sha(prepared) == comparison['plot_ready_sha256']
    headline_properties = {}
    threshold = 10**12.5
    with np.load(prepared, allow_pickle=False) as arrays:
        cohort = (arrays['z0__mass__old'] >= threshold) & (arrays['z0__mass__new'] >= threshold)
        for name, prop in comparison['epochs']['z0']['properties'].items():
            change = arrays[f'z0__{name}__change']
            valid = cohort & arrays[f'z0__{name}__valid'] & np.isfinite(change)
            values = change[valid]
            headline_properties[name] = dict(valid_pairs=len(values),
                difference_definition=prop['difference'],
                percentiles_16_50_84=np.percentile(values, [16, 50, 84]).tolist() if len(values) else None)
    z0 = physics['z0']
    baseline = next(r for r in validation['records'] if r['variant'] == 'baseline' and r['threads'] == 64)
    old_baseline = next(r for r in validation['records'] if r['variant'] == 'old' and r['z'] == 0)
    timing = dict(new_seconds=baseline['elapsed_seconds'], old_seconds=old_baseline['elapsed_seconds'],
                  new_over_old=baseline['elapsed_seconds']/old_baseline['elapsed_seconds'], threads=64,
                  scope='Single standard standalone run of each finder on the same z0 snapshot; '
                        'includes read/density/finder/output, excludes Python input hashing and diagnostic dumps')
    controlled_times = {threads: [r['seconds'] for r in control['native_finder_timings']
        if r['tag'].startswith(f'd64-t{threads}-')] for threads in [32, 64]}
    assert all(len(values) == 2 and all(value > 0 for value in values) for values in controlled_times.values())
    controlled_medians = {threads: float(np.median(values)) for threads, values in controlled_times.items()}
    speedup = controlled_medians[32]/controlled_medians[64]
    parallel = dict(seconds=controlled_times, median_seconds=controlled_medians,
                    speedup_32_to_64=speedup, efficiency_relative_to_doubling=speedup/2,
                    scope='Two repeats per thread count on the same first immutable density field '
                        'within the diagnostic allocation; includes BDM(0) only, with probe memory overhead')
    summary = dict(completed=False, created_at_utc=datetime.now(timezone.utc).isoformat(),
        generator_sha256=sha(__file__), input_sha256={k: sha(p) for k, p in inputs.items()},
        physics=physics, normal_density_comparisons=normal,
        fixed_density_control=dict(job_id=control['job_id'], passed=True,
            between_density_comparison=control['between_density_comparison'],
            density_difference=control['density_difference']),
        headline_property_cohort=dict(z=0, minimum_old_and_new_mass_msun_h=threshold,
            statement='Conservative positional matches, with both old and new catalogued masses '
                      'above the z0 literature guidance line; property-specific validity/shape filters apply.'),
        headline_properties=headline_properties, standard_finder_timing=timing,
        controlled_thread_timing=parallel,
        active_jobs=active_jobs, total_billed_core_hours=accounting['total_billed_core_hours'],
        limitations=comparison['limitations'])

    rows = '\n'.join(f'{z} & {number(physics[f"z{z}"]["old_rows"])} & '
        f'{number(physics[f"z{z}"]["refined_rows"])} & 0 & 0 \\\\' for z in [2, 1, 0])
    changed = [r['changed_rows'] for r in normal if r.get('same_shape')]
    if all(r['byte_identical_to_inline'] for r in normal):
        normal_statement = 'Ordinary replays also match the inline outputs in this realization.'
    elif all(r.get('same_shape') for r in normal):
        normal_statement = (f'Ordinary density recomputation changes up to {number(max(changed))} '
            'printed rows per replay; see the numerical appendix.')
    else:
        normal_statement = 'Ordinary density recomputation also changes selection or order; '
        normal_statement += 'the separate replay differences are retained explicitly.'
    between = control['between_density_comparison']
    if between.get('first_rows') == between.get('second_rows') and 'changed_rows' in between:
        assert len(between['rows']) == between['changed_rows']
        count_changes = sum(row['first'][13] != row['second'][13] for row in between['rows'])
        count_statement = ('one bound-particle count' if count_changes == 1 else
                           f'{number(count_changes)} bound-particle counts')
        numerical_appendix = (f'$1024^3$: changing the density field changes '
            f'{number(between["changed_rows"])} of {number(between["first_rows"])} rows, '
            f'including {count_statement}. Each fixed-field 32/64-thread group is identical.')
    else:
        numerical_appendix = ('At $1024^3$, each fixed-field 32/64-thread group is identical; '
                              'different density fields produce different selections or ordering.')
    def median(name):
        values = headline_properties[name]['percentiles_16_50_84']
        assert values is not None and headline_properties[name]['valid_pairs'] >= 20
        return values[1]
    main_job = accounting['jobs'][str(simulation['job_id'])]
    macros = dict(
        ValidationHeadline='Membership and host checks pass',
        ValidationRows=rows,
        StateResult='Eight $z=0$ calls preserve all six particle arrays bit for bit, including sizes and counts.',
        ThreadResult='At $z=0$, each tested fixed density field gives identical catalogues at 32 and 64 threads.',
        NormalDensityResult=normal_statement,
        ConclusionsHeadline='Abundance and halo properties change',
        CountConclusion=f'At $z=0$: {number(z0["old_rows"])} to '
            f'{number(z0["refined_rows"])} haloes ({z0["count_change_percent"]:+.1f}\\%). '
            'Several physics and configuration repairs contribute.',
        PropertyConclusion=f'For matched $z=0$ haloes with both masses above $10^{{12.5}}\\,\\msunh$, '
            f'median changes are {median("mass"):+.1f}\\% in mass, {median("radius"):+.1f}\\% '
            f'in radius and {median("vmax"):+.2f}\\% in $V_\\mathrm{{max}}$.',
        NumericalConclusion='No identical bound sets or host-exclusion violations survive in the checked refined catalogues.',
        MainNumericalAppendix=numerical_appendix,
        ResourceResult=f'Shared COSMA8 allocations. Simulation job: {main_job["elapsed_seconds"]/60:.1f} min on 64 cores; '
            f'batch peak {main_job["maxrss_gib"]:.1f} GiB, allocated 192 GiB.',
        TimingResult=f'Standard $z=0$ replay at 64 threads: pre-audit {timing["old_seconds"]:.1f} s, '
            f'refined {timing["new_seconds"]:.1f} s ({timing["new_over_old"]:.2f}$\\times$). '
            'The refined finder performs additional SO and iterative-unbinding work.',
        ParallelResult=f'Controlled $32\\rightarrow64$ threads: {speedup:.2f}$\\times$ speedup '
            f'({100*speedup/2:.0f}\\% of ideal doubling); two repeats per thread count on one fixed field.',
        RefinedCommit=simulation['spec']['new_finder_commit'][:12],
        JobResult='Simulation, replay, fixed-field and plot jobs: '+', '.join(active_jobs)+'.')
    text = '% Generated only from completed, verified evidence by summarize_validation.py\n'
    text += '\n'.join('\\newcommand{\\'+key+'}{'+value+'}' for key, value in macros.items())+'\n'
    atomic_text(ROOT/'slides/validation_results.tex', text)
    summary.update(completed=True, results_tex_sha256=sha(ROOT/'slides/validation_results.tex'))
    atomic_text(ROOT/'validation-summary.json', json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(physics={k: {n: r[n] for n in ['old_rows', 'refined_rows', 'matched_pairs',
          'count_change_percent']} for k, r in physics.items()}, headline_properties=headline_properties,
          timing=timing), indent=2))


if __name__ == '__main__':
    main()
