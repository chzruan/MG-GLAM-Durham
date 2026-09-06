"""Freeze and build the refined finder after production-pilot follow-up repairs.

Use the Intel2024.2 environment in validation.sbatch, then execute with
micromamba run -n cosemu python3 -B. This is a serial native build, not a simulation.
The original pilot binaries and receipts are preserved in their own directory.
"""
import importlib.util
import json
import shutil

import run_validation as v


def main():
    # Restore historical executables as labelled references, not as a claim
    # that their source equals the current repaired working tree.
    v.prepare(512., verify_current_sources=False)
    definition = importlib.util.spec_from_file_location('native_builder', v.REPAIRS/'build_native.py')
    builder = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(builder)
    builder.ROOT = v.ROOT
    builder.REPO = v.REPO
    builder.BUILD = v.WORK/'final-native-build'
    if (v.ROOT/'native-build.json').exists():
        raise FileExistsError('Refuse to overwrite the frozen final build receipt')
    builder.main()
    build = json.loads((v.ROOT/'native-build.json').read_text())
    assert build['completed']
    shutil.copy2(v.WORK/'bin/PMP2BDM.old.exe', builder.BUILD/'PMP2BDM.old.exe')
    experiment = json.loads((v.ROOT/'preparation.json').read_text())
    experiment['pilot_source_commit'] = experiment['source_commit']
    experiment.update(source_commit=build['source_commit'], new_finder_commit=build['source_commit'],
                      sources_sha256=build['source_sha256'],
                      binaries_sha256={**build['binaries_sha256'],
                          'PMP2BDM.old.exe': v.sha(builder.BUILD/'PMP2BDM.old.exe')},
                      binary_directory=str(builder.BUILD), native_build_sha256=v.sha(v.ROOT/'native-build.json'),
                      pilot_preparation_sha256=v.sha(v.ROOT/'preparation.json'),
                      final_builder_sha256=v.sha(__file__),
                      shape_repair_commit='1d871a40084ee2458a61f6c78774162f96f30eaa')
    v.write_json(v.ROOT/'main-preparation.json', experiment)
    print('Frozen final native build and main preparation', flush=True)


if __name__ == '__main__':
    main()
