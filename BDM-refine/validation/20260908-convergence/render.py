"""Run frozen lightweight plotting/presentation sources after measured analysis."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

from common import ROOT, WORK, now, sha, write_json


def main():
    data=json.loads((ROOT/'convergence.json').read_text())
    if not data['comparisons']:raise RuntimeError('No measured comparison to render')
    started=now();bundle=Path(sys.argv[0])
    with tempfile.TemporaryDirectory(prefix='render-',dir=WORK) as folder:
        frozen=Path(folder)
        names=['common.py','plots/plot_convergence.py','plots/house_style.py','plots/chz-paper.mplstyle',
               'slides/make_slides.py','slides/beamerthemeStanford.sty','slides/beamercolorthemestanford.sty',
               'slides/beamerouterthemesimplefooter.sty','slides/durham_logo.png']
        sources={}
        for name in names:
            target=frozen/name;target.parent.mkdir(parents=True,exist_ok=True)
            if zipfile.is_zipfile(bundle):
                with zipfile.ZipFile(bundle) as archive:target.write_bytes(archive.read(name))
            else:shutil.copyfile(ROOT/name,target)
            sources[name]=sha(target)
        env=dict(os.environ,BDM_CONVERGENCE_ROOT=str(ROOT),MPLCONFIGDIR=str(frozen/'mplconfig'))
        for name in ['plots/plot_convergence.py','slides/make_slides.py']:
            command=['micromamba','run','-n','cosemu','python3','-B',str(frozen/name)]
            subprocess.run(command,env=env,check=True)
    plots=json.loads((ROOT/'plot-manifest-n300.json').read_text())
    presentation=json.loads((ROOT/'presentation-validation.json').read_text())
    write_json(ROOT/'render-validation.json',dict(completed=True,complete_campaign=data['completed'],
        started_at_utc=started,completed_at_utc=now(),source_sha256=sources,
        plot_manifest_sha256=sha(ROOT/'plot-manifest-n300.json'),
        presentation_manifest_sha256=sha(ROOT/'presentation-validation.json'),
        figure_pdf=plots['pdf'],presentation_pdf=presentation['pdf'],
        visual_review='Automated build and text/layout-log checks only; newly completed results still need visual scientific review.'))


if __name__=='__main__':main()
