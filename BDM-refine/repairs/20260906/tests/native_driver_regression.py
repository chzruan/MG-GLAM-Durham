"""Adversarial controls for the native validation driver's reviewed failures."""
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
from unittest import mock

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('native_validation',HERE.parent/'native_validation.py')
driver=importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)


def main():
    checks=[]
    reference=driver.np.ones((2,24))
    driver.catalogue_agreement(reference,reference.copy())
    for name,changed in [('different population',reference[:1]),
                          ('different physics',reference*1.02)]:
        try:driver.catalogue_agreement(reference,changed)
        except AssertionError:checks.append(name+' rejected')
        else:raise AssertionError(name+' incorrectly accepted')
    changed=reference.copy();changed[0,13]+=1
    try:driver.catalogue_agreement(reference,changed)
    except AssertionError:checks.append('different particle count rejected')
    else:raise AssertionError('changed particle count accepted')

    with tempfile.TemporaryDirectory(prefix='bdm-native-driver-control-') as tmp:
        root=Path(tmp);(root/'work').mkdir();(root/'snapshot').mkdir();build=root/'build';build.mkdir()
        executable=build/'PMP2BDM.exe'
        executable.write_text('#!/bin/bash\nsleep 30 &\nprintf "CHILD:%s\\n" "$!"\nwait\n')
        executable.chmod(0o700)
        original_popen=driver.subprocess.Popen
        def short_timeout(*args,**kwargs):
            process=original_popen(*args,**kwargs)
            original_communicate=process.communicate
            def communicate(input=None,timeout=None):
                return original_communicate(input,timeout=.2 if timeout is not None else None)
            process.communicate=communicate
            return process
        records=[];log=io.StringIO()
        with mock.patch.object(driver,'ROOT',root),mock.patch.object(driver,'BUILD',build), \
             mock.patch.object(driver.subprocess,'Popen',short_timeout), \
             mock.patch.dict(os.environ,{'BDM_AUDIT_NATIVE_LIBS':os.environ.get('LD_LIBRARY_PATH','')}):
            record,data=driver.one_replay(dict(directory=root/'snapshot',step=1),64,'baseline',1,0,log,records)
        assert record['timed_out'] and record['returncode']!=0 and data is None
        assert records==[record] and 'CHILD:' in log.getvalue()
        child=int(next(line.split(':')[1] for line in record['stdout_tail'] if line.startswith('CHILD:')))
        # A killed orphan can await init's reaping briefly; it must not run.
        proc=Path(f'/proc/{child}/stat')
        if proc.exists():assert proc.read_text().split(') ',1)[1].split()[0]=='Z'
        checks.extend(['timeout stops wrapper and finder child','timeout preserves failed record and output'])
    print(json.dumps(dict(passed=len(checks),checks=checks),indent=2))


if __name__=='__main__':main()
