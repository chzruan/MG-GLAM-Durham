"""Replay the compiled finder predicate/component rule on real stored catalogue fields."""
import argparse
import ctypes
import json
from pathlib import Path
import re
import subprocess
import sys

import h5py
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from catalogue_core import read_ascii,pairs_within,edge_flags,components,sha256

p=argparse.ArgumentParser()
p.add_argument('source',type=Path)
p.add_argument('sidecar',type=Path)
p.add_argument('result',type=Path)
a=p.parse_args()
root=Path(__file__).resolve().parents[4]
finder=root/'PMP2linker.f90'
code=finder.read_text()
module=re.search(r'^module BdmDuplicateRules\b.*?^end module BdmDuplicateRules',code,re.I|re.M|re.S).group()
wrapper='''
module replay
use iso_c_binding
use BdmDuplicateRules
contains
subroutine flags(n,fields,out) bind(C)
integer(c_int),value :: n
real(c_double),intent(in) :: fields(8,n)
integer(c_int),intent(out) :: out(n)
integer :: i
do i=1,n
out(i)=0
if(StrictDuplicate(fields(1,i),fields(2,i),fields(3,i),fields(4,i), &
fields(5,i),fields(6,i),fields(7,i),fields(8,i)))out(i)=1
end do
end subroutine
subroutine roots(n,ne,edges,out) bind(C)
integer(c_int),value :: n,ne
integer(c_int),intent(in) :: edges(2,ne)
integer(c_int),intent(out) :: out(n)
integer :: i,ri,rj
do i=1,n
out(i)=i
end do
do i=1,ne
ri=DuplicateRoot(out,edges(1,i)+1);rj=DuplicateRoot(out,edges(2,i)+1)
out(max(ri,rj))=min(ri,rj)
end do
do i=1,n
ri=DuplicateRoot(out,i)
out(i)=ri
end do
out=out-1
end subroutine
end module
'''
build=Path(__file__).resolve().parent/'cases'
build.mkdir(exist_ok=True)
(build/'policy_replay.f90').write_text(module+wrapper)
subprocess.run(['gfortran','-O0','-g','-fcheck=all','-fPIC','-shared','policy_replay.f90','-o','policy_replay.so'],cwd=build,check=True)
lib=ctypes.CDLL(str(build/'policy_replay.so'))
double_array=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
int_array=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
lib.flags.argtypes=[ctypes.c_int,double_array,int_array]
lib.roots.argtypes=[ctypes.c_int,ctypes.c_int,int_array,int_array]
lib.flags.restype=None;lib.roots.restype=None
data,_=read_ascii(a.source)
edges,distance,dr=pairs_within(data,.5)
i,j=edges.T
fields=np.column_stack([np.einsum('ij,ij->i',dr,dr),data['Mbound'][i],data['Mbound'][j],
                       data['Nparticles'][i],data['Nparticles'][j],
                       data['Mtot'][i],data['Mtot'][j],
                       sum((data[k][i]-data[k][j])**2 for k in ['vx','vy','vz'])])
flags=np.zeros(len(edges),dtype=np.int32)
lib.flags(len(edges),fields,flags)
python_flags=edge_flags(data,edges,distance)['strict']
np.testing.assert_array_equal(flags.astype(bool),python_flags)
selected=np.ascontiguousarray(edges[flags.astype(bool)],dtype=np.int32)
parent=np.empty(len(data['x']),dtype=np.int32)
lib.roots(len(parent),len(selected),selected,parent)
drop=parent!=np.arange(len(parent))
with h5py.File(a.sidecar,'r') as f:
    np.testing.assert_array_equal(drop,f['drop_mask'][...])
    receipt=json.loads(f['receipt_json'][...].tobytes())
assert receipt['bitwise_validation']['all_surviving_columns_byte_identical']
assert sha256(receipt['output'])==receipt['output_sha256']
result=dict(source=str(a.source),source_sha256=sha256(a.source),
            finder_sha256=sha256(finder),compiled_policy_edges=len(selected),
            source_rows=len(parent),removed_rows=int(drop.sum()),
            compiled_finder_mask_equals_python_mask=True,
            all_surviving_reference_fields_byte_identical=True,
            scope='catalogue-stage replay at stored precision; no particle rerun')
a.result.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
