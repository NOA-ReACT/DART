#!/usr/bin/env python3
"""Independent native-query comparisons, cache failures and MPI lifecycle tests.

Usage: verify_search.py DART_BUILD REFERENCE_BUILD TRIANGULATION_FIXTURES OUTPUT
Requires ncgen/ncdump and the NetCDF C library, but no Python NetCDF packages.
"""
import ctypes as C
import ctypes.util
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import time

build,reference,fixtures,output=map(lambda x:Path(x).resolve(),sys.argv[1:])
exe=build/'bin/horizontal_search_check'
ref=reference/'icon_search_reference'
output.mkdir(parents=True,exist_ok=False)
summary=[]

def run(args,directory,label,expect=None):
    start=time.monotonic()
    with (directory/(label+'.log')).open('w') as log:
        p=subprocess.Popen(list(map(str,args)),cwd=directory,stdout=log,stderr=subprocess.STDOUT)
        deadline=start+240
        while True:
            pid,status,usage=os.wait4(p.pid,os.WNOHANG)
            if pid:
                p.returncode=os.waitstatus_to_exitcode(status)
                break
            if time.monotonic()>deadline:
                p.kill();p.wait();raise AssertionError('Timeout: '+label)
            time.sleep(.05)
    text=(directory/(label+'.log')).read_text()
    if expect:
        assert expect in text and 'horizontal_search_check PASS' not in text,(label,text[-2000:])
    else:
        assert p.returncode==0,(label,text[-3000:])
    return dict(wall_seconds=time.monotonic()-start,peak_rss_KiB=usage.ru_maxrss)

def lonlat(x):
    return math.atan2(x[1],x[0]),math.atan2(x[2],math.hypot(x[0],x[1]))

def read_tri(path):
    lines=path.read_text().splitlines()
    n,t=map(int,lines[0].split())
    xyz=[tuple(map(float,s.split()[1:])) for s in lines[1:n+1]]
    tri=[tuple(int(v)-1 for v in s.split()) for s in lines[n+1:]]
    return xyz,tri

def queries(path,values):
    path.write_text(str(len(values))+'\n'+''.join(f'{a:.17g} {b:.17g}\n' for a,b in values))

def results(path):
    lines=path.read_text().splitlines()
    return [(list(map(int,lines[i].split())),list(map(float,lines[i+1].split()))) for i in range(0,len(lines),2)]

def base_nml(directory):
    (directory/'input.nml').write_text('&utilities_nml\n/\n&mpi_utilities_nml\n/\n')

names=['tetrahedron','octahedron','icosahedron','global_fibonacci','regional_symmetric',
       'regional_shuffled','dateline','bora']
for name in names:
    directory=output/name;directory.mkdir();base_nml(directory)
    original=fixtures/name
    xyz,tri=read_tri(original/'reference.txt')
    qs=[lonlat(x) for x in xyz]
    selected=tri if len(tri)<1000 else random.Random(33).sample(tri,600)
    for t in selected:
        corners=[xyz[i] for i in t]
        qs.append(lonlat([sum(p[k] for p in corners) for k in range(3)]))
        edge=lonlat([corners[0][k]+corners[1][k] for k in range(3)])
        qs.extend([edge,(edge[0]-1.e-12,edge[1]),(edge[0]+1.e-12,edge[1]),
                   (edge[0]-1.e-8,edge[1]),(edge[0]+1.e-8,edge[1])])
    qs.extend([(0.,0.),(math.pi,0.),(0.,math.pi/2),(0.,-math.pi/2)])
    query=directory/'queries.txt';queries(query,qs)
    cache=directory/'search.nc'
    cold=run([exe,original/'points.txt',cache,query,directory/'cold'],directory,'cold')
    assert 'Constructing horizontal' in (directory/'cold.log').read_text()
    before=cache.read_bytes();stamp=cache.stat().st_mtime_ns
    warm=run([exe,original/'points.txt',cache,query,directory/'warm',directory/'roundtrip.nc'],directory,'warm')
    assert 'Constructing horizontal' not in (directory/'warm.log').read_text()
    assert 'Loaded horizontal' in (directory/'warm.log').read_text()
    assert cache.stat().st_mtime_ns==stamp and cache.read_bytes()==before
    assert (directory/'roundtrip.nc').read_bytes()==before,'Non-exact NetCDF round trip'
    assert (directory/'cold.0').read_bytes()==(directory/'warm.0').read_bytes()
    native=run([ref,original/'reference.txt',query,directory/'reference.txt'],directory,'reference')
    actual=results(directory/'cold.0');expected=results(directory/'reference.txt')
    assert len(actual)==len(expected)==len(qs)
    largest=0.;success=0;retries=0;notfound=0;overflows=0
    for i,((ai,aw),(ei,ew)) in enumerate(zip(actual,expected)):
        assert ai==ei,(name,i,qs[i],ai,ei)
        err=max(abs(a-b) for a,b in zip(aw,ew));largest=max(largest,err)
        assert err<=2.e-13,(name,i,aw,ew,err)
        if ai[1]==0:
            assert min(aw)>=0 and max(aw)<=1 and abs(sum(aw)-1)<1.e-14
            success+=1;retries+=ai[6]>0
            # Independent affine-reproduction check for unclipped interiors:
            # the weighted vertices must reproduce orthogonal projection.
            if min(aw)>1.e-6:
                a,b,c=[xyz[j] for j in tri[ai[2]-1]]
                ab=[b[k]-a[k] for k in range(3)];ac=[c[k]-a[k] for k in range(3)]
                normal=[ab[1]*ac[2]-ab[2]*ac[1],ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]]
                lo,la=qs[i];x=(math.cos(lo)*math.cos(la),math.sin(lo)*math.cos(la),math.sin(la))
                scale=sum((x[k]-a[k])*normal[k] for k in range(3))/sum(v*v for v in normal)
                projected=[x[k]-scale*normal[k] for k in range(3)]
                reconstructed=[sum(w*v[k] for w,v in zip(aw,[a,b,c])) for k in range(3)]
                assert max(abs(x-y) for x,y in zip(projected,reconstructed))<1.e-11
        notfound+=ai[1]==1;overflows+=ai[1]==3
    summary.append(dict(case=name,queries=len(qs),success=success,retries=retries,notfound=notfound,
                        candidate_overflows=overflows,max_weight_difference=largest,cache_bytes=len(before),
                        cold=cold,warm=warm,reference=native))
    print('PASS native/fresh/cached',name,summary[-1],flush=True)

# Mutate only copies of a valid small cache via the C API.
lib=C.CDLL(ctypes.util.find_library('netcdf'))
def mutate(path,var=None,value=None,att=None):
    nc=C.c_int();assert lib.nc_open(str(path).encode(),1,C.byref(nc))==0
    if att:
        x=C.c_int(value)
        assert lib.nc_put_att_int(nc,-1,att.encode(),4,C.c_size_t(1),C.byref(x))==0
    else:
        vid=C.c_int();assert lib.nc_inq_varid(nc,var.encode(),C.byref(vid))==0
        rank=C.c_int();assert lib.nc_inq_varndims(nc,vid,C.byref(rank))==0
        idx=(C.c_size_t*rank.value)(*([0]*rank.value))
        x=C.c_double(value) if var in ['cell_lon','xyz','obj_min'] else C.c_int(value)
        fn=lib.nc_put_var1_double if isinstance(x,C.c_double) else lib.nc_put_var1_int
        assert fn(nc,vid,idx,C.byref(x))==0
    assert lib.nc_close(nc)==0

small=output/'tetrahedron';points=fixtures/'tetrahedron/points.txt';query=small/'queries.txt'
failures=[('schema',None,99,'schema_version','Unsupported search cache schema'),
          ('algorithm',None,99,'algorithm_version','Unsupported search algorithm'),
          ('depth',None,9,'octree_depth','Incompatible octree depth'),
          ('coordinates','cell_lon',0.,None,'ordered cell coordinates mismatch'),
          ('vertex','vertex',0,None,'Invalid cached triangle index'),
          ('mapping','point_cell',0,None,'Invalid cached cell mapping'),
          ('object','object_idx',0,None,'Invalid cached object index'),
          ('box','box',0,None,'Invalid cached box ID'),
          ('finite','xyz',float('nan'),None,'Nonfinite cached xyz'),
          ('bounds','obj_min',float('nan'),None,'Nonfinite cached object bounds')]
for name,var,value,att,error in failures:
    case=output/('invalid_'+name);case.mkdir();base_nml(case)
    cache=case/'search.nc';shutil.copyfile(small/'search.nc',cache)
    mutate(cache,var,value,att)
    run([exe,points,cache,query,case/'result'],case,'run',error)
# Schema/offset failures, and a trusted but intentionally dense search tree
# to exercise the native 99-candidate limit without inventing a fallback.
def cdl_case(name,source,transform):
    directory=output/name;directory.mkdir();base_nml(directory)
    cdl=subprocess.check_output(['ncdump','-p','17,17',str(source)],text=True)
    (directory/'cache.cdl').write_text(transform(cdl))
    subprocess.run(['ncgen','-o',str(directory/'search.nc'),str(directory/'cache.cdl')],check=True)
    return directory
case=cdl_case('invalid_shape',small/'search.nc',lambda s:s.replace('pair = 2','pair = 3'))
run([exe,points,case/'search.nc',query,case/'result'],case,'run','Invalid pair dimension')
case=cdl_case('invalid_offsets',small/'search.nc',
             lambda s:re.sub(r'\n box =.*?;', '\n box = 1, 299594, 0, 5;',s,flags=re.S))
run([exe,points,case/'search.nc',query,case/'result'],case,'run','Invalid cached offset endpoints')
def dense(s):
    s=re.sub(r'box_entry = \d+', 'box_entry = 2',s)
    for name,values in [('box','1,299594,1,253'),('obj_min',','.join(['-1']*756)),
                        ('obj_max',','.join(['1']*756))]:
        s=re.sub(r'\n '+name+r' =.*?;', '\n '+name+' = '+values+';',s,flags=re.S)
    return s
case=cdl_case('candidate_limit',output/'global_fibonacci/search.nc',dense)
queries(case/'queries.txt',[(0.,0.)])
run([exe,fixtures/'global_fibonacci/points.txt',case/'search.nc',case/'queries.txt',case/'result'],case,'run')
assert results(case/'result.0')[0][0][1]==3
# Same UUID but reordered input coordinates.
case=output/'reordered';case.mkdir();base_nml(case)
lines=points.read_text().splitlines();lines[1],lines[2]=lines[2],lines[1]
(case/'points.txt').write_text('\n'.join(lines)+'\n')
run([exe,case/'points.txt',small/'search.nc',query,case/'result'],case,'run','ordered cell coordinates mismatch')
# Truncated/partial data must not trigger reconstruction or replacement.
case=output/'truncated';case.mkdir();base_nml(case)
(case/'search.nc').write_bytes((small/'search.nc').read_bytes()[:100])
run([exe,points,case/'search.nc',query,case/'result'],case,'run','ERROR')
assert (case/'search.nc').stat().st_size==100
# An abandoned temporary sibling must never be mistaken for the final cache.
case=output/'abandoned_temp';case.mkdir();base_nml(case)
(case/'search.nc.tmp.abandoned').write_bytes(b'partial')
run([exe,points,case/'search.nc',query,case/'result'],case,'run')
assert (case/'search.nc.tmp.abandoned').read_bytes()==b'partial'
# A nonexistent parent makes publication fail without losing in-memory results.
case=output/'publication_failure';case.mkdir();base_nml(case)
run([exe,points,case/'missing/search.nc',query,case/'result'],case,'run')
assert 'Could not publish cache' in (case/'run.log').read_text()
assert (case/'result.0').read_bytes()==(small/'cold.0').read_bytes()
# Native invalid construction remains a controlled failure and never publishes.
case=output/'fine_rejection';case.mkdir();base_nml(case)
run([exe,fixtures/'regional_fine/points.txt',case/'search.nc',query,case/'result'],case,'run','invalid triangle count')
assert not (case/'search.nc').exists()
# Invalid query arguments do not enter octree traversal.
case=output/'invalid_queries';case.mkdir();base_nml(case)
queries(case/'queries.txt',[(float('nan'),0.),(0.,2.),(7.,0.)])
run([exe,points,small/'search.nc',case/'queries.txt',case/'result'],case,'run')
assert all(row[0][1]==2 for row in results(case/'result.0'))

home=Path.home()
grid=home/'datasets/grids/icon/GRC_Rec_R17B5_Av1/raw/grid/DOM01/Bora_DOM01.nc'
vertical=home/'experiments/ICON/Storm_Bora/icon-dart-tests/output/static_geometry/geometry_20241126T000000Z_ML.nc'
def geometry_nml(case):
    base_nml(case)
    with (case/'input.nml').open('a') as f:
        f.write(f"&icon_geometry_nml\n horizontal_grid_file='{grid}',\n vertical_geometry_file='{vertical}',\n nudge_zone_width=10\n/\n")
case=output/'mpi_geometry';case.mkdir();geometry_nml(case)
cache=case/'search.nc';query=output/'bora/queries.txt'
for mode in ['cold','warm']:
    run(['mpiexec','--bind-to','none','-n','2',exe,'--geometry',cache,query,case/mode],case,mode)
    assert (case/(mode+'.0')).read_bytes()==(case/(mode+'.1')).read_bytes()==(output/'bora/cold.0').read_bytes()
    log=(case/(mode+'.log')).read_text()
    assert log.count('Constructing horizontal')==(1 if mode=='cold' else 0)
# Incompatible cache UUID, even when ordered coordinates match.
case2=output/'uuid';case2.mkdir();geometry_nml(case2)
run([exe,'--geometry',output/'bora/search.nc',query,case2/'result'],case2,'run','horizontal UUID mismatch')
# Two independent cold starts publishing to the same path (Bora ensures overlap).
case=output/'concurrent';case.mkdir();geometry_nml(case)
processes=[]
for i in range(2):
    log=(case/f'job{i}.log').open('w')
    processes.append((subprocess.Popen([str(exe),'--geometry',str(case/'search.nc'),str(query),str(case/f'job{i}')],
                                      cwd=case,stdout=log,stderr=subprocess.STDOUT),log))
for proc,log in processes:
    assert proc.wait(timeout=120)==0;log.close()
assert (case/'job0.0').read_bytes()==(case/'job1.0').read_bytes()==(output/'bora/cold.0').read_bytes()
assert not list(case.glob('search.nc.tmp.*'))
assert sum('Saved horizontal search:' in (case/f'job{i}.log').read_text() for i in range(2))==1
summary.append(dict(cache_failure_cases=len(failures)+9,mpi_tasks=2,concurrent_publishers=2,passed=True))
(output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('PASS all search/cache/MPI checks',flush=True)
