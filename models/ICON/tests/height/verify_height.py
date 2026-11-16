#!/usr/bin/env python3
"""Coherent spherical grid fixtures and public DART/MPI operator tests.

Usage: verify_height.py BUILD OUTPUT. Run inside a two-task Slurm allocation.
Only Python's standard library and the site's ncgen utility are required.
"""
import itertools
import math
from pathlib import Path
import random
import subprocess
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from test_state_metadata import FIELDS, table

build, output = map(lambda p: Path(p).resolve(), sys.argv[1:])
output.mkdir(parents=True, exist_ok=False)

def norm(v):
    length = math.sqrt(sum(x*x for x in v))
    return tuple(x/length for x in v)
def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def ll(v): return math.atan2(v[1],v[0]),math.atan2(v[2],math.hypot(v[0],v[1]))
def vals(seq): return ','.join(format(x,'.17g') if isinstance(x,float) else str(x) for x in seq)

def mesh(global_grid, hole=False):
    rng=random.Random(7331)
    if global_grid:
        phi=(1+math.sqrt(5))/2
        vertices=[]
        for a,b in itertools.product([-1,1],repeat=2):
            vertices += [(0,a,b*phi),(a,b*phi,0),(b*phi,0,a)]
        vertices=[norm(tuple(x+rng.uniform(-.01,.01) for x in v)) for v in vertices]
        faces=[]
        for face in itertools.combinations(range(12),3):
            a,b,c=(vertices[i] for i in face)
            normal=cross(sub(b,a),sub(c,a))
            signs=[dot(normal,sub(v,a)) for i,v in enumerate(vertices) if i not in face]
            if min(signs)>0 or max(signs)<0: faces.append(face)
        assert len(faces)==20
    else:
        nx,ny=8,7
        vertices=[]
        for j in range(ny):
            for i in range(nx):
                lon=math.radians(10+i+.12*rng.random())
                lat=math.radians(30+j+.12*rng.random())
                vertices.append((math.cos(lat)*math.cos(lon),math.cos(lat)*math.sin(lon),math.sin(lat)))
        faces=[]
        for j in range(ny-1):
            for i in range(nx-1):
                # Interior hole and an open notch: the auxiliary convex hull
                # covers points that the native union must reject.
                if hole and ((i==3 and j==2) or (i>=5 and j>=4)): continue
                a=j*nx+i; b=a+1; c=a+nx; d=c+1
                faces.extend([(a,b,c),(b,d,c)])
    centers=[norm(tuple(sum(vertices[i][k] for i in f) for k in range(3))) for f in faces]
    neighbors=[[-1]*3 for f in faces]; edges={}
    for i,f in enumerate(faces):
        for k in range(3):
            edge=tuple(sorted((f[k],f[(k+1)%3])))
            if edge in edges:
                other,side=edges[edge]; neighbors[i][k]=other+1; neighbors[other][side]=i+1
            else: edges[edge]=(i,k)
    return vertices,faces,centers,neighbors

def fixture(directory,global_grid,hole,selection):
    directory.mkdir()
    vertices,faces,centers,neighbors=mesh(global_grid,hole)
    nc,nv=len(faces),len(vertices)
    rows=[0]*nc
    if not global_grid:
        rows[0]=4; rows[-1]=1; rows[nc//2]=2
    vll=list(map(ll,vertices)); cll=list(map(ll,centers))
    grid=f'''netcdf grid {{
dimensions: cell={nc}; vertex={nv}; nv=3;
variables: double clon(cell); clon:units="radian";
double clat(cell); clat:units="radian";
double vlon(vertex); vlon:units="radian";
double vlat(vertex); vlat:units="radian";
int vertex_of_cell(nv,cell); int neighbor_cell_index(nv,cell); int refin_c_ctrl(cell);
:uuidOfHGrid="height-fixture-{directory.name}"; :global_grid={int(global_grid)}; :semi_major_axis=6371229.;
data:
clon={vals(x[0] for x in cll)}; clat={vals(x[1] for x in cll)};
vlon={vals(x[0] for x in vll)}; vlat={vals(x[1] for x in vll)};
vertex_of_cell={vals(f[k]+1 for k in range(3) for f in faces)};
neighbor_cell_index={vals(f[k] for k in range(3) for f in neighbors)};
refin_c_ctrl={vals(rows)};
}}'''
    # Alternating terrain creates distinct per-column brackets around 1000 m.
    terrain=[float((i%3)*350) for i in range(nc)]
    z=[[3000.,1800.,600+t,100+t] for t in terrain]
    zi=[[3500.,2400.,1600.,350+t,t] for t in terrain]
    geometry=f'''netcdf geometry {{
dimensions: ncells={nc}; height=4; height_2=5;
variables: double z_mc(height,ncells); z_mc:units="m";
double z_ifc(height_2,ncells); z_ifc:units="m";
:uuidOfHGrid="height-fixture-{directory.name}"; :uuidOfVGrid="height-vertical";
data: z_mc={vals(v[k] for k in range(4) for v in z)};
z_ifc={vals(v[k] for k in range(5) for v in zi)};
}}'''
    declarations=[]
    for name,qty,unit in FIELDS:
        axis='height_2' if name=='w' else 'height'
        declarations.append(f'float {name}(time,{axis},ncells); {name}:units="{unit}";')
    history=f'''netcdf history {{
dimensions: ncells={nc}; height=4; height_2=5; time=UNLIMITED;
variables: double time(time); time:units="minutes since 2024-11-26 00:00:00";
{' '.join(declarations)}
:uuidOfHGrid="height-fixture-{directory.name}"; :uuidOfVGrid="height-vertical";
data: time=0;
}}'''
    for stem,text in [('grid',grid),('geometry',geometry),('history',history)]:
        cdl=directory/(stem+'.cdl'); cdl.write_text(text)
        subprocess.run(['ncgen','-o',str(directory/(stem+'.nc')),str(cdl)],check=True)
    (directory/'input.nml').write_text('''&utilities_nml
/
&ensemble_manager_nml
/
&distributed_state_nml
/
&location_nml
/
&xyz_location_nml
/
&obs_kind_nml
/
&icon_geometry_nml
horizontal_grid_file='grid.nc',vertical_geometry_file='geometry.nc',nudge_zone_width=0
/
&model_nml
 template_file='history.nc',
'''+table(selection)+'\n/\n')

summary=[]
coverage_total=[0]*9
for name,global_grid,hole,selection in [
    ('global_all',True,False,FIELDS),
    ('lam_reordered',False,False,[FIELDS[i] for i in [1,4,3,0,2]]),
    ('nonconvex_subset',False,True,[FIELDS[3]])]:
    directory=output/name
    fixture(directory,global_grid,hole,selection)
    saved=None
    for run,tasks in [('cold',1),('warm',1),('mpi',2)]:
        result=subprocess.run(['mpiexec','--bind-to','none','-n',str(tasks),str(build/'bin/height_interpolation_check')],
                              cwd=directory,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=240)
        (directory/(run+'.log')).write_text(result.stdout)
        assert result.returncode==0 and result.stdout.count('height_interpolation_check PASS task')==tasks,(name,run,result.stdout[-5000:])
        content=(directory/'operator_results.0').read_bytes()
        if run=='cold': saved=content
        if run=='warm': assert content==saved,'cold/warm mismatch'
        totals=[0]*9
        for task in range(tasks):
            path=directory/f'operator_results.{task}'
            content=path.read_text(); (directory/f'{run}.{task}').write_text(content)
            counts=list(map(int,content.split('STATUSES')[1].split()))
            coverage=list(map(int,content.split('COVERAGE')[1].split('STATUSES')[0].split()))
            coverage_total=[a+b for a,b in zip(coverage_total,coverage)]
            totals=[a+b for a,b in zip(totals,counts)]
        assert totals[0]>0 and totals[6]>0,(name,totals)
        if not global_grid: assert totals[4]>0,(name,totals)
        if hole: assert totals[3]>0,(name,totals)
        summary.append(f'{name} {run} tasks={tasks} statuses={totals}')
        print(summary[-1],flush=True)
assert all(n>0 for n in coverage_total),coverage_total
summary.append(f'Coverage nearest mask/corner mask/zero-weight masked corner/different brackets/zero data/nonzero data/partial brackets/zero-weight bracket failure/remote samples: {coverage_total}')
(output/'summary.txt').write_text('\n'.join(summary)+'\n')
