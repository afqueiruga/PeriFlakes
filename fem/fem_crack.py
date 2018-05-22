from fenics import *
from SimDataDB import *
import os
import tempfile
import shutil
import numpy as np

# This script calls gmsh each time, but needs v3+
#gmshloc = "gmsh"
gmshloc = "~/Documents/gmsh-3.0.4-Linux/bin/gmsh"
# gmshloc = "/Applications/Gmsh.app/Contents/MacOS/gmsh"
dbloc = "femdata.db"

# Switch to a scratch spack and make a new mesh
sdb = SimDataDB(dbloc)
geo = os.path.abspath("crack.geo")

@sdb.Decorate("lfm",
              [("clscale","FLOAT")],
              [("x","array"),("u","array")])
def sim(clscale):
    global dx,ds
    tmpdir = tempfile.mkdtemp()
    #os.chdir(tmpdir)
    os.system("cp "+geo+" "+tmpdir+"/crack.geo")
    print tmpdir
    #from IPython import embed ; embed()
    os.system(gmshloc+" "+tmpdir+"/"+"crack.geo -clscale {0} - ".format(clscale))
    os.system("dolfin-convert {0}/crack2.msh {0}/crack2.xml".format(tmpdir))
    
    mesh = Mesh(tmpdir+"/crack2.xml")
    cellids = MeshFunctionSizet(mesh,tmpdir+"/crack2_physical_region.xml")
    facetids = MeshFunctionSizet(mesh,tmpdir+"/crack2_facet_region.xml")

    dx = dx(subdomain_data=cellids)
    ds = ds(subdomain_data=facetids)
    
    V = VectorFunctionSpace(mesh,"CG",1)
    u,tu,Du = Function(V,name='u'),TestFunction(V),TrialFunction(V)
    
    n = FacetNormal(mesh)
    p = 1.0
    K = 1.0
    G = 1.0
    sigma = (K-2.0/3.0*G)*tr(sym(grad(u)))*Identity(2) + 2.0*G*sym(grad(u))
    f_R = inner(grad(tu),sigma)*dx - inner(tu,-n*p)*ds(2) - inner(tu,-n*p)*ds(4)
    bcs = [
        DirichletBC(V,Constant((0.0,0.0)), facetids,3)
    ]
    solve(f_R==0, u, bcs=bcs)
    
    shutil.rmtree(tmpdir)

    return mesh.coordinates(), u.vector().get_local()

import multiprocessing as multi
import itertools
resolutions = np.linspace(0.5,2.0,10)
def f(t):
    print t
    return sim(*t)
poo = multi.Pool(processes=3)
map(f, itertools.product( resolutions) ) 
