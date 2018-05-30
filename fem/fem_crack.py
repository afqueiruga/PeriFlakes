from fenics import *
from SimDataDB import *
import os
import tempfile
import shutil
import numpy as np
import timeit

# This script calls gmsh each time, but needs v3+
#gmshloc = "gmsh"
gmshloc = "~/Documents/gmsh-3.0.4-Linux/bin/gmsh"
# gmshloc = "/Applications/Gmsh.app/Contents/MacOS/gmsh"
dbloc = "femdata.db"

# Switch to a scratch space and make a new mesh
sdb = SimDataDB(dbloc)
geo = os.path.abspath("crack.geo")
onum = 0
# Some global properties
E = 1.0
nu= 0.25
P = 0.1
a = 0.1
C = 0.0 # for the analytical solution

# We use the analytical solution
class LFM_Anal(Expression):
    def eval(self, value, x):
        Zhat = lambda z : np.sqrt( -a**2 + z**2.0 ) - z - C *z
        Z = lambda z : z/np.sqrt(-a**2 + z**2.0 ) - 1.0 - C
        ux = lambda z : P/(2.0*E/(2.0*(1.0+nu))) * (
            (2.0-4.0*nu)/2.0 * np.real(Zhat(z)) - np.imag(z) * np.imag(Z(z)))
        uy = lambda z : P/(2.0*E/(2.0*(1.0+nu))) * (
            (4.0-4.0*nu)/2.0 * np.imag(Zhat(z)) - np.imag(z) * np.real(Z(z)))
        value[0] = np.sign(x[0])*ux( np.abs(x[0]) + 1j*(np.abs(x[1])+0.0e-12) )
        value[1] = np.sign(x[1])*uy( np.abs(x[0]) + 1j*(np.abs(x[1])+0.0e-12) )
    def value_shape(self):
        return (2,)

@sdb.Decorate("lfm",
              [("clscale","FLOAT"),("poly","INT")],
              [("h","FLOAT"),("error","FLOAT"),("time","FLOAT"),("x","array"),("u","array")])
def sim(clscale,poly):
    global dx,ds, onum
    # Make a temporary directory to make a new set of meshes
    tmpdir = tempfile.mkdtemp()
    os.system("cp "+geo+" "+tmpdir+"/crack.geo")
    os.system(gmshloc+" "+tmpdir+"/"+"crack.geo -clscale {0} - ".format(clscale))
    os.system("dolfin-convert {0}/crack2.msh {0}/crack2.xml".format(tmpdir))
    # Open up the mesh files and define integrals
    mesh = Mesh(tmpdir+"/crack2.xml")
    cellids = MeshFunctionSizet(mesh,tmpdir+"/crack2_physical_region.xml")
    facetids = MeshFunctionSizet(mesh,tmpdir+"/crack2_facet_region.xml")
    dx = dx(subdomain_data=cellids)
    ds = ds(subdomain_data=facetids)
    # Define the function space and solution/test/trial functions
    V = VectorFunctionSpace(mesh,"CG",poly)
    u,tu,Du = Function(V,name='u'),TestFunction(V),TrialFunction(V)
    expr_lfm_anal = LFM_Anal(degree=0) #element=V.ufl_element())
    # Define the weak form
    n = FacetNormal(mesh)
    K = E/(3.0-6.0*nu)
    G = E/(2.0+2.0*nu)
    Dsigma = (K-2.0/3.0*G)*tr(sym(grad(Du)))*Identity(2) + 2.0*G*sym(grad(Du))
    a = inner(grad(tu),Dsigma)*dx
    L = - inner(tu, n*P)*ds(2) - inner(tu, n*P)*ds(4)
    bcs = [
        DirichletBC(V,expr_lfm_anal, facetids,3)
    ]
    # Solve it and time it
    def work():
        solve(a==L, u, bcs=bcs)
    runtime = timeit.timeit(work, number=10)
    # We need to nudge of the crack nodes so that the analytical solution can be evaluated
    # on the crack exactly.
    for l in facets(mesh):
        if facetids[l.index()]==2:
            mesh.coordinates()[l.entities(0),:] -= 1.0e-8
        if facetids[l.index()]==4:
            mesh.coordinates()[l.entities(0),:] += 1.0e-8
    # Calculate the error
    error = np.sqrt( assemble( inner( u-expr_lfm_anal, u-expr_lfm_anal)*dx ) \
                     / assemble( inner(expr_lfm_anal, expr_lfm_anal )*dx(mesh) ) )
    # Write vtk files if we want
    if False:
        File('outs/u_{0}.pvd'.format(onum)) << u
        en = project(u-expr_lfm_anal, V)
        en.rename('e','e')
        File('outs/en_{0}.pvd'.format(onum)) << en
        onum += 1
    # Delete the meshes
    shutil.rmtree(tmpdir)        
    return mesh.rmin(), error, runtime, mesh.coordinates(), u.vector().get_local()

import multiprocessing as multi
import itertools
resolutions = np.linspace(0.25,2.0,10)
orders = [1,2,3]
def f(t):
    print t
    return sim(*t)
poo = multi.Pool(processes=3)
map(f, itertools.product( resolutions,orders) ) 
