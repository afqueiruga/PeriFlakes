from fenics import *
import numpy as np
import os
import tempfile
import shutil

# gmshloc = "gmsh"
gmshloc = "~/Documents/gmsh-3.0.4-Linux/bin/gmsh"
geo = os.path.abspath("crack.geo")
E = 1.0
nu= 0.25
resolutions = np.linspace(0.25,4.0,10)

nnzs = []
for clscale in resolutions:
    tmpdir = tempfile.mkdtemp()
    print "Storing temporary files in ", tmpdir
    os.system("cp "+geo+" "+tmpdir+"/crack.geo")
    os.system(gmshloc+" "+tmpdir+"/"+"crack.geo -clscale {0} - ".format(clscale))
    os.system("dolfin-convert {0}/crack2.msh {0}/crack2.xml".format(tmpdir))
    # Open up the mesh files and define integrals
    mesh = Mesh(tmpdir+"/crack2.xml")
    # Define the function space and solution/test/trial functions
    V = VectorFunctionSpace(mesh,"CG",1)
    u,tu,Du = Function(V,name='u'),TestFunction(V),TrialFunction(V)
    # Define the weak form
    n = FacetNormal(mesh)
    K = E/(3.0-6.0*nu)
    G = E/(2.0+2.0*nu)
    Dsigma = (K-2.0/3.0*G)*tr(sym(grad(Du)))*Identity(2) + 2.0*G*sym(grad(Du))
    a = inner(grad(tu),Dsigma)*dx
    K = assemble(a)
    nnzs.append(( clscale,K.nnz() ))
    # Delete the meshes
    shutil.rmtree(tmpdir)
print nnzs
