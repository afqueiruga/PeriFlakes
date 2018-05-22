from fenics import *
import os

clscale = 0.5

os.system("~/Documents/gmsh-3.0.4-Linux/bin/gmsh crack.geo -clscale {0} - ".format(clscale))
os.system("dolfin-convert crack2.msh crack2.xml".format(clscale))


mesh = Mesh("crack2.xml".format(clscale))
cellids = MeshFunctionSizet(mesh,"crack2_physical_region.xml".format(clscale))
facetids = MeshFunctionSizet(mesh,"crack2_facet_region.xml".format(clscale))

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
# f_K = derivative(f_R,u,Du)
bcs = [
    DirichletBC(V,Constant((0.0,0.0)), facetids,3)
]
solve(f_R==0, u, bcs=bcs)

File('u.pvd')<<u
File('faces.pvd')<<facetids
