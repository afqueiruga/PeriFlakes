import numpy as np
import cornflakes as cf
import scipy.sparse.linalg as splin

import util
import husk_peridynamics as hp

# Params
gdim = 2
Nside = 25
L = 1.0
Ndesired = 8
cutoff = L/float(Nside)*2.5

# Place the particles and make the mesh
x = cf.PP.init_grid(Nside,Nside, [-L,-L], [2*L,0.0], [0.0,2*L])
particle_Vol = L**2/float(Nside**2)
NPart = x.shape[0]
HBond = cf.Graphers.Build_Pair_Graph(x, cutoff)
NBond = len(HBond)
HBond.Add_Edge_Vertex(NPart)
HAdj = util.Make_Bond_Adjacency(HBond)

# Allocate the data and make the cornflakes structures
y = x.copy()
alpha = np.ones(NBond,dtype=np.double)
delta = np.ones(NPart,dtype=np.double)
dm_PtVec = cf.Dofmap_Strided(gdim)
dm_PtSca = cf.Dofmap_Strided(1)
dm_BondSca = cf.Dofmap_Strided(1,-NPart)
dm_GlobalSca = cf.Dofmap_Strided(1,stride=0)
data = {
    'x':(x, dm_PtVec),
    'y':(y, dm_PtVec),
    'alpha':(alpha,dm_BondSca),
    'delta':(delta,dm_PtSca),
    'p_E':(np.array([1.0]),dm_GlobalSca),
    'p_nu':(np.array([0.25]),dm_GlobalSca),
}

# Mark boundaries
eps = 1.0e-10
right = cf.select_nodes(x, lambda a: a[0]<-L+eps )
left  = cf.select_nodes(x, lambda a: a[0]> L-eps )
bottom= cf.select_nodes(x, lambda a: a[1]<-L+eps )
top   = cf.select_nodes(x, lambda a: a[1]> L-eps )

loaddofs = dm_PtVec.Get_List(top)[1::2]
dirrdofs = np.array([
    dm_PtVec.Get_List(right)[0::2],
    dm_PtVec.Get_List(left)[0::2],
    dm_PtVec.Get_List(bottom)[1::2]]).flatten()
Nbc = len(dirrdofs)
ubc = np.zeros(Nbc)

# Make the linear system
F,K = cf.Assemble(hp.kernel_silling,HAdj,
                  data,
                  {'F':(cf.Dofmap_Strided(gdim),),
                   'K':(cf.Dofmap_Strided(gdim),)},
                  gdim*NPart) 
cf.Apply_BC(dirrdofs,ubc, K,F)
F[loaddofs]-= 1.0
