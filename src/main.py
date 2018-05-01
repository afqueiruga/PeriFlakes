import numpy as np
import cornflakes as cf
import util

import husk_peridynamics as hp

# Params
gdim = 2
Nside = 25
L = 1.0
Ndesired = 8
cutoff = L/float(Nside)*2.5
x = cf.PP.init_grid(Nside,Nside, [-L,-L], [2*L,0.0], [0.0,2*L])
particle_Vol = L**2/float(Nside**2)
eps = 1.0e-10
boundary = cf.select_nodes(x, lambda a:
                        a[0]<-L+eps or a[0]>L-eps or 
                        a[1]<-L+eps or a[1]>L-eps )
Nbc = len(boundary)
ubc = np.zeros(Nbc)
NPart = x.shape[0]
HBond = cf.Graphers.Build_Pair_Graph(x, cutoff)
NBond = len(HBond)
HBond.Add_Edge_Vertex(NPart)
HAdj = util.Make_Bond_Adjacency(HBond)

y = x.copy()
alpha = np.ones(NBond,dtype=np.double)
delta = np.ones(NPart,dtype=np.double)

data = {
    'x':(x,cf.Dofmap_Strided(gdim)),
    'y':(y,cf.Dofmap_Strided(gdim)),
    'alpha':(alpha,cf.Dofmap_Strided(1,-NPart)),
    'delta':(delta,cf.Dofmap_Strided(1)),
    'p_E':(np.array([1.0]),cf.Dofmap_Strided(1,stride=0)),
    'p_nu':(np.array([0.25]),cf.Dofmap_Strided(1,stride=0)),
}

oo = cf.Assemble(hp.kernel_silling,HAdj,
                  data,
                  {'F':(cf.Dofmap_Strided(gdim),),
                   'K':(cf.Dofmap_Strided(gdim),)},
                  gdim*NPart) 
