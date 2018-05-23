import numpy as np
import cornflakes as cf
import scipy.sparse.linalg as splin

import util
import husk_peridynamics as hp
import husk_bonds as hb
gdim = 2

class PeriBlock():
    """
    This is a base class for making a square peridynamics domain
    """
    def __init__(self, L,Nside, cutoff):
        self.x = cf.PP.init_grid(Nside,Nside, [-L,-L], [2*L,0.0], [0.0,2*L])
        particle_Vol = L**2/float(Nside**2)
        self.NPart = self.x.shape[0]
        self.HPair = cf.Graphers.Build_Pair_Graph(self.x, cutoff*1.1)
        self.NBond = len(self.HPair)
        self.HBond = cf.Graphers.Build_Pair_Graph(self.x, cutoff*1.1)
        self.HBond.Add_Edge_Vertex(self.NPart)
        self.HAdj = util.Make_Bond_Adjacency(self.HBond)

        y = self.x.copy()
        alpha = np.ones(self.NBond,dtype=np.double)
        delta = np.ones(self.NPart,dtype=np.double)
        self.dm_PtVec = cf.Dofmap_Strided(gdim)
        self.dm_PtSca = cf.Dofmap_Strided(1)
        self.dm_BondSca = cf.Dofmap_Strided(1,-self.NPart)
        self.dm_GlobalSca = cf.Dofmap_Strided(1,stride=0)
        self.data = {
            'x':(self.x, self.dm_PtVec),
            'y':(y, self.dm_PtVec),
            'alpha':(alpha,self.dm_BondSca),
            'delta':(delta,self.dm_PtSca),
            'p_E':(np.array([1.0]),self.dm_GlobalSca),
            'p_nu':(np.array([0.25]),self.dm_GlobalSca),
            'p_Vol':(np.array([particle_Vol]),self.dm_GlobalSca),
            'p_stab':(np.array([1.0]),self.dm_GlobalSca),
        }
        # Mark boundaries
        eps = 1.0e-10
        self.right = cf.select_nodes(self.x, lambda a: a[0]<-L+eps )
        self.left  = cf.select_nodes(self.x, lambda a: a[0]> L-eps )
        self.bottom= cf.select_nodes(self.x, lambda a: a[1]<-L+eps )
        self.top   = cf.select_nodes(self.x, lambda a: a[1]> L-eps )
        
    def setbcs(self, diri=None,neum=None):
        if neum is None:
            self.loaddofs = self.dm_PtVec.Get_List(self.top)[1::2]
        else:
            self.loaddofs = np.array([self.dm_PtVec.Get_List(A)[B::2]
                                      for A,B in neum]).flatten()
        if diri is None:
            self.diridofs = np.array([
                self.dm_PtVec.Get_List(self.right)[0::2],
                self.dm_PtVec.Get_List(self.left)[0::2],
                self.dm_PtVec.Get_List(self.bottom)[1::2]]).flatten()
        else:
            self.diridofs = np.array([self.dm_PtVec.Get_List(A)[B::2]
                                      for A,B in diri]).flatten()
        Nbc = len(self.diridofs)
        self.ubc = np.zeros(Nbc)

    def cutbonds(self, x0,y0,x1,y1):
        damage, = cf.Assemble(hb.kernel_line_intersection,
                              self.HBond,
                              [self.data,
                               {'pts':(np.array([x0,y0,x1,y1]),
                                       cf.Dofmap_Strided(4,stride=0))}],
                              {'test':(self.dm_PtSca,)},
                              ndof=self.NBond)
        return damage
                              
    def solve(self, method, weight):
        K,R = cf.Assemble(hp.__dict__['kernel_{0}_{1}'.format(method,weight)],
                          self.HAdj, self.data,
                          {'R':(self.dm_PtVec,),
                           'K':(self.dm_PtVec,)},
                          gdim*self.NPart)
        cf.Apply_BC(self.diridofs,self.ubc, K,R)
        R[self.loaddofs]-= 1.0
        u = splin.spsolve(K,R)
        return u

    def output(self, fname, u=None):
        cf.GraphIO.write_graph(fname,self.HPair,self.x,
                       [('x',self.x)] +
                       ([('u',u.reshape((-1,2)))] if not u is None else [])
                       )
