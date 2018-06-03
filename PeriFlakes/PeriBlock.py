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
    def __init__(self, L,Nside, cutoff, E=1.0, nu=0.0):
        """
        Initialize a square domain of half-side-length L with Nside 
        particles along its side. cutoff is the support of the influence
        function. 
        A Peridynamics scheme has 3 hyperparameters:
        1. method: The force law
        2. weight: The form of the influence function
        3. cutoff: The scaling of the support radius
        
        The cutoff is the only hyperparamter that needs to be set
        at first; the same data and graph can be applied to many schemes in 
        solve(). The cutoff radius is used to build the neighbor list.
        """
        self.x = cf.PP.init_grid(Nside,Nside, [-L,-L], [2*L,0.0], [0.0,2*L])
        particle_Vol = (2.0*L)**2/float(Nside**2)
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
            'p_E':(np.array([E]),self.dm_GlobalSca),
            'p_nu':(np.array([nu]),self.dm_GlobalSca),
            'p_Vol':(np.array([particle_Vol]),self.dm_GlobalSca),
            'p_stab':(np.array([1.0]),self.dm_GlobalSca),
        }
        # Mark boundaries
        eps = 1.0e-10
        self.left  = cf.select_nodes(self.x, lambda a: a[0]<-L+eps )
        self.right = cf.select_nodes(self.x, lambda a: a[0]> L-eps )
        self.bottom= cf.select_nodes(self.x, lambda a: a[1]<-L+eps )
        self.top   = cf.select_nodes(self.x, lambda a: a[1]> L-eps )
        
    def setbcs(self, diri=None,neum=None):
        """
        Set the boundary conditions that will be applied by solve()
        """
        if neum is None:
            self.loaddofs = self.dm_PtVec.Get_List(self.top)[1::2]
        else:
            self.loaddofs = np.array([self.dm_PtVec.Get_List(A)[B::2]
                                      for A,B in neum],dtype=np.int).flatten()
        if diri is None:
            self.diridofs = np.array([
                self.dm_PtVec.Get_List(self.right)[0::2],
                self.dm_PtVec.Get_List(self.left)[0::2],
                self.dm_PtVec.Get_List(self.bottom)[1::2]]).flatten()
        else:
            self.diridofs = np.array([ (self.dm_PtVec.Get_List(A)[B::2] if B>=0
                                          else self.dm_PtVec.Get_List(A))
                                      for A,B in diri]).flatten()
        Nbc = len(self.diridofs)
        self.ubc = np.zeros(Nbc)

    def cutbonds(self, x0,y0,x1,y1):
        """
        Cut all of the bonds manually. 

        PeriFlakes does not actually implement the damage evolution, as
        it is meant to illustrate the properties that affect it.
        """
        hcut,hintact = cf.Filter(hb.kernel_line_intersection,
                                 self.HBond,
                                [self.data,
                                {'pts':(np.array([x0,y0,x1,y1]),
                                        cf.Dofmap_Strided(4,stride=0))}] )
        dofs = self.dm_BondSca.Get_List([e[2] for e in hcut])
        self.data['alpha'][0][dofs] = 0.0
        self.HCut = hcut

    def solve(self, method, weight, P=1.0, smoothing=""):
        """
        Solves the deformation of the block matrix given the method name and influence 
        function. The influence support is decided at initialization of the PeriBlock
        object.

        The method and weight strings key directly into the peridynamics kernels in the 
        husk. There is an additional smoothing argument, which, if not false-valued, will
        call the smoothing kernel after assembling and solving the system.

        If there are cut bonds present, self.HCut, the P argument will be used to apply
        the fluid-filled pressure force on the bonds.
        """
        # Assemble the matrix and load for the given peridynamics law
        K,R = cf.Assemble(hp.__dict__['kernel_{0}_{1}'.format(method,weight)],
                          self.HAdj, self.data,
                          {'R':(self.dm_PtVec,),
                           'K':(self.dm_PtVec,)},
                          gdim*self.NPart)
        # Assemble the fluid pressure on bonds if present
        try:
            Rp, = cf.Assemble(hb.kernel_bond_pressure,
                              self.HCut,
                              [self.data,{'p':(np.array([P]),self.dm_GlobalSca)}],
                              {'R':(self.dm_PtVec,),},
                              gdim*self.NPart)
            R += Rp
        except AttributeError:
            pass
        # Apply boundary conditions and then solve the matrix system
        cf.Apply_BC(self.diridofs,self.ubc, K,R)
        R[self.loaddofs]-= 1.0 * self.data['p_Vol'][0]**0.5/self.data['p_Vol'][0]
        u = splin.spsolve(K,R)
        # If we specified a smoothing function, post process the solution
        if smoothing:
            us, = cf.Assemble(hp.__dict__['kernel_smooth_{0}'.format(weight)],
                              self.HAdj, [self.data,{'y':(u,self.dm_PtVec)}],
                              {'ys':(self.dm_PtVec,)},
                              gdim*self.NPart)
            u = us
        # The PeriBlock object is stateless w.r.t. solution, so we return the result
        return u

    def output(self, fname, u=None):
        """
        Writes a vtk file with the bond states and, optionally, the displacement solution.
        """
        cf.GraphIO.write_graph(fname,self.HPair,self.x,
                       [('x',self.x)] +
                       ([('u',u.reshape((-1,2)))] if not u is None else []),
                       [('alpha', self.data['alpha'][0])]     
                       )
