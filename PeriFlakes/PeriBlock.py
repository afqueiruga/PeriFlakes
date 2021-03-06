import numpy as np
import cornflakes as cf
import scipy.sparse.linalg as splin

import util
import husk_peridynamics as hp
import husk_bonds as hb
import husk_ficticious as hf
gdim = 2
eps = 1.0e-12

class PeriBlock():
    """
    This is a base class for making a square peridynamics domain
    """
    def __init__(self, L,Nside, delta, E=1.0, nu=0.0, ficticious=False, ficticious_size=1.0):
        """
        Initialize a square domain of half-side-length L with Nside
        particles along its side. delta is the support of the influence
        function.
        A Peridynamics scheme has 3 hyperparameters:
        1. method: The force law
        2. weight: The form of the influence function
        3. delta: The scaling of the support radius

        The delta is the only hyperparamter that needs to be set
        at first; the same data and graph can be applied to many schemes in
        solve(). The delta radius is used to build the neighbor list.
        """
        if ficticious:
            # The band shouldn't mess up the gridding
            h = 2.0*L/float(Nside-1)
            # ceil would give us one more, but it wouldn't be inside
            # any particles' support. Then it's BC wouldn't be obvious
            Ndel = int(np.floor(delta/h)*ficticious_size)
            band = Ndel * h
            Ltot = L+band
            Nside = Nside + Ndel*2
        else:
            Ltot = L
        self.ficticious = ficticious

        # Make a grid and the graphs
        self.x = cf.PP.init_grid(Nside,Nside, [-Ltot,-Ltot], [2*Ltot,0.0], [0.0,2*Ltot])
        particle_Vol = (2.0*L)**2/float(Nside**2)
        self.NPart = self.x.shape[0]
        self.HPair = cf.Graphers.Build_Pair_Graph(self.x, delta*1.1)
        self.NBond = len(self.HPair)
        self.HBond = cf.Graphers.Build_Pair_Graph(self.x, delta*1.1)
        self.HBond.Add_Edge_Vertex(self.NPart)
        self.HAdj = util.Make_Bond_Adjacency(self.HBond)

        # The ficticious nodes are set aside in a different graph, because
        # they will apply a different operator than the PD stencil
        if ficticious:
            self.FictNodes = set()
            self.HFict = cf.Hypergraph()
            Hnew = cf.Hypergraph()
            for e in self.HAdj:
                xi = self.x[e[0],:]
                if xi[0]<L+eps and xi[0]>-L-eps and xi[1]>-L-eps and xi[1]<L+eps:
                    Hnew.Push_Edge(e)
                else:
                    self.HFict.Push_Edge(e)
                    self.FictNodes.add(e[0])
            self.HAdj = Hnew

            # Make the FDM stencils
            self.HFictStencil4 = cf.Hypergraph()
            self.HFictStencil3 = cf.Hypergraph()
            for e in self.HFict:
                x0 = self.x[e[0],:]
                bodyvertices = []
                for i in e[1:(len(e)/2+1)]:
                    if x0[1] > L:
                        if self.x[i,1] < L+eps:
                            bodyvertices.append(i)
                    elif x0[1] <- L:
                        if self.x[i,1] > -L-eps:
                            bodyvertices.append(i)
                    elif x0[0] > L:
                        if self.x[i,0] < L+eps:
                            bodyvertices.append(i)
                    elif x0[0] <- L:
                        if self.x[i,0] > -L-eps:
                            bodyvertices.append(i)

                bodyvertices.sort( key=lambda i: np.linalg.norm(self.x[i,:]-x0) )
                if len(bodyvertices)==2:
                    self.HFictStencil3.Push_Edge([e[0]]+bodyvertices)
                elif len(bodyvertices)==3:
                    self.HFictStencil4.Push_Edge([e[0]]+bodyvertices)

        # Make the data arrays, vertex-to-data mappings, and the data dictionary
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
            'load':(np.zeros(self.x.shape),self.dm_PtVec)
        }

        # Mark boundaries
        if self.ficticious:
            self.fleft  = cf.select_nodes(self.x, lambda a: a[0]<-L-h+eps )
            self.fright = cf.select_nodes(self.x, lambda a: a[0]> L+h-eps )
            self.fbottom= cf.select_nodes(self.x, lambda a: a[1]<-L-h+eps )
            self.ftop   = cf.select_nodes(self.x, lambda a: a[1]> L+h-eps )
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
            #self.loaddofs = np.array([self.dm_PtVec.Get_List(A)[B::2]
            #                          for A,B in neum],dtype=np.int).flatten()
            for A,B in neum:
                if self.ficticious:
                    for a in A:
                        if a in self.FictNodes:
                            self.data['load'][0][a,:] = B
                else:
                    self.data['load'][0][A,:] = B
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

    def _assemble_KR(self, method, weight, stab=0.0):
        K,R = cf.Assemble(hp.__dict__['kernel_{0}_{1}'.format(method,weight)],
                          self.HAdj,
                          [self.data,{'p_stab':(np.array([stab]),self.dm_GlobalSca)}],
                          {'R':(self.dm_PtVec,), 'K':(self.dm_PtVec,)},
                          gdim*self.NPart)
        return K,R
    def _assemble_KR_fict(self, K,R, fictmet, method, weight, stab=0.0):
        if fictmet == "trivial":
            Kf,Rf = cf.Assemble(hp.__dict__['kernel_{0}_{1}'.format(method,weight)],
                              self.HFict,
                              [self.data,{'p_stab':(np.array([stab]),self.dm_GlobalSca)}],
                              {'R':(self.dm_PtVec,), 'K':(self.dm_PtVec,)},
                              gdim*self.NPart)
            Rf -= self.data['load'][0].ravel() * self.data['p_Vol'][0]**0.5/self.data['p_Vol'][0]
            return K+Kf,R+Rf
        elif fictmet=="bobaru":
            K3,R3 = cf.Assemble(hf.kernel_bobaru_n3, self.HFictStencil3,
                                self.data,
                                {'R':(self.dm_PtVec,), 'K':(self.dm_PtVec,)},
                                gdim*self.NPart)
            K4,R4 = cf.Assemble(hf.kernel_bobaru_n, self.HFictStencil4,
                                self.data,
                                {'R':(self.dm_PtVec,), 'K':(self.dm_PtVec,)},
                                gdim*self.NPart)
            Kb = K3 + K4
            Rb = R3 + R4
            Kt = K.copy()
            Rt = R.copy()
            for e in self.FictNodes:
                Kt[self.dm_PtVec.Get_List([e]),:]=Kb[self.dm_PtVec.Get_List([e]),:]
                Rt[self.dm_PtVec.Get_List([e])] = Rb[self.dm_PtVec.Get_List([e])]
            return Kt,Rt
        elif fictmet=="bobaru_F":
            K3,R3 = cf.Assemble(hf.kernel_bobaru_F3, self.HFictStencil3,
                                self.data,
                                {'R':(self.dm_PtVec,), 'K':(self.dm_PtVec,)},
                                gdim*self.NPart)
            K4,R4 = cf.Assemble(hf.kernel_bobaru_F, self.HFictStencil4,
                                self.data,
                                {'R':(self.dm_PtVec,), 'K':(self.dm_PtVec,)},
                                gdim*self.NPart)
            Kb = K3 + K4
            Rb = R3 + R4
            Kt = K.copy()
            Rt = R.copy()
            for e in self.FictNodes:
                Kt[self.dm_PtVec.Get_List([e]),:]=Kb[self.dm_PtVec.Get_List([e]),:]
                Rt[self.dm_PtVec.Get_List([e])] = Rb[self.dm_PtVec.Get_List([e])]
            return Kt,Rt
        elif fictmet =="both":
            Kp,Rp = cf.Assemble(hp.__dict__['kernel_{0}_{1}'.format(method,weight)],
                              self.HFict,
                              [self.data,{'p_stab':(np.array([stab]),self.dm_GlobalSca)}],
                              {'R':(self.dm_PtVec,), 'K':(self.dm_PtVec,)},
                              gdim*self.NPart)
            K3,R3 = cf.Assemble(hf.kernel_bobaru_n3, self.HFictStencil3,
                                self.data,
                                {'R':(self.dm_PtVec,), 'K':(self.dm_PtVec,)},
                                gdim*self.NPart)
            K4,R4 = cf.Assemble(hf.kernel_bobaru_n, self.HFictStencil4,
                                self.data,
                                {'R':(self.dm_PtVec,), 'K':(self.dm_PtVec,)},
                                gdim*self.NPart)
            Kb = K3+K4
            Rb = R3+R4
            Kt = K + Kp
            Rt = R + Rp
            for e in self.FictNodes:
                Kt[self.dm_PtVec.Get_List([e]),:]=Kb[self.dm_PtVec.Get_List([e]),:]
                Rt[self.dm_PtVec.Get_List([e])] = Rb[self.dm_PtVec.Get_List([e])]
            return Kt,Rt
        return K,R
    def solve(self, method, weight, P=1.0, smoothing="", stab=2.0,fictmet="trivial"):
        """
        Solves the deformation of the block matrix given the method name and
        influence function. The influence support is decided at initialization
        of the PeriBlock object.

        The method and weight strings key directly into the peridynamics kernels
        in the husk. There is an additional smoothing argument, which, if not
        false-valued, will call the smoothing kernel after assembling and
        solving the system.

        If there are cut bonds present, self.HCut, the P argument will be used
        to apply the fluid-filled pressure force on the bonds.
        """
        # Assemble the matrix and load for the given peridynamics law
        K,R = self._assemble_KR(method,weight,stab)
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
        # Add the ficticious domain component
        if self.ficticious:
            K,R = self._assemble_KR_fict(K,R, fictmet,method,weight,stab)
        else:
            R-= self.data['load'][0].ravel() * self.data['p_Vol'][0]**0.5/self.data['p_Vol'][0]

        # Apply boundary conditions and then solve the matrix system
        cf.Apply_BC(self.diridofs,self.ubc, K,R)

        u = splin.spsolve(K,R)
        # If we specified a smoothing function, post process the solution
        if smoothing:
            us, = cf.Assemble(hp.__dict__['kernel_smooth_{0}'.format(weight)],
                              self.HAdj, [self.data,{'y':(u,self.dm_PtVec)}],
                              {'ys':(self.dm_PtVec,)},
                              gdim*self.NPart)
            u = us
        # The PeriBlock object is stateless w.r.t. solution, so
        return u

    def output(self, fname, u=None):
        """
        Writes a vtk file with the bond states and, optionally, the displacement
        solution.
        """
        cf.GraphIO.write_graph(fname,self.HPair,self.x,
                       [('x',self.x)] +
                       ([('u',u.reshape((-1,2)))] if not u is None else []),
                       [('alpha', self.data['alpha'][0])]
                       )
