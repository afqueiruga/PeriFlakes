from popcorn import *
import operator
# This is going to be a variable length kernel, but it's a little funky.
# It has both Particle and Bond DOFs, and in the edge they'll be ordered like this
#
# P0  P1 P2 P3 ... Pn  B1 B2 B3 ... Bn 
#
# to represent this
#
# (P2)~~B2~~ ((P0)) ~~B1~~(P1)
#                `.
#                  `~B3~~(P3)
#
gdim = 2
II = Symbol('i')
JJ = Symbol('j')
l_edge = Symbol("l_edge")
Npart = (l_edge+1)/2

# Typedefs
PointVec = DofSpace(gdim, 0,Npart,1)
BondSca = DofSpace(1, Npart, l_edge)
PointSca = DofSpace(1, 0, Npart)
PointTen = DofSpace(gdim*gdim, 0,Npart)
Param = DofSpace(1, 0,1)

# Inputs
i_y = Input("y",PointVec)
i_x = Input("x",PointVec)
i_p = Input("p",PointSca)
i_r = Input("r",PointSca)
i_alpha = Input("alpha",BondSca)
i_pfrac= Input("pfrac",BondSca)
# Parameter inputs
i_delta = Input("delta", Param)
i_E = Input("p_E", Param)
i_nu = Input("p_nu", Param)
i_Vol = Input("p_Vol", Param)
# Outputs
o_R  = Output("R",  [PointVec], 1)
o_K  = Output("K",  [PointVec], 2 )

# Do some symbolic math
c_K = i_E[0]/(3.0-6.0*i_nu[0])
c_G = i_E[0]/(2.0+2.0*i_nu[0])
y0,yI,yJ = i_y.Vertex_Handles(0,II,JJ)
x0,xI,xJ = i_x.Vertex_Handles(0,II,JJ)
alpha_I = i_alpha.Vertex_Handle(II-1)[0,0]
pfrac_I = i_pfrac.Vertex_Handle(II-1)[0,0]
alpha_J = i_alpha.Vertex_Handle(JJ-1)[0,0]
pfrac_J = i_pfrac.Vertex_Handle(JJ-1)[0,0]

def norm(a):
    return sqrt((a.T*a)[0,0])
rxI = xI-x0
ryI = yI-y0
rxabs = norm(rxI)
ryabs = norm(ryI)
xn = rxI/rxabs
yn = ryI/ryabs
strainI = (ryabs-rxabs)

m     = PopcornVariable("m",1,1)
theta = PopcornVariable("theta",1,1)

class StateBased():
    def __init__(self, name, w_of_r):
        self.name = name
        self.w_of_r = w_of_r
        self.w0I = w_of_r( rxabs )
    def m_expr(self):
        return self.w0I* alpha_I * rxabs*rxabs
    def theta_expr(self):
        return gdim/m[0] *  self.w0I*alpha_I * ((ryabs-rxabs)* rxabs)
    def force(self):
        A =  alpha_I * self.w0I/m[0] *(  ( 2.0*(c_K-c_G/3.0) * theta[0] )* rxabs
                                          + 8.0*c_G * (strainI - theta[0]*rxabs/3.0) )
        return A * xn
    def _sum_prgm(self, var, expr):
        return [
            var,
            Loop(II,1,Npart,[
                Asgn(var, expr, "+=")
                ])
            ]
    def _tang_prgm(self, var, expr):
        return []
    def list_states(self):
        return [(theta,self.theta_expr())]
    def list_moments(self):
        return [(m, Matrix([self.m_expr()]) )]
    def kernel(self):
        theta_expr = self.theta_expr()
        tA = self.force()
        theta_Dy0 = PopcornVariable("theta_Dy0",gdim,1)
        theta_exprJ = theta_expr.subs(
            zip(yI,yJ)+zip(xI,xJ)+[(alpha_I,alpha_J)])
        
        theta_DyJ = Matrix([theta_exprJ]).jacobian(yJ).T
        # And tangent matrices
        tA_Dtheta = tA.diff(theta[0])
        tAI_Dy0 =  tA.jacobian(y0) + tA_Dtheta*theta_Dy0.T 
        tAI_DyI =  tA.jacobian(yI)
        tAI_DyJ = tA_Dtheta * theta_DyJ.T
        prgm = reduce(operator.add,[self._sum_prgm(P,E) for P,E in self.list_moments()]) + \
            self._sum_prgm(theta, Matrix([theta_expr])) + \
            self._sum_prgm( theta_Dy0, Matrix([theta_expr]).jacobian(y0).T ) + \
        [
            Loop(II,1,Npart,[
                Asgn(o_R.View((0,)), tA,"+="),
                Asgn(o_R.View((gdim*II,)), tA,"-="),                
            ]),
            Loop(II,1,Npart,[
                Asgn(o_K.View((0,0)),             tAI_Dy0,"+="),
                Asgn(o_K.View((0,gdim*II)),       tAI_DyI,"+="),
                Asgn(o_K.View((gdim*II,0)),       tAI_Dy0,"-="),
                Asgn(o_K.View((gdim*II,gdim*II)), tAI_DyI,"-="),
                Loop(JJ,1,Npart,[
                    Asgn(o_K.View((0,gdim*JJ)),   tAI_DyJ,"+="),
                    Asgn(o_K.View((gdim*II,gdim*JJ)),   tAI_DyJ,"-="),
                ])
            ])
        ]
        return Kernel(self.name, listing=prgm)

a_integ = PopcornVariable('a_integ',1,1)
b_denom = PopcornVariable('b_denom',1,1)
class Oterkus(StateBased):
    def m_expr(self):
        return self.w0I * alpha_I * rxabs
    def a_expr(self):
        return self.w0I * alpha_I * rxabs**2
    def b_expr(self):
        return self.w0I * alpha_I * (xn[0] * xn[1] )**2 * rxabs**2
    def list_moments(self):
        return [(m, Matrix([self.m_expr()]) ),
                (a_integ,Matrix([self.a_expr()])),
                (b_denom,Matrix([self.b_expr()])) ]
    def theta_expr(self):
        return 2/m[0] * self.w0I * alpha_I * (ryabs-rxabs) *  (yn.T*xn)[0,0]
    def force(self):
        a = (c_K - 2*c_G*a_integ[0]/b_denom[0])/2
        b = c_G /( 2* b_denom[0])
        A = 2 * self.w0I* alpha_I  * ( a*2/m[0]*(yn.T*xn)[0,0]*theta[0] + b*(ryabs-rxabs) )
        return A * yn
    
mK = PopcornVariable("mK",gdim,2)
mKi = PopcornVariable("mKi",gdim,2)
mN = PopcornVariable("mN",gdim,2)
class Fbased(StateBased):
    def K_expr(self):
        return self.w0I * alpha_I * rxI * rxI.T * i_Vol[0]
    def N_expr(self):
        return self.w0I * alpha_I * ryI * rxI.T * i_Vol[0]
    def stress(self):
        F = mN * mKi
        E = 0.5 * ( F.T + F ) - eye(gdim)
        return ( (c_K-2.0*c_G/3.0)*trace(E)*eye(gdim) + 2.0*c_G*E )
    def force(self):
        P = self.stress()
        return self.w0I * alpha_I * P * mKi.T * rxI * i_Vol[0]
    def kernel(self):
        N_expr = self.N_expr()
        NJ_expr = N_expr.subs(zip(yI,yJ)+zip(xI,xJ)+[(alpha_I,alpha_J)])
        N_dy0 =  [ [PopcornVariable("N00_dy0",gdim,1),
                    PopcornVariable("N01_dy0",gdim,1)],
                   [PopcornVariable("N10_dy0",gdim,1),
                    PopcornVariable("N11_dy0",gdim,1)] ]
        N_dy0_expr = [ [ Matrix([N_expr[j,i]]).jacobian(y0).T
                         for j in xrange(gdim) ]
                       for i in xrange(gdim) ]
        N_dyJ_expr = [ [ Matrix([NJ_expr[j,i]]).jacobian(yJ).T
                         for j in xrange(gdim) ]
                       for i in xrange(gdim) ]
        # Derivatives of the force
        tI = self.force()
        tI_DN = [ [ tI.diff(mN[i,j]) for j in xrange(gdim) ]
                  for i in xrange(gdim) ]
        #from IPython import embed ; embed()
        tI_Dy0 = tI.jacobian(y0) + reduce(operator.add,
                [ tI_DN[i][j]*N_dy0[i][j].T
                  for i in xrange(gdim) for j in xrange(gdim) ])
        tI_DyI = tI.jacobian(yI)
        tI_DyJ = tI.jacobian(yI)+ reduce(operator.add,
                [ tI_DN[i][j]*N_dyJ_expr[i][j].T
                  for i in xrange(gdim) for j in xrange(gdim) ])
        # The final program
        prgm = self._sum_prgm(mK,self.K_expr()) + \
               [ mKi, Asgn(mKi,mK.inv()) ] + \
               self._sum_prgm(mN, N_expr) + \
               reduce(operator.add,[
                   self._sum_prgm(N_dy0[i][j], N_dy0_expr[i][j])
                   for i in xrange(gdim) for j in xrange(gdim)
                   ]) + \
        [
            Loop(II,1,Npart,[
                Asgn(o_R.View((0,)), tI,"+="),
                Asgn(o_R.View((gdim*II,)), tI,"-="),                
            ]),
            Loop(II,1,Npart,[
                Asgn(o_K.View((0,0)),             tI_Dy0,"+="),
                Asgn(o_K.View((0,gdim*II)),       tI_DyI,"+="),
                Asgn(o_K.View((gdim*II,0)),       tI_Dy0,"-="),
                Asgn(o_K.View((gdim*II,gdim*II)), tI_DyI,"-="),
                Loop(JJ,1,Npart,[
                    Asgn(o_K.View((0,gdim*JJ)),   tI_DyJ,"+="),
                    Asgn(o_K.View((gdim*II,gdim*JJ)),   tI_DyJ,"-="),
                ])
            ])
        ]
        return Kernel(self.name, listing=prgm)

        
formulations = {
    'silling' : StateBased,
    'Oterkus2': Oterkus,
    'Fbased':Fbased
}
delta = i_delta[0]
weight_funcs = {
    'const' : lambda r : Piecewise( (1.0    , Le(r,delta)), (0.0,True)),
    'inv'   : lambda r : Piecewise( (delta/r, Le(r,delta)), (0.0,True)),
    'linear': lambda r : Piecewise( ((1.0-r / delta)**1.0, Le(r,delta)), (0.0,True)),
    'quadr' : lambda r : Piecewise( ((1.0-r / delta)**2.0, Le(r,delta)), (0.0,True)),
    'cubic' : lambda r : Piecewise( ((1.0-r / delta)**3.0, Le(r,delta)), (0.0,True)),
    #'cubicspline':lambda r : Piecewise(
    #    (2.0/3.0 - 4.0*r**2 + 4.0* r**3, Le(r,0.5)  ),
    #    (4.0/3.0 - 4.0*r + 4.0*r**2 - 4.0/3.0*r**3, Lt(r,1.0)),
    #    (0.0,True) ),
    'quarticA':lambda r: Piecewise( (1.0 - 6.0*(r/delta)**2 + 8.0*(r/delta)**3 - 3.0*(r/delta)**4, Le(r,delta)), (0.0,True))

}

for n,c in formulations.iteritems():
    for w,f in weight_funcs.iteritems():
       c(n+"_"+w,weight_funcs['cubic']).kernel()


Husk('peridynamics')
