from popcorn import *

#
# This is a kernel for the ficticious node method.
# It is more like a finite difference stencil than the PD point
# cloud, so we generate it in a different file.
#
# The hyperedge layout is:
#  Pf  P1 P2 P3
#
# to represent this:
#
#       (Pf)
#
#
#  (P1) (P2) (P3)
#

gdim = 2
PointSca = DofSpace(1, 0, 4,1)
PointVec = DofSpace(gdim, 0,4,1)
Param = DofSpace(1, 0,1)
GlobalVec = DofSpace(gdim, -1)

i_y = Input("y",PointVec)
i_x = Input("x",PointVec)
i_E = Input("p_E", Param)
i_nu = Input("p_nu", Param)
i_t = Input("load", GlobalVec)

o_R = Output("R", [PointVec],1)
o_K = Output("K", [PointVec],2)

yf, yp, yo, ym = i_y.Vertex_Split()
xf, xp, xo, xm = i_x.Vertex_Split()

uf,up,uo,um = yf-xf, yp-xp, yo-xo, ym-xm

# Let's make the axis-aligned kernel

Dy = xf[1] - xo[1]
Dx = xp[0] - xm[0]

rowx = (uf[0]-uo[0])/Dy + i_nu[0]*(up[1]-um[1])/(2*Dx) - i_t[0]*(1-i_nu[0]**2)/i_E[0]
rowy = (uf[1]-uo[1])/Dy + i_nu[0]*(up[0]-um[0])/(2*Dx) - i_t[1]*(1+i_nu[0])/i_E[0]


R_expr = Matrix([ rowx, rowy, 0,0, 0,0, 0,0])
K_expr = R_expr.jacobian(i_y)

Kernel("ficticious_bobaru_y",
       listing=[
           Asgn(o_R,R_expr,'='),
           Asgn(o_K,K_expr,'=')
       ])



# Now let's generalize it

def norm(a):
    return sqrt((a.T*a)[0,0])
def difference(ua,ub, xa,xb):
    return (ua-ub).dot(xa-xb) / (xa-xb).dot(xa-xb)

en = difference(uf,uo, xf,xo)
nn = (xf-xo)/norm(xf-xo)

et_plus = difference(up,uo, xp,xo)
et_minus = difference(uo,um, xo,xm)
nt_plus = (xp-xo)/norm(xp-xo)
nt_minus = (xo-xm)/norm(xo-xm)
et = (et_plus+et_minus)/2
nt = (nt_plus+nt_minus)/2

KnuE = 1/i_E[0] * Matrix([[ 1+i_nu[0],    1-i_nu[0]**2],
                          [ 1-i_nu[0]**2, 1+i_nu[0]   ]])
constraint = en*nn + i_nu[0]*et*nt - i_t.multiply_elementwise( KnuE * nn )

R_expr = Matrix([ constraint[0], constraint[1], 0,0, 0,0, 0,0])
K_expr = R_expr.jacobian(i_y)

Kernel("ficticious_bobaru_n",
       listing=[
           Asgn(o_R,R_expr,'='),
           Asgn(o_K,K_expr,'=')
       ])


# Write out the kernels
Husk('ficticious')
