import cornflakes as cf
import numpy as np
from PeriBlock import PeriBlock

L = 1.0
N = 10
delta = 1.5*2.0*L/float(N)

PB = PeriBlock(L, N, delta, ficticious=True)
PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.top,1)])

eps = 1.0e-12
inside = lambda x : x[0]<L+eps and x[0]>-L-eps or x[1]>-L-eps or x[1]<L+eps

HStencil4 = cf.Hypergraph()
HStencil3 = cf.Hypergraph()
for e in PB.HFict:
    x0 =  PB.x[e[0],:]
    bodyvertices = []
    print
    print x0
    for i in e[1:(len(e)/2+1)]:
        print PB.x[i,:]
        if x0[0] > L:
            if PB.x[i,0] < L+eps:
                bodyvertices.append(i)
        elif x0[0] <- L:
            if PB.x[i,0] > -L-eps:
                bodyvertices.append(i)
        elif x0[1] > L:
            if PB.x[i,1] < L+eps:
                bodyvertices.append(i)
        elif x0[1] <- L:
            if PB.x[i,1] > -L-eps:
                bodyvertices.append(i)
    bodyvertices.sort( key=lambda i: np.linalg.norm(PB.x[i,:]-x0) )
    if len(bodyvertices)==2:
        HStencil3.Push_Edge([e[0]]+bodyvertices)
    elif len(bodyvertices)==3:
        HStencil4.Push_Edge([e[0]]+bodyvertices)
