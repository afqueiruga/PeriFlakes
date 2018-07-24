from SimDataDB import *
import numpy as np
from PeriBlock import PeriBlock

NS = [10,20,25,35,50] #,60,70,80,90,100]
RFS = [1.5,2.5,3.5] #,2.5,3.5] #,2.5,3.01,3.5]
methods = ['Fbased','Silling'] #,'Oterkus2','Fbased']#,'Fstab_Littlewood','Fstab_Silling']
weights = ['cubic']#['const','inv','linear','quadr','cubic','quarticA']
surface_methods = ['none','trivial',]
sdb = SimDataDB('results_dirr.db')

onum = 0

E = 1.0
nu = 0.25
T = 1.0
H = 2.0
sol = lambda x : ((x[:,1]+H/2.0)* \
                 np.array([[ 0, \
                    T*(2.0*nu**2+nu-1.0)/(E*(nu-1)) ]]).T).T

@sdb.Decorate('uniaxial',
              [('surface','TEXT'),('method','TEXT'),('weight','TEXT'),('RF','FLOAT'),('N','INT')],
              [('x','array'),('u','array')])
def sim_uni(surf,met,wei,RF,N):
    print "Solving", surf, " ", met," ",wei," ",RF," ",N
    global onum
    PB = PeriBlock(1.0, N, RF*2.0/float(N), E,nu, ficticious=surf!="none")
    if surf=="none":
        PB.setbcs([(PB.right,-1),(PB.left,-1),(PB.bottom,-1),(PB.top,-1)])
    else:
        PB.setbcs([(PB.right,-1),(PB.left,-1),(PB.bottom,-1),(PB.top,-1)])
    # Set all of them and see if it "acccepts" the solution
    PB.data['y'][0][:,:] += sol(PB.x)
    # from IPython import embed; embed()
    Dy=PB.solve(met,wei,fictmet = surf)
    ysolved = PB.data['y'][0].ravel() + Dy
    u = ysolved - PB.x.ravel()
    return PB.x,u

@sdb.Decorate('uniaxial_try_harder',
              [('surface','TEXT'),('method','TEXT'),('weight','TEXT'),('RF','FLOAT'),('N','INT')],
              [('x','array'),('u','array')])
def sim_uni(surf,met,wei,RF,N):
    print "Solving", surf, " ", met," ",wei," ",RF," ",N
    global onum
    PB = PeriBlock(1.0, N, RF*2.0/float(N), E,nu, ficticious=surf!="none")
    if surf=="none":
        PB.setbcs([(PB.right,-1),(PB.left,-1),(PB.bottom,-1),(PB.top,-1)])
    else:
        PB.setbcs([(PB.right,-1),(PB.left,-1),(PB.bottom,-1),(PB.top,-1)])
    # Only set the edges
    PB.data['y'][0].ravel()[PB.diridofs] += sol(PB.x).ravel()[PB.diridofs]
    # from IPython import embed; embed()
    Dy=PB.solve(met,wei,fictmet = surf)
    ysolved = PB.data['y'][0].ravel() - Dy
    u = ysolved - PB.x.ravel()
    return PB.x,u

for N in NS:
    for RF in RFS:
        for met in methods:
            for wei in weights:
                for surf in surface_methods:
                    sim_uni(surf,met,wei,RF,N)
