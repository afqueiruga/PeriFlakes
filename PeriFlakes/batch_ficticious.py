from SimDataDB import *
from PeriBlock import PeriBlock

NS = [10,15,20,25,35,50,75,100,125,150,175,200]
RFS = [1.5]
methods = ['Silling','Fbased','Fstab_Silling','Oterkus2']#,'Fstab_Littlewood','Fstab_Silling']
weights = ['cubic','const']#['const','inv','linear','quadr','cubic','quarticA']
surface_methods = ['none','trivial','bobaru','both']
sdb = SimDataDB('results_ficticious_nu0.db')

E = 1.0
nu = 0.0

@sdb.Decorate('uniaxial',
              [('surface','TEXT'),('method','TEXT'),('weight','TEXT'),('RF','FLOAT'),('N','INT')],
              [('x','array'),('u','array')])
def sim_uni(surf,met,wei,RF,N):
    print "Solving", surf, " ", met," ",wei," ",RF," ",N
    global onum
    PB = PeriBlock(1.0, N, RF*2.0/float(N), E,nu, ficticious=surf!="none")
    if surf=="none":
        PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.top,[0,1])])
    else:
        PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.ftop,[0,1])])

    u=PB.solve(met,wei,fictmet = surf)
    return PB.x,u

@sdb.Decorate('shear',
              [('surface','TEXT'),('method','TEXT'),('weight','TEXT'),('RF','FLOAT'),('N','INT')],
              [('x','array'),('u','array')])
def sim_shear(surf,met,wei,RF,N):
    print "Solving", surf, " ", met," ",wei," ",RF," ",N
    global onum
    PB = PeriBlock(1.0, N, RF*2.0/float(N), E,nu, ficticious=surf!="none")
    if surf=="none":
        PB.setbcs([(PB.right,1),(PB.left,1),(PB.bottom,0)], [(PB.top,[1,0])])
    else:
        PB.setbcs([(PB.right,1),(PB.left,1),(PB.bottom,0)], [(PB.ftop,[1,0])])

    u=PB.solve(met,wei,fictmet = surf)
    return PB.x,u


@sdb.Decorate('isotropic',
              [('surface','TEXT'),('method','TEXT'),('weight','TEXT'),('RF','FLOAT'),('N','INT')],
              [('x','array'),('u','array')])
def sim_iso(surf,met,wei,RF,N):
    print "Solving", surf, " ", met," ",wei," ",RF," ",N
    global onum
    PB = PeriBlock(1.0, N, RF*2.0/float(N), E,nu, ficticious=surf!="none")
    if surf=="none":
        PB.setbcs([(PB.left,0),(PB.bottom,1)], [(PB.top,[0,1]),(PB.right,[1,0])])
    else:
        PB.setbcs([(PB.left,0),(PB.bottom,1)], [(PB.ftop,[0,1]),(PB.fright,[1,0])])

    u=PB.solve(met,wei,fictmet = surf)
    return PB.x,u


for N in NS:
    for RF in RFS:
        for met in methods:
            for wei in weights:
                for surf in surface_methods:
                    for sim in [sim_uni,sim_shear,sim_iso]:
                        sim(surf,met,wei,RF,N)
