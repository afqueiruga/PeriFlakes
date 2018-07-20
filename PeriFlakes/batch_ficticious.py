from SimDataDB import *
from PeriBlock import PeriBlock

NS = [10,15,20,25,35,50,75,100]
RFS = [1.5]
methods = ['Silling','Fbased']#,'Fstab_Littlewood','Fstab_Silling']
weights = ['cubic','const']#['const','inv','linear','quadr','cubic','quarticA']
surface_methods = ['none','trivial','bobaru','both']
sdb = SimDataDB('results_ficticious.db')

E = 1.0
nu = 0.25

@sdb.Decorate('uniaxial',
              [('surface','TEXT'),('method','TEXT'),('weight','TEXT'),('RF','FLOAT'),('N','INT')],
              [('x','array'),('u','array')])
def sim(surf,met,wei,RF,N):
    print "Solving", surf, " ", met," ",wei," ",RF," ",N
    global onum
    PB = PeriBlock(1.0, N, RF*2.0/float(N), E,nu, ficticious=surf!="none")
    if surf=="none":
        PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.top,[0,1])])
    else:
        PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.ftop,[0,1])])

    u=PB.solve(met,wei,fictmet = surf)
    return PB.x,u
for N in NS:
    for RF in RFS:
        for met in methods:
            for wei in weights:
                for surf in surface_methods:
                    sim(surf,met,wei,RF,N)
