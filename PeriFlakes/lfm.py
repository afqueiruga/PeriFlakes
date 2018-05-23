from SimDataDB import *

from PeriBlock import PeriBlock

NS = [25]
RFS = [1.5,2.5]
methods = ['Silling','Oterkus2','Fbased','Fstab_Littlewood','Fstab_Silling']
weights = ['const','inv','linear','quadr','cubic','quarticA']
sdb = SimDataDB('lfm_data.db')

onum = 0

a = 0.1

@sdb.Decorate('uniaxial',
              [('method','TEXT'),('weight','TEXT'),('RF','FLOAT'),('N','INT')],
              [('u','array')])
def sim(met,wei,RF,N):
    print "Solving", met," ",wei," ",RF," ",N
    global onum
    u=PB.solve(met,wei)
    PB.output('./outs/poo_{0}.vtk'.format(onum),u)
    onum += 1
    return PB.solve(met,wei),
for N in [25]:
    for RF in [1.5]:
        PB = PeriBlock(1.0,25, RF*2.0/float(N))
        PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.top,1)])
        PB.cutbonds(-a,0,a,0)
        for met in methods:
            for wei in weights:
                sim(met,wei,RF,N)
