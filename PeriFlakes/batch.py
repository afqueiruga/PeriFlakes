from SimDataDB import *

from PeriBlock import PeriBlock

onum = 0

sdb = SimDataDB('results.db')
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

NS = [25]
RFS = [1.5,2.5]
methods = ['silling','Oterkus2','Fbased','Fstab_Littlewood','Fstab_Silling']
weights = ['const','inv','linear','quadr','cubic','quarticA']

for N in [25]:
    for RF in [1.5]:
        PB = PeriBlock(1.0,25, RF*2.0/25.0)
        PB.setbcs()
        for met in methods:
            for wei in weights:
                sim(met,wei,RF,N)
