from PeriBlock import PeriBlock

onum = 0
def sim(met,wei,RF,N):
    global onum
    u=PB.solve(met,wei)
    PB.output('./outs/poo_{0}.vtk'.format(onum),u)
    onum += 1
    return PB.solve(met,wei),

for N in [25]:
    for RF in [1.5]:
        PB = PeriBlock(1.0,25, RF*2.0/25.0)
        PB.setbcs()
        for met in ['silling','Oterkus2','Fbased']:
            for wei in ['const','inv','linear','quadr','cubic','quarticA']:
                sim(met,wei,RF,N)
