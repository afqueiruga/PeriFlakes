from PeriBlock import PeriBlock

PB = PeriBlock(1.0,25, 1.5*2.0/25.0)

PB.setbcs()

u=PB.solve('silling','const')

PB.output('./poo.vtk',u)
