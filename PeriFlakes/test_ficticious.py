from PeriBlock import PeriBlock

PB = PeriBlock(1.0, 10, 0.5, ficticious=True)

PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.top,1)])

K,R = PB._assemble_KR('Fbased','cubic')
