from PeriBlock import PeriBlock

import unittest as ut


class construct(ut.TestCase):
    def onecase(self,L):
        PB = PeriBlock(L, 10, 0.5, ficticious=True)
        PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.top,1)])

        self.assertTrue(len(PB.HAdj)+len(PB.HFict) == PB.x.shape[0])
    
        for e in PB.HFict:
            xi = PB.x[e[0],:]
            self.assertTrue(xi[0] > L or xi[0] < -L or xi[1] < -L or xi[1] > L)
        K,R = PB._assemble_KR('Fbased','cubic')

    def test(self):
        self.onecase(1.0)
        self.onecase(2.0)

    
if __name__=='__main__':
    ut.main()
