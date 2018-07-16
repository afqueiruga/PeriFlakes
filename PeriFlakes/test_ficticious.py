from PeriBlock import PeriBlock

import unittest as ut


class construct(ut.TestCase):
    def onecase(self,L,N,delta):
        PB = PeriBlock(L, N, delta, ficticious=True)
        PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.top,1)])

        # Assert all nodes have a computation
        self.assertTrue(len(PB.HAdj)+len(PB.HFict) == PB.x.shape[0])

        # Assert the ficticious graph is only the exterior
        for e in PB.HFict:
            xi = PB.x[e[0],:]
            self.assertTrue(xi[0] > L or xi[0] < -L or xi[1] < -L or xi[1] > L)

        # Assert that there's enough buffer
        self.assertTrue(PB.x[:,0].max() >  L + delta)
        self.assertTrue(PB.x[:,0].min() < -L - delta)
        self.assertTrue(PB.x[:,1].max() >  L + delta)
        self.assertTrue(PB.x[:,1].min() < -L - delta)

        # Verify some of the nodes are exactly on the boundary
        
    def test(self):
        self.onecase(1.0,10,0.5)
        self.onecase(1.0,10,0.25)
        self.onecase(2.0,10,0.5)
        self.onecase(1.0,20,0.15)


#K,R = PB._assemble_KR('Fbased','cubic')

    
if __name__=='__main__':
    ut.main()
