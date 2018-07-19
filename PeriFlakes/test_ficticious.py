

import unittest as ut


class construct(ut.TestCase):
    def onecase(self,L,N,delta):
        from PeriBlock import PeriBlock
        PB = PeriBlock(L, N, delta, ficticious=True)
        PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.top,1)])

        # Assert all nodes have a computation
        self.assertTrue(len(PB.HAdj)+len(PB.HFict) == PB.x.shape[0])

        # All of the nodes on the interior should have the same horizon size
        self.assertTrue(len(PB.HAdj.view())==1)

        # Assert the ficticious graph is only the exterior
        for e in PB.HFict:
            xi = PB.x[e[0],:]
            self.assertTrue(xi[0] > L or xi[0] < -L or xi[1] < -L or xi[1] > L)

        # Assert that there's enough buffer
        h = 2.0*L/float(N)
        self.assertTrue(PB.x[:,0].max() >  L + delta - h)
        self.assertTrue(PB.x[:,0].min() < -L - delta + h)
        self.assertTrue(PB.x[:,1].max() >  L + delta + h)
        self.assertTrue(PB.x[:,1].min() < -L - delta - h)

        # Verify some of the nodes are exactly on the boundary

    def test(self):
        self.onecase(1.0,10,0.5)
        self.onecase(1.0,10,0.25)
        self.onecase(2.0,10,0.5)
        self.onecase(1.0,20,0.15)

class assemble(ut.TestCase):
    " Just make sure it doesn't crash'"
    def onecase(self,L,N,delta):
        from PeriBlock import PeriBlock
        PB = PeriBlock(L, N, delta, ficticious=True)
        PB.setbcs([(PB.right,0),(PB.left,0),(PB.bottom,1)], [(PB.top,1)])
        K,R = PB._assemble_KR_fict('standard','Fbased','cubic')
    def test(self):
        self.onecase(1.0,10,0.5)
        self.onecase(1.0,10,0.25)
        self.onecase(2.0,10,0.5)
        self.onecase(1.0,20,0.15)


class stencil(ut.TestCase):
    def test_y_aligned(self):
        import numpy as np
        import husk_ficticious
        import cornflakes as cf
        ke = husk_ficticious.kernel_ficticious_bobaru_y
        # No load
        dat = np.array([0,0, 1,0.25,
                        0,1, -1,0, 0,0, 1,0 ,
                        0,1, -1,0, 0,0, 1,0 ])
        out = cf.cornflakes_library.call_kernel(ke,4,dat)
        R = out[0:8]
        K = out[8:72].reshape(8,8)
        self.assertTrue(R.sum() == 0.0)

        print "Normal Load:"
        dat = np.array([0,1, 1,0.25,
                        0,1, -1,0, 0,0, 1,0 ,
                        0,1, -1,0, 0,0, 1,0 ])
        out = cf.cornflakes_library.call_kernel(ke,4,dat)
        K = out[0:64].reshape(8,8)
        R = out[64:72]
        print R
        print K

        print "Shear Load:"
        dat = np.array([1,0, 1,0.25,
                        0,1, -1,0, 0,0, 1,0 ,
                        0,1, -1,0, 0,0, 1,0 ])
        out = cf.cornflakes_library.call_kernel(ke,4,dat)
        K = out[0:64].reshape(8,8)
        R = out[64:72]

        print R
        print K
    def test_n_version(self):
        import numpy as np
        import husk_ficticious
        import cornflakes as cf
        ken = husk_ficticious.kernel_ficticious_bobaru_n
        key = husk_ficticious.kernel_ficticious_bobaru_n
        # No load
        dat = np.array([0,0, 1,0.25,
                        0,1, -1,0, 0,0, 1,0 ,
                        0,1, -1,0, 0,0, 1,0 ])
        out = cf.cornflakes_library.call_kernel(ken,4,dat)
        R = out[0:8]
        K = out[8:72].reshape(8,8)
        self.assertTrue(R.sum() == 0.0)

        print "Normal Load:"
        dat = np.array([0,1, 1,0.25,
                        0,1, -1,0, 0,0, 1,0 ,
                        0,1, -1,0, 0,0, 1,0 ])
        outn = cf.cornflakes_library.call_kernel(ken,4,dat)
        outy = cf.cornflakes_library.call_kernel(key,4,dat)
        self.assertTrue( np.linalg.norm(outn-outy) < 1.0e-12 )
        K = outn[0:64].reshape(8,8)
        R = outn[64:72]
        print R
        print K

        print "Shear Load:"
        dat = np.array([1,0, 1,0.25,
                        0,1, -1,0, 0,0, 1,0 ,
                        0,1, -1,0, 0,0, 1,0 ])
        outn = cf.cornflakes_library.call_kernel(ken,4,dat)
        outy = cf.cornflakes_library.call_kernel(key,4,dat)
        self.assertTrue( np.linalg.norm(outn-outy) < 1.0e-12 )
        K = outn[0:64].reshape(8,8)
        R = outn[64:72]
        print R
        print K

        print "Shear Load along x:"
        dat = np.array([0,1, 1,0.25,
                        1,0, 0,-1, 0,0, 0,1 ,
                        1,0, 0,-1, 0,0, 0,1 ])
        outn = cf.cornflakes_library.call_kernel(ken,4,dat)
        K = outn[0:64].reshape(8,8)
        R = outn[64:72]
        print R
        print K

#K,R = PB._assemble_KR('Fbased','cubic')


if __name__=='__main__':
    ut.main()
