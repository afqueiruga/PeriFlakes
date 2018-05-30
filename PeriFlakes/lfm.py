from SimDataDB import *

from PeriBlock import PeriBlock
import numpy as np
NS = [25]
RFS = [1.5,2.5]
methods = ['Silling','Oterkus2','Fbased','Fstab_Littlewood','Fstab_Silling']
weights = ['const','inv','linear','quadr','cubic','quarticA']
sdb = SimDataDB('lfm_data.db')

onum = 0

P = 1.0
a = 0.1
E = 1.0
nu= 0.25

def analytical(X):
    C = 0.0
    Zhat = lambda z : np.sqrt( -a**2 + z**2.0 ) - z - C *z
    Z = lambda z : z/np.sqrt(-a**2 + z**2.0 ) - 1 - C
    ux = lambda z : P/(2.0*E/(2.0*(1.0+nu))) * (
        (2.0-4.0*nu)/2.0 * np.real(Zhat(z)) - np.imag(z) * np.imag(Z(z)))
    uy = lambda z : P/(2.0*E/(2.0*(1.0+nu))) * (
        (4.0-4.0*nu)/2.0 * np.imag(Zhat(z)) - np.imag(z) * np.real(Z(z)))
    u = np.empty(X.shape)
    u[:,0] = np.sign(X[:,0])*ux( np.abs(X[:,0]) + 1j*np.abs(X[:,1]) )
    u[:,1] = np.sign(X[:,1])*uy( np.abs(X[:,0]) + 1j*np.abs(X[:,1]) )
    return u

@sdb.Decorate('uniaxial',
              [('method','TEXT'),('weight','TEXT'),('RF','FLOAT'),('N','INT')],
              [('u','array')])
def sim(met,wei,RF,N):
    global onum
    print "Solving", met," ",wei," ",RF," ",N

    u=PB.solve(met,wei)
    if False:
        PB.output('./outs/poo_{0}.vtk'.format(onum),u)
        onum += 1
    
    return u,

for N in [24,50,74,100,]:
    for RF in [1.5]:
        PB = PeriBlock(1.0, N, RF*2.0/float(N))
        # Make far-field bcs by evaluating the analytical solution
        farfieldidx = np.array([A
                                for A in [PB.right,PB.left,PB.bottom,PB.top]
                               ],dtype=np.intc).flatten()
        X = PB.x[A,:]
        U = analytical(X)
        PB.setbcs([(farfieldidx,-1)], [])
        PB.ubc = U.ravel()
        PB.cutbonds(-a,0,a,0)
        PB.output('./outs/alpha.vtk')
        # from IPython import embed ; embed()
        for met in methods:
            for wei in weights:
                sim(met,wei,RF,N)
