from SimDataDB import *

from PeriBlock import PeriBlock
import numpy as np
import timeit

NS = [25]
RFS = [1.5]
stabs = np.linspace(0.1,10.0,10)
methods = ['Silling','Oterkus2','Fbased'] #,'Fstab_Littlewood','Fstab_Silling']
weights = ['cubic'] #['const','inv','linear','quadr','cubic','quarticA']
smoothing_weights = ['const','inv','linear','quadr','cubic','quarticA']
sdb = SimDataDB('lfm_data.db')
smooth_methods = ['Fbased']
stab_methods = ['Fstab_Littlewood','Fstab_Silling']
onum = 0

P = 0.1
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

@sdb.Decorate('lfm',
              [('method','TEXT'),('weight','TEXT'),('smoothing','TEXT'),
               ('stab','FLOAT'),('RF','FLOAT'),('N','INT')],
              [("error","FLOAT"),("time","FLOAT"),('u','array')])
def sim(met,wei,smo,stab,RF,N):
    global onum
    print "Solving", met," ",wei," ",smo," ",stab," ",RF," ",N
    def work():
        work.u=PB.solve(met,wei, P=P, smoothing=smo, stab=stab)
    runtime = timeit.timeit(work,number=1)
    if False:
        PB.output('./outs/poo_{0}.vtk'.format(onum),work.u)
        onum += 1
    ua = analytical(PB.x)
    en=work.u-ua.ravel()
    error = np.linalg.norm(en) / np.linalg.norm(ua.flatten())
    return error,runtime,work.u,

for N in [24,38,50,62,74,86,100,112,124,150,162,174,186,200,212,224,238,250]:
    for RF in [1.5]:
        PB = PeriBlock(1.0, N, RF*2.0/float(N), E=E,nu=nu)
        # Make far-field bcs by evaluating the analytical solution
        farfieldidx = np.array([A
                                for A in [PB.right,PB.left,PB.bottom,PB.top]
                               ],dtype=np.intc).flatten()
        X = PB.x[farfieldidx,:]
        U = analytical(X)
        # from IPython import embed ; embed()
        PB.setbcs([(farfieldidx,-1)], [])
        PB.ubc = U.ravel()
        PB.cutbonds(-a,0,a,0)
        PB.output('./outs/alpha.vtk')

        # The big grid search
        for met in methods:
            for wei in weights:
                sim(met,wei,"",0.0,RF,N)
        # Only try smoothing methods on one method
        for met in smooth_methods:
            for smo in smoothing_weights:
                sim(met,"cubic",smo,0.0,RF,N)
        # And only two of the methods use the stabilizer parameter
        for met in stab_methods:
            for s in stabs:
                sim(met,"cubic","",s,RF,N)
