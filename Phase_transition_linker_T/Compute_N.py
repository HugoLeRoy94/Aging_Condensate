import numpy as np
import sys
sys.path.append('/home/hcleroy/PostDoc/Simulation/Aging_Condensates/System_backend/')
import System_backend as backend
from multiprocessing import Pool

Emin,Emax = -10,0.1

ell_tot = 1000.
#distance_anchor = 1000.
rho0 = 5.*10**-4
teq = 10000
t_compute = 1000


def get_N(BindingEnergy):
    S = backend.System(ell_tot,rho0,BindingEnergy,seed=np.random.randint(100000))
    for t in range(teq):
        if t%t_compute//10==0 and t!=0:
            S.reset_crosslinkers()
        S.evolve()
        Nav = 0
    ttot = 0
    for t in range(t_compute):
        Nloop = S.get_N_loop()
        bind,dt = S.evolve()
        Nav += Nloop*dt
        ttot+=dt
    return Nav/ttot

with Pool(12) as p:
    N_av = p.map(get_N,np.linspace(Emin,Emax,10)) 

N_E = np.array([N_av,np.linspace(Emin,Emax,10)]).transpose()
np.save('N_E_rho_1_25E-4.npy',N_E)