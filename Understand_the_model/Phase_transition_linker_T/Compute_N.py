import numpy as np
import sys
sys.path.append('/home/hcleroy/PostDoc/Simulation/Aging_Condensates/Gillespie_backend/')
import Gillespie_backend as backend
from multiprocessing import Pool

Emin,Emax = -50.,0
Npoints = 100
ell_tot = 50
#distance_anchor = 1000.
rho0 = 6*10**-4
teq = 1000
t_compute = 1000
Nreplica = 100 # number of copy of the gillespie we generate to overcome nonergodicity
reset_linker_time = 10 # number of time per simulation the crosslinkers are reset
filename = 'rho/N_E_L50_rho6_E-4.npy'

def get_N(BindingEnergy,seed):
    S = backend.Gillespie(ell_tot,rho0,BindingEnergy,seed=seed,sliding=False)
    for t in range(teq):
        if t%t_compute//reset_linker_time==0 and t!=0:
            S.reset_crosslinkers()
        S.evolve()
    Nav = 0
    ttot = 0
    #print('process : '+"{:.2f}".format(BindingEnergy)+' equilibration over')
    for t in range(t_compute):
        if t%(t_compute//reset_linker_time)==0:
            #print('process : '+"{:.2f}".format(BindingEnergy)+' time : '+str(t))
            S.reset_crosslinkers()
        Nloop = S.get_N_loop()
        bind,dt = S.evolve()
        if np.isnan(Nloop):
            print('loop')
        if np.isnan(dt):
            print('dt')
        Nav += Nloop*dt
        ttot+=dt
    if ttot!=0:
        return Nav/ttot
    else:
        return 0.

N_E = np.zeros((Npoints,2),dtype=float)
N_E[:,1] = np.linspace(Emax,Emin,Npoints)
for replica in range(Nreplica):
    with Pool(10) as p:
        N_av = p.starmap(get_N,zip(np.linspace(Emax,Emin,Npoints),
                                np.array([np.random.randint(10000000) for _ in range(Npoints)]))
                                ) 
    N_av = np.array(N_av)
    if any(np.isnan(N_av)):
        print(N_av)
        input()
    N_E[:,0]+=N_av/ell_tot/Nreplica
np.save(filename,N_E)