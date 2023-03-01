import numpy as np
import sys
sys.path.append('/home/hcleroy/PostDoc/Simulation/Aging_Condensate/')
import Gillespie_backend as backend


temperature = 0.15
ell_tot = 2000.
#distance_anchor = 1000.
rho00 = 0.001
teq = 100
t_compute = 10000

ell_rho = list()
N_rho = list()

for rho0 in np.linspace(rho00,100*rho00,100):
    print(rho0)
    S = backend.Gillespie(ell_tot,rho0,temperature)
    for t in range(teq):
        S.evolve()
    Nav = 0
    lav = 0
    for t in range(t_compute):
        S.evolve()
        Nav += S.get_N_loop()/(t_compute-1)
        #lav += np.mean(np.array([S.get_ell()[i+1]-S.get_ell()[i]  for i in range(S.get_ell().shape[0]-1)]))/t_compute
        lav +=np.mean(S.get_ell())/(t_compute-1)
    ell_rho.append([lav,rho0])
    N_rho.append([Nav,rho0])
np.save('ell_rho.npy',ell_rho)
np.save('N_rho.npy',N_rho)
