import numpy as np
import sys
sys.path.append('/home/hcleroy/PostDoc/Simulation/Aging_Condensate')
import System_backend as backend

Tmin =0.001
Tmax =0.25

ell_tot = 2000.
distance_anchor = 1000.
rho0 = 0.0001
teq = 100
t_compute = 10000

N_T = list()

for Temp in np.linspace(Tmin,Tmax,100):
    print(Temp)
    S = backend.System(ell_tot,distance_anchor,rho0,Temp)
    for t in range(teq):
        S.evolve()
        Nav = 0
    for t in range(t_compute):
        S.evolve()
        Nav += S.get_N_loop()/t_compute
    N_T.append([Temp,Nav])

N_T = np.array(N_T)
np.save('N_T_rho_0_0001.npy',N_T)