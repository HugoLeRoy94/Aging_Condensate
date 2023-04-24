import numpy as np
import sys

# import the C++ module
sys.path.append('/home/hcleroy/PostDoc/aging_condensates/Simulation/Gillespie/Gillespie_backend')
import Gillespie_backend as backend
import Simulate_Gillespie as SimSys

R = 50
L = 100
Eb = -10.
kdiff = 0.000207352
size = 100
dangling = True

Sys =  backend.Gillespie(ell_tot=L,rho0=0.,BindingEnergy=sys.argv[1],kdiff=kdiff*5,seed = np.random.randint(0,1000000),Nlinker=2)
Sim = SimSys.Simulation(step_tot = 10**5,size=100,Gillespie = Sys)

