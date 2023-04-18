import numpy as np
import sys
from multiprocessing import Pool
sys.path.append('../Module_Analysis/')
import Pair_Corelation_Function as PCF
#import Clustering as Clust
sys.path.append('/home/hcleroy/PostDoc/aging_condensates/Simulation/Gillespie/Gillespie_backend')
import Gillespie_backend as backend
import Simulate_System as SimSys

# define the unit of the system :
dimension = 3

L = 1000
NLinker = 20

EbLow = -1.
EbHigh = -10.

kdiff = 0.000207352*5

seed = np.random.randint(0,1000000)

def Make_Simulation(a0,a1,a2,a3,a4,a5,a6,step_tot):
    Sys = backend.Gillespie(ell_tot=a0,
                                    rho0=a1,
                                    BindingEnergy=a2,
                                    kdiff=a3,
                                    seed = a4,
                                    Nlinker=a5,dimension=a6)
    Sim = SimSys.Simulation(step_tot = step_tot,Gillespie = Sys)
    Sim.simulate_eq_distribution()
    return Sim.R,Sim.L

#args = [(L,0.,eb,kdiff,np.random.randint(0,100000),NLinker,dimension,10**6) for eb in np.linspace(-10,0,10,endpoint=False)]
args = (L,0.,-9,kdiff,np.random.randint(0,100000),NLinker,dimension,10**6)

#pool = Pool(12)
#Sims = pool.starmap(Make_Simulation,args)

Rs,Ls = Make_Simulation(*args)

#Rs = np.array([Sims[i][0] for i in range(10)])
#Ls = np.array([Sims[i][1] for i in range(10)])

np.save('step_10_E6_E_10_1_Rs.npy',Rs,allow_pickle=True)
np.save('step_10_E6_E_10_1_Ls.npy',Ls,allow_pickle=True)