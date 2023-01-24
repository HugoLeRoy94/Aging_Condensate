import numpy as np
import sys
#sys.path.append('/home/hcleroy/PostDoc/aging_condensates/Simulation/Gillespie/System_backend')
#import System_backend as Sys
class Simulation:
    def __init__(self,step_tot,size,System):
        self.step_tot = step_tot
        self.size = size
        self.System = System
        self.dr = (2*np.sqrt(self.System.ell_tot)) / self.size
        self.X = np.array([i*self.dr for i in range(self.size)])
    def average_r(self,r_array):
        """
        This function compute the distance between the [0.,0.,0.] linker 
        """
        av_dist = 0.
        # remove the [0.,0.,0.] linker :
        r_array = r_array[np.argwhere([np.all((r_array-np.array([0.,0.,0.]))!=0, axis=1)])[:,1]]
        N_linker = r_array.shape[0]
        # iterate over every linker
        for r in r_array:
            # compute the distance with [0.,0.,0.]
            av_dist+=np.linalg.norm(r)/N_linker
        return av_dist
    def I(self,r):
            #return the PR index corresponding to a position
            # the maximum posision 2*sqrt(System.ell)
            return int(r/(2*np.sqrt(self.System.ell_tot)) * self.size)
    def get_X(self,size_max):
        return np.array([i*self.dr for i in range(size_max)])
    def simulate_eq_distribution(self):
        self.PR = np.zeros(self.size,dtype=float)
        self.dt = np.zeros(self.step_tot,dtype=float)
        self.move = np.zeros(4,dtype=float)
        prev_R = self.average_r(self.System.get_r())
        for i in range(self.step_tot):
            movetype,Dt = self.System.evolve()
            self.dt[i] = Dt
            if np.isnan(Dt):
                raise ValueError
            try:
                self.PR[self.I(prev_R)]+=Dt/max(prev_R**2,1)
            except IndexError:
                self.PR.resize(self.I(prev_R)+1)
                self.PR[self.I(prev_R)]+=Dt/max(prev_R**2,1)
            prev_R = self.average_r(self.System.get_r())
            self.move[movetype] +=1            
        self.time = np.cumsum(self.dt)
        self.PR = self.PR/self.time[-1]
        self.move = self.move/self.step_tot           
        self.Mean_distance = np.sum(self.PR*self.get_X(self.PR.shape[0]))
    def Pmeet(r,l,):
        raise NotImplementedError
    def compute_average_Pmeet(self):
        """
        This function compute the average meeting probability within the whole cubic
        volume that define the region where linkers are generated
        """
        raise NotImplementedError
