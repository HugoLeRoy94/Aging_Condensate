import numpy as np
import sys
#sys.path.append('/home/hcleroy/PostDoc/aging_condensates/Simulation/Gillespie/System_backend')
#import System_backend as Sys
def Vol(R) : 
    return np.arccos(-1)*4./3. * R**3
def rescale_with_r(H,R):
    """
    R.shape[0] = H.sape[0] +1 
    R is the list of the bar's extremities 
    """
    Hres = np.zeros(H.shape,dtype=float)
    Rres = np.zeros(R.shape,dtype=float)
    for i in range(H.shape[0]):
        Hres[i] = H[i] / (Vol(R[i+1])-Vol(R[i]))
        Rres[i] = 0.5*(R[i+1]+R[i])
    return Hres,Rres
class Simulation:
    def __init__(self,step_tot,size,L_size,System):
        """
        step_tot is the total number of step in the system
        size is the number of discretization of space
        """
        self.step_tot = step_tot
        self.size = size
        self.L_size = L_size
        self.System = System
        #self.dr = (2*np.sqrt(self.System.ell_tot)) / self.size
        self.dr = np.sqrt(System.ell_tot/2)/self.size
        self.X = np.array([i*self.dr for i in range(self.size)])
    def average_distance(self,r_array):
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
    def average_ell_coordinate(self):
        return np.mean(self.System.get_ell_coordinates())
    def I(self,r):
            #return the PR index corresponding to a position
            # the maximum posision 2*sqrt(System.ell)
            return int(r/self.dr)
    def get_X(self,size_max):
        return np.array([i*self.dr for i in range(size_max)])
    def simulate_eq_distribution(self):
        """
        This function computes the probability density to be at a certain distance
        from the centerlinker. 
        """
        self.PL = np.zeros(self.L_size,dtype=float)
        self.PR = np.zeros(self.size,dtype=float)
        self.dt = np.zeros(self.step_tot,dtype=float)
        self.move = np.zeros(4,dtype=float)
        prev_R = self.average_distance(self.System.get_r())
        tot_bound_time=0
        for i in range(self.step_tot):
            movetype,Dt = self.System.evolve()
            self.dt[i] = Dt
            if np.isnan(Dt):
                raise ValueError
            try:
                if self.I(prev_R)!=0:
                    self.PR[self.I(prev_R)]+=Dt#/(4*np.arccos(-1)*self.dr*(self.I(prev_R)*self.dr)**2)
                else :
                    self.PR[self.I(prev_R)]+=Dt#/(4/3*np.arccos(-1)*self.dr**3)            
            except IndexError:
                self.PR.resize(self.I(prev_R)+1)
                self.PR[self.I(prev_R)]+=Dt/(4*np.arccos(-1)*self.dr*(self.I(prev_R)*self.dr)**2)
            if self.System.move_types[movetype] == 'unbind':
                try:
                    self.PL[int(self.average_ell_coordinate()/(self.System.ell_tot/self.L_size))-1]+=Dt
                except IndexError:
                    print(self.System.get_ell())
                    print(self.System.get_ell()/(self.System.ell_tot/self.L_size))
                    raise
                tot_bound_time += Dt
            prev_R = self.average_distance(self.System.get_r())
            self.move[movetype] +=1
        self.PL = self.PL/tot_bound_time
        self.time = np.cumsum(self.dt)
        self.PR = self.PR/self.time[-1]
        self.PR,self.IX = rescale_with_r(self.PR,np.array([i*self.dr for i in range(self.PR.shape[0]+1)]))
        self.move = self.move/self.step_tot           
        self.Mean_distance = self.integrate_density(self.PR*self.get_X(self.PR.shape[0]))#np.sum(self.PR*self.get_X(self.PR.shape[0]))
    def integrate_density(self,P):
        return np.sum([P[i]*(Vol(self.dr*(i+1))-Vol(self.dr*i)) for i in range(P.shape[0])])
    def density_to_probability(self):
        return np.array([self.PR[i]*(4*np.arccos(-1)*self.dr*(i*self.dr)**2) for i in range(self.PR.shape[0])])
    def Pmeet(r,l,):
        raise NotImplementedError
    def compute_average_Pmeet(self):
        """
        This function compute the average meeting probability within the whole cubic
        volume that define the region where linkers are generated
        """
        raise NotImplementedError
