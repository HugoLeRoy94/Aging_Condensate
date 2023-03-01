import numpy as np
from numpy.linalg import norm

pi = np.arcos(-1)

class Gillespie:
    def __init__(self,xy,L,*args,**kwargs):
        """
        The gillespie is a polymer chain that is bound at its two ends.
        The state of the gillespie is defined by a list of 3d vectors that defined
        the coordinates of binding points, and a distance L that correspond to
        the length of polymer associated to each distance.
        We start with a single R and L.
        xy : coordinates of the second anchoring point. the first is (0,0)
        L   : length of polymer in between these two points
        """
        self.chains = list(Chain(L,xy))
        self._beta = 0.
        try:
            self.DE = kwargs.get('DE')
        except:
            pass
    def _ComputeProcessesList(self):
        """
        Private function to compute the list of available process and their
        associated rates.
        """
        self._process_list = list()
        if self.chains.__len__()!=1:
            for i in range(self.chains.__len__()-1):
                # compute the breaking rates
                DS = 1.5*np.log(0.75*pi*self.chains[i].L)-
                    ((self.chains[i].R[0]-self.chains[i+1].R[0])**2+
                    (self.chains[i].R[0]-self.chains[i+1].R[0])**2)/
                    2*self.chains[i].L
                self._process_list.append()
                # compute the binding rates
    def _DrawProcess(self):
        """
        Private function to draw a random process within the list of available
        ones. change the gillespie accordingly, and return the time increment.
        """
    def Evolv(self):
        self._ComputeProcessesList()
        time_increment = _DrawProcess()
        return time_increment
    def UpdateTemperature(self,Beta):
        self._beta = Beta

class Chain:
    def __init__(self,L,R):
    """
    Chaine is a polymeric chain that is anchored at two points : in (0,0), and
    in R.
    """
    if norm(R) > L:
        print('length of the polymer smaller than the two anchoring points')
        raise ValueError
    self.L = L
    self.R = R
    # compute the volume occupied :
    self.V = 2 * pi*np.sqrt(L**2 - norm(R)**2)*L
