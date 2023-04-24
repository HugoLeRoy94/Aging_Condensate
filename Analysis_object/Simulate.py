import numpy as np
import sys
sys.path.append('/home/hcleroy/PostDoc/aging_condensates/Simulation/Gillespie/Gillespie_backend')
import Gillespie_backend as Gil

class ParallelSimulation:
    """
    This class is used to quickly simulate a system, extract some parameters on several nodes.
    
    output :    is a dictionnary, to obtain the list of possible outputs call possible_outputs. the last
                argument is a list that tells if each of the previous output should be averaged or not.
                
    """
    def __init__(self,
                 gillespie_param,
                 simulation_param,
                 output,
                 cpu_param):
        self.gillespie_param = gillespie_param
        self.step_tot = simulation_param['step_tot']
       
        self.Nnodes = cpu_param['Nnodes']

        self.output = output
    def possible_outputs(self):
        print('R, r, ell average')
        return {'R':False,'r':False,'ell':False,'average':[False,False,False]}
    def actualize_outputs(self,output):
        self.output = output
    
    def extract_parameter(self,gillespie):
        """
        given a dictionnary of outputs, and a gillespie object, this function
        extract the relevant parameters
        """
        for key,value in self.output:
            if key == 'R':
                
