import numpy as np
import pathlib
from ctypes import cdll
from ctypes import c_double
from ctypes import c_int
from ctypes import POINTER
from ctypes import c_void_p
from ctypes import c_char_p
lib = cdll.LoadLibrary(str(pathlib.Path(__file__).parent.absolute())+'/lib.so')

lib.create_system.argtypes=[c_double,c_double,c_double,c_double]
lib.create_system.restype=POINTER(c_void_p)

lib.evolve.argtypes=[POINTER(c_void_p)]
lib.evolve.restype=c_double

lib.get_r_size.argtypes=[POINTER(c_void_p)]
lib.get_r_size.restype=c_int
lib.get_r.argtypes = [POINTER(c_void_p),POINTER(c_double),c_int]
lib.get_N.argtypes=[POINTER(c_void_p)]
lib.get_N.restype=c_int
lib.get_R.argtypes=[POINTER(c_void_p),POINTER(c_double),c_int]
lib.get_ell.argtypes=[POINTER(c_void_p),POINTER(c_double),c_int]
lib.Print_Loop_positions.argtypes=[POINTER(c_void_p)]

class System:
    def __init__(self,ell_tot,distance,rho0,temperature):
        self.ell_tot,self.D,self.rho0,self.T = ell_tot,distance,rho0,temperature
        self.Address = lib.create_system(self.ell_tot,self.D,self.rho0,self.T)
    def evolve(self):
        return lib.evolve(self.Address)
    def get_N_loop(self):
        return lib.get_N(self.Address)
    def get_R(self):
        size = (self.get_N_loop()+1)*3
        R = np.zeros(size,dtype=np.double)
        lib.get_R(self.Address,R.ctypes.data_as(POINTER(c_double)),size)
        return np.reshape(R, (-1, 3))
    def get_ell(self):
        size = self.get_N_loop()
        ell = np.zeros(size,dtype=np.double)
        lib.get_ell(self.Address,ell.ctypes.data_as(POINTER(c_double)),size)
        return ell
    def get_r(self):
        size = lib.get_r_size(self.Address)
        r = np.zeros(size,dtype=np.double)
        lib.get_r(self.Address,r.ctypes.data_as(POINTER(c_double)),size)
        return np.reshape(r,(-1,3))
    def get_r_size(self):
        return lib.get_r_size(self.Address)
    def Print_loop_positions(self):
        lib.Print_Loop_positions(self.Address)
