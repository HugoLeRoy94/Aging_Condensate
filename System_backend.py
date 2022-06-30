import numpy as np
import ellipse as ell # module to plot the system
import matplotlib.pyplot as plt
import pathlib
from ctypes import cdll
from ctypes import c_double
from ctypes import c_int
from ctypes import POINTER
from ctypes import c_void_p
from ctypes import c_char_p
from ctypes import c_bool
from ctypes import byref
lib = cdll.LoadLibrary(str(pathlib.Path(__file__).parent.absolute())+'/lib.so')

lib.create_system.argtypes=[c_double,c_double,c_double,c_double,c_int]
lib.create_system.restype=POINTER(c_void_p)

lib.evolve.argtypes=[POINTER(c_void_p),POINTER(c_bool)]
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
    def __init__(self,ell_tot,distance,rho0,temperature,seed=19874,rho_adjust=True):
        self.ell_tot,self.D,self.rho0,self.T = ell_tot,distance,rho0,temperature
        self.Address = lib.create_system(self.ell_tot,self.D,self.rho0,self.T,seed,rho_adjust)
    def evolve(self):
        bind = c_bool(False)
        time = lib.evolve(self.Address,byref(bind))
        return bind.value, time
    def get_N_loop(self):
        return lib.get_N(self.Address)
    def get_R(self):
        size = (self.get_N_loop()+1)*3
        R = np.zeros(size,dtype=np.double)
        lib.get_R(self.Address,R.ctypes.data_as(POINTER(c_double)),size)
        R = np.reshape(R, (-1, 3))
        #return R[np.argsort(R[:,0])]
        return R
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
    def get_links_linear_position(self, bins=None):
        """
        This function returns the position of the links along the polymer
        """
        position = np.zeros(self.ell_tot,dtype=float)
        for n,l in enumerate(self.get_ell()):
            position[n] = position[n-1]+l
        return position
    def Plot3DSystem(self,*arg,**kwargs):
        fig = plt.figure(*arg,**kwargs)
        ax = fig.add_subplot(projection='3d')
        Ell = self.get_ell()
        R = self.get_R()
        r = self.get_r()
        for i in range(0,R.shape[0]-1):
            #axes = ell.construct_axes_from_main_axe((R[i+1]-R[i])/2,np.sqrt(Ell[i])/2)
            #print(axes)
            if np.linalg.norm(R[i+1]-R[i]) < Ell[i]*0.1:
                a = np.sqrt(Ell[i])*0.5
                b = np.sqrt(Ell[i])*0.5

            else:
                a = np.linalg.norm(R[i+1]-R[i])/2
                b = np.sqrt(Ell[i])*0.5
            xel,yel,zel = ell.ellipse_from_main_ax(-(R[i+1]-R[i])/2,a,b,[0.5*(R[i,0]+R[i+1,0]),0.5*(R[i+1,1]+R[i,1]),(R[i+1,2]+R[i,2])*0.5])
            ax.plot_wireframe(xel, yel, zel,  rstride=4, cstride=4, color='#2980b9', alpha=0.2)
        ax.set_xlim(min(R[:,0]),max(R[:,0]))
        ax.set_ylim(-25,25)
        ax.set_zlim(-25,25)
        ax.scatter(r[:,0],r[:,1],r[:,2],color='black')
        ax.scatter(R[:,0],R[:,1],R[:,2],s=50,color='red')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        return fig,ax
