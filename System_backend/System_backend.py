from ctypes.wintypes import POINT
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
lib = cdll.LoadLibrary(str(pathlib.Path(__file__).parent.absolute())+'/../System_cpp/lib.so')

lib.create_system.argtypes=[c_double,c_double,c_double,c_double,c_int,c_bool,c_int,c_int]
lib.create_system.restype=POINTER(c_void_p)
lib.CopySystem.argtypes = [POINTER(c_void_p)]
lib.CopySystem.restype = POINTER(c_void_p)

lib.evolve.argtypes=[POINTER(c_void_p),POINTER(c_int)]
lib.evolve.restype=c_double
lib.delete_System.argtypes = [POINTER(c_void_p)]
lib.reset_crosslinkers.argypes = POINTER(c_void_p)


lib.get_F.argtypes = [POINTER(c_void_p)]
lib.get_F.restype = c_double
lib.get_r_size.argtypes=[POINTER(c_void_p)]
lib.get_r_size.restype=c_int
lib.get_r_system_size.argtypes=[POINTER(c_void_p)]
lib.get_r_system_size.restype=c_int
lib.get_r.argtypes = [POINTER(c_void_p),POINTER(c_double),c_int]
lib.get_r_system.argtypes = [POINTER(c_void_p),POINTER(c_double),c_int]
lib.get_N_strand.argtypes=[POINTER(c_void_p)]
lib.get_N_strand.restype=c_int
lib.get_R.argtypes=[POINTER(c_void_p),POINTER(c_double),c_int]
lib.get_ell_coordinates.argtypes=[POINTER(c_void_p),POINTER(c_double),c_int]
lib.get_ell.argtypes=[POINTER(c_void_p),POINTER(c_double),c_int]

lib.Print_Loop_positions.argtypes=[POINTER(c_void_p)]
lib.print_random_stuff.argtypes=[POINTER(c_void_p)]

class System:
    def __init__(self,ell_tot=100,rho0=0.1,BindingEnergy=-1.,kdiff=1.,seed=19874,sliding=False,Nlinker=0,old_system=None,dimension=3):
        self.move_types = {0 : 'unbind', 1:'diffuse', 2:'slide', 3:'bind'}
        if old_system is None:
            self.ell_tot,self.rho0,self.binding_energy,self.k_diff = ell_tot,rho0,BindingEnergy,kdiff
            self.Nlinker,self.dimension = Nlinker,dimension
            self.Address = lib.create_system(self.ell_tot,self.rho0,self.binding_energy,self.k_diff,seed,sliding,self.Nlinker,self.dimension)            
        else:
            self.copy(old_system)

    def copy(self,old_system):
        raise NotImplemented
        self.ell_tot = old_system.ell_tot
        self.rho0 = old_system.rho0
        self.binding_energy = old_system.binding_energy
        self.k_diff = old_system.k_diff
        self.Address = lib.CopySystem(old_system.Address)
    
    def __del__(self):
        lib.delete_System(self.Address) # deleting pointers is important in c++
    
    def evolve(self,steps=1):
        """
        Make the system evolves 
        """
        if steps>1:
            binds,time = np.zeros(steps,dtype=int),np.zeros(steps,dtype=float)
            for dt in range(steps):
                bind = c_int(0)
                time[dt] = lib.evolve(self.Address,byref(bind))
                binds[dt] = bind.value
            return binds,time
        else:
            bind = c_int(0)
            time = lib.evolve(self.Address,byref(bind))
            return bind.value, time
    def get_N_loop(self):
        return lib.get_N_strand(self.Address)
    
    def get_F(self):
        return lib.get_F(self.Address)
    
    def get_R(self):
        size = (self.get_N_loop())*3
        R = np.zeros(size,dtype=np.double)
        lib.get_R(self.Address,R.ctypes.data_as(POINTER(c_double)),size)
        R = np.reshape(R, (-1, 3))
        return R
    
    def get_ell_coordinates(self):
        size = self.get_N_loop()
        ell_coordinates = np.zeros(size,dtype=np.double)
        lib.get_ell_coordinates(self.Address,ell_coordinates.ctypes.data_as(POINTER(c_double)),size)
        return ell_coordinates
    
    def get_ell(self):
        size = self.get_N_loop()
        ell = np.zeros(size,dtype=np.double)
        lib.get_ell(self.Address,ell.ctypes.data_as(POINTER(c_double)),size)
        return ell

    def get_r(self):    
        size = lib.get_r_system_size(self.Address)
        r = np.zeros(size,dtype=np.double)
        # r_system access the linkers that are owned by the system
        lib.get_r_system(self.Address,r.ctypes.data_as(POINTER(c_double)),size)
        return np.reshape(r,(-1,3))
    
    def get_r_from_loops(self):
        size = lib.get_r_size(self.Address)
        r = np.zeros(size,dtype=np.double)
        lib.get_r(self.Address,r.ctypes.data_as(POINTER(c_double)),size)
        return np.reshape(r,(-1,3))
    
    def get_r_size(self):
        return lib.get_r_size(self.Address)
    
    def get_r_system_size(self):
        return lib.get_r_system_size(self.Address)
    
    def Print_loop_positions(self):
        lib.Print_Loop_positions(self.Address)
    
    def get_links_linear_position(self, bins=None):
        """
        This function returns the position of the links along the polymer
        """
        position = np.zeros(self.get_N_loop(),dtype=float)
        for n,l in enumerate(self.get_ell()):
            position[n] = position[n-1]+l
        return position
    
    def Plot3DSystem(self,*arg,**kwargs):
        fig = plt.figure(*arg,**kwargs)
        ax = fig.add_subplot(projection='3d')
        Ell = self.get_ell()
        R = self.get_R()
        r = self.get_r()
        # plot all the loops
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
        # plot the dangling end :
        Ell_dangling = self.ell_tot -sum( Ell )
        xel,yel,zel = ell.ellipse_from_main_ax([0,0,0],np.sqrt(Ell_dangling),np.sqrt(Ell_dangling),R[-1])
        ax.plot_wireframe(xel, yel, zel,  rstride=4, cstride=4, color='r', alpha=0.2)
        ax.set_xlim(min(np.append(R,r)),max(np.append(r,R)))
        ax.set_ylim(min(np.append(R,r)),max(np.append(r,R)))
        ax.set_zlim(min(np.append(R,r)),max(np.append(r,R)))
        ax.scatter(r[:,0],r[:,1],r[:,2],color='black')
        ax.scatter(R[:,0],R[:,1],R[:,2],s=50,color='red')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        return fig,ax

    def print_random_stuff(self):
        lib.print_random_stuff(self.Address)

    def reset_crosslinkers(self):
        lib.reset_crosslinkers(self.Address)
    
    def compute_correlation_function(self,Array,distmax=None,bins=10):
        """
        This function return the radial distribution function in any dimension

        Parameter: Array : 
                SystemSize:
                bins:
                distmax:
        return:
                P:
                R:
        """
        #--------------------------------------------------
        # construct a list of all the neighboring distances
        #--------------------------------------------------
        distances = list()
        for r in Array:
                for R in Array:
                    if distmax:
                        if np.linalg.norm(r-R)<distmax:
                            distances.append(np.linalg.norm(r-R))
                    else:        
                        distances.append(np.linalg.norm(r-R))
        distances = np.array(distances)
        #--------------------------------------------------
        # Compute the histogram of distances
        #--------------------------------------------------
        H,R  = np.histogram(distances,bins=bins,density=True)
        R = np.array([(R[i]+R[i+1])/2 for i in range(R.__len__()-1)]) # re-center the values of the bins
        #--------------------------------------------------
        try:
            iter(Array[0])
        except TypeError:
            Norm=np.array([self.rho0 for _ in range(R.__len__())])
        else:
            if Array[0].__len__()==2:
                Norm = np.array([2*np.arccos(-1)*r*self.rho0 for r in R])
            if Array[0].__len__()==3:
                Norm = np.array([4*np.arccos(-1)*r**2*self.rho0 for r in R])
        return H/Norm,R