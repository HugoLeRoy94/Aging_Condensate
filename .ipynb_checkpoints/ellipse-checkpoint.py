import numpy as np

def ellipse(ax,ctr):
    """
    ax a 3x3 array/list with the coordinate of each demi ax
    ctr is a 3D array/list with the coordinate of the center of mass
    """
    # points on unit sphere
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    z = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    x = np.outer(np.ones_like(u), np.cos(v))
    
    # transform points to ellipsoid
    for i in range(len(x)):
        for j in range(len(x)):
               x[i,j], y[i,j], z[i,j] = ctr + np.dot(ax,
                                                      [x[i,j],y[i,j],z[i,j]])
    
    return x,y,z
def construct_axes_from_main_axe(u,b):
    """
    Given an main ax for the ellipse : u
    we construct two other normalize orthogonal axes.
    The function returns all the three axes.
    b is the demi-small axe of the ellipse
    """
    if u[2]!=0:
        v = np.array([1,1,(-u[0]-u[1])/u[2]])
        w = np.array([u[1]*(-u[0]-u[1])/u[2]-u[2], 
                  u[2]-u[0]*(-u[0]-u[1])/u[2],
                 u[0]-u[1]])
    else:
        v = np.array([0,0,1])
        w = np.array([u[1],u[0],0])
    
    v = v/np.linalg.norm(v)*b
    w = w/np.linalg.norm(w)*b
    return [u,v,w]