import numpy as np

def ellipse(ax,ctr,OmZ,OmY):
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
               x[i,j], y[i,j], z[i,j] = ctr + np.dot(OmY,np.dot(OmZ,np.dot(ax,
                                                      [x[i,j],y[i,j],z[i,j]])))
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
    if np.dot(u,v)>10**-10 :
        print("u,v not perpendicular")
        print("u.v = "+str(np.dot(u,v)))
    if np.dot(u,w)>10**-10 :
        print("u,v not perpendicular")
        print("u.v = "+str(np.dot(u,w)))
    if np.dot(w,v)>10**-10 :
        print("u,v not perpendicular")
        print("u.v = "+str(np.dot(w,v)))    
    return [u,v,w]

def omegaY(theta):
    res = np.zeros((3,3),dtype=float)
    res[0,0] = np.cos(theta)
    res[0,2] = np.sin(theta)
    res[1,1] = 1
    res[2,0] = -np.sin(theta)
    res[2,2] = np.cos(theta)
    return res
def omegaZ(theta):
    res = np.zeros((3,3),dtype=float)
    res[0,0] = np.cos(theta)
    res[0,1] = -np.sin(theta)
    res[2,2] = 1
    res[1,0] = np.sin(theta)
    res[1,1] = np.cos(theta)
    return res

def ellipse_from_main_ax(main_ax,b,ctr):
    """
    main_ax is the vector of the main ax of the ellipse.
    it gives botjh the orientation and the norm.
    """
    theta = np.arctan2(main_ax[1],main_ax[0])
    phi = np.arctan2(main_ax[0],main_ax[2])-np.arccos(-1)/2
    
    OmZ = omegaZ(theta)
    OmY = omegaY(phi)
    
    axes = construct_axes_from_main_axe(main_ax,b)
    
    return ellipse(axes,ctr,OmZ,OmY)