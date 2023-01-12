import numpy as np

def compute_adjacency_matrix(linear_position):
    """
    This function transform a set of linear coordinate into an adjacency coordinates
    """
    A = np.zeros((linear_position.shape[0],linear_position.shape[0]),dtype=float)
    for i,pos1 in enumerate(linear_position):
        for j,pos2 in enumerate(linear_position):
            if np.linalg.norm(pos1-pos2)!=0:
                A[i,j] = 1/np.linalg.norm(pos1-pos2)
            else :
                A[i,j] = 0.
    return A
def compute_clustering(A):
    return np.matrix.trace(np.dot(A,np.dot(A,A)))/np.sum(np.sum(A,axis=1)*(np.sum(A,axis=1)-1))