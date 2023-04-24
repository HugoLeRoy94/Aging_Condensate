import numpy as np

import math
from scipy.sparse import csr_matrix


import sknetwork as skn


from matplotlib.colors import LinearSegmentedColormap
from collections import defaultdict
import random
cdict = {'blue':   ((0.0,  0.9, 0.9),
                    (0.5,  0.4, 0.4),
                    (1.0,  0.1, 0.1)),

         'green': ((0.0,  0.5, 0.5),
                   (0.5, 1, 1),
                   (1.0,  0.3, 0.3)),

         'alpha': ((0.0,  1, 1),
                   (0.5, 0.8, 0.8),
                   (1.0,  1, 1)),

         'red':  ((0.0,  0.4, 0.4),
                  (0.5,  0.5, 0.5),
                  (1.0,  0.9, 0.9)),
         }
def weight(distance):
    return 1/distance
def make_adjacency_matrix(array):
    adjacency = np.zeros((array.shape[0], array.shape[0]),dtype=float)
    for i in range(array.shape[0]):
        for j in range(i+1,array.shape[0]):
            adjacency[i,j] = weight(np.linalg.norm(array[i]-array[j]))
            adjacency[j,i] = adjacency[i,j]
    return adjacency
class Cluster:
    def __init__(self, array):
        """
        This function takes  an array of NxD number, where N is
        the number particles, and D the dimension. It first create
        and adjacency matrix, where the weight is proportional to
        the inverse of the distance between each points.
        """

        self.Nparticle = array.shape[0]
        self.dimension = array.shape[1]

        self.adjacency = make_adjacency_matrix(array)
    def make_cluster(self,resolution=1.):
        louvain = skn.clustering.Louvain(resolution=resolution)
        labels = louvain.fit_transform(self.adjacency)
        #propagation = skn.clustering.PropagationClustering()
        #labels = propagation.fit_transform(AdjacencyMatrix)
        return labels
