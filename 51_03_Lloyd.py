# python3
import sys
from math import *
import numpy as np
from copy import deepcopy

'''
Implement the Lloyd algorithm for k-means clustering.
     Input: Integers k and m followed by a set of points Data in m-dimensional space.
     Output: A set Centers consisting of k points (centers) resulting from applying the Lloyd algorithm to Data and Centers,
     where the first k points from Data are selected as the first k centers.
Sample Input:
2 2
1.3 1.1
1.3 0.2
0.6 2.8
3.0 3.2
1.2 0.7
1.4 1.6
1.2 1.0
1.2 1.1
0.6 1.5
1.8 2.6
1.2 1.3
1.2 1.0
0.0 1.9
Sample Output:
1.800 2.867
1.060 1.140

The Lloyd algorithm is one of the most popular clustering heuristics for the k-Means Clustering Problem. It first chooses k 
arbitrary points Centers from Data as centers and then iteratively performs the following two steps (see figure below):

Centers to Clusters: After centers have been selected, assign each data point to the cluster corresponding to its nearest 
center; ties are broken arbitrarily. 
Clusters to Centers: After data points have been assigned to clusters, assign each cluster’s center of gravity to be the 
cluster’s new center.
'''

class Lloyd:
    def __init__(self):
        data, k = self.readFromFile()
        centers = self.findCenters(data, k)
        self.saveResult(centers)

    def readFromFile(self):
        f = open('input.txt', 'r')
        raw = f.read().strip().split()
        f.close()
        k, m = int(raw[0]), int(raw[1])
        raw = raw[2:]
        data = np.zeros((len(raw)//m, m))
        for i in range(len(raw)//m):
            data[i, ] = [float(d) for d in raw[i*m:(i+1)*m]]
        return data, k

    def getDist(self, p1, p2):
        # calculate the distance between two points
        return sqrt(sum((p1-p2)**2))

    def findCenters(self, data, k, TOL = 1e-6, MAXITER = 100):
        n, m = data.shape
        centers = data[:k, ]
        for _ in range(MAXITER):
            prevCenters = deepcopy(centers)
            # centers to clusters
            clusters = [np.argmin([self.getDist(data[i, ], prevCenters[ii, ]) for ii in range(k)]) for i in range(n)]
            # clusters to centers
            centers = np.zeros((k, m))
            counts = np.zeros(k, dtype = int)
            for i in range(n):
                centers[clusters[i], ] += data[i, ]
                counts[clusters[i]] += 1
            for i in range(k):
                centers[i, ] /= counts[i]
            # has converged?
            if np.linalg.norm(centers - prevCenters) < TOL:
                break
        return centers
    
    def saveResult(self, centers):
        k = centers.shape[0]
        f = open('result.txt', 'w')
        for i in range(k):
            print(' '.join([str(p) for p in centers[i, ]]))
            f.write(' '.join([str(p) for p in centers[i, ]]) + '\n')
        f.close()

if __name__ == '__main__':
    Lloyd()