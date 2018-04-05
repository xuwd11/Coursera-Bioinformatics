# python3
import sys
from math import *
import numpy as np
from copy import deepcopy

'''
Implement the expectation maximization algorithm for soft k-means clustering.
     Input: Integers k and m, followed by a stiffness parameter β, followed by a set of points Data in m-dimensional space.
     Output: A set Centers consisting of k points (centers) resulting from applying the expectation maximization algorithm for soft
     k-means clustering. Select the first k points from Data as the first centers for the algorithm and run the algorithm for 100
    E-steps and 100 M-steps. Results should be accurate up to three decimal places.
Sample Input:
2 2
2.7
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
1.662 2.623
1.075 1.148
'''

class SoftClustering:
    def __init__(self):
        data, k, stiffness = self.readFromFile()
        centers = self.findCenters(data, k, stiffness)
        self.saveResult(centers)

    def readFromFile(self):
        f = open('input.txt', 'r')
        raw = f.read().strip().split()
        f.close()
        k, m, stiffness = int(raw[0]), int(raw[1]), float(raw[2])
        raw = raw[3:]
        data = np.zeros((len(raw)//m, m))
        for i in range(len(raw)//m):
            data[i, ] = [float(d) for d in raw[i*m:(i+1)*m]]
        return data, k, stiffness

    def getDist(self, p1, p2):
        # calculate the distance between two points
        return sqrt(sum((p1-p2)**2))

    def findCenters(self, data, k, stiffness, MAXITER = 100):
        n, m = data.shape
        centers = data[:k, ]
        hiddenMatrix = np.zeros((k, n))
        for _ in range(MAXITER):
            # Centers to Soft Clusters (E-step)
            for i in range(k):
                for j in range(n):
                    hiddenMatrix[i, j] = np.exp(-stiffness*self.getDist(data[j, :], centers[i, :]))
            for j in range(n):
                hiddenMatrix[:, j] /= sum(hiddenMatrix[:, j])
            # Soft Clusters to Centers (M-step)
            centers = np.dot(hiddenMatrix, data)
            for i in range(k):
                centers[i, ] /= sum(hiddenMatrix[i, ])
        return centers

    def saveResult(self, centers):
        k = centers.shape[0]
        f = open('result.txt', 'w')
        for i in range(k):
            print(' '.join([str(p) for p in centers[i, ]]))
            f.write(' '.join([str(p) for p in centers[i, ]]) + '\n')
        f.close()      

if __name__ == '__main__':
    SoftClustering()