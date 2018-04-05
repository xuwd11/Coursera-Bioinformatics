# python3
import sys
from math import *
import numpy as np

'''
Implement the FarthestFirstTraversal clustering heuristic.
     Input: Integers k and m followed by a set of points Data in m-dimensional space.
     Output: A set Centers consisting of k points (centers) resulting from applying FarthestFirstTraversal(Data, k),
     where the first point from Data is chosen as the first center to initialize the algorithm.
Sample Input:
3 2
0.0 0.0
5.0 5.0
0.0 5.0
1.0 1.0
2.0 2.0
3.0 3.0
1.0 2.0
Sample Output:
0.0 0.0
5.0 5.0
0.0 5.0

Although the k-Center Clustering Problem is easy to state, it is NP-Hard. The Farthest First Traversal heuristic, whose 
pseudocode is shown below, selects centers from the points in Data (instead of from all possible points in m-dimensional 
space). It begins by selecting an arbitrary point in Data as the first center and iteratively adds a new center as the point 
in Data that is farthest from the centers chosen so far, with ties broken arbitrarily.

FarthestFirstTraversal(Data, k) 
    Centers ← the set consisting of a single randomly chosen point from Data
    while |Centers| < k 
        DataPoint ← the point in Data maximizing d(DataPoint, Centers) 
        add DataPoint to Centers 
    return Centers 
'''

class FarthestFirstTraversal:
    def __init__(self):
        data, k = self.readFromFile()
        centers = self.findCenters(data, k)
        self.saveResult(centers)
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        raw = f.read().strip().split()
        k, m = int(raw[0]), int(raw[1])
        raw = raw[2:]
        data = np.zeros((len(raw)//m, m))
        for i in range(len(raw)//m):
            data[i, ] = [float(d) for d in raw[i*m:(i+1)*m]]
        return data, k

    def getDist(self, p1, p2):
        # calculate the distance between two points
        return sqrt(sum((p1-p2)**2))
    
    def updateClusters(self, newCenter, data, assignedCenters, dists):
        for i in range(data.shape[0]):
            newDist = self.getDist(newCenter, data[i, :])
            if newDist < dists[i]:
                assignedCenters[i, ] = newCenter
                dists[i] = newDist
    
    def findCenters(self, data, k):
        n, m = data.shape
        centers = [data[0, ]]
        assignedCenters = np.tile(centers[0], (n, 1))
        dists = np.zeros(n)
        for i in range(n):
            dists[i] = self.getDist(data[i, ], centers[0])
        for _ in range(k - 1):
            i = np.argmax(dists)
            newCenter = data[i, :]
            self.updateClusters(newCenter, data, assignedCenters, dists)
            centers.append(newCenter)
        return centers

    def saveResult(self, centers):
        f = open('result.txt', 'w')
        for center in centers:
            print(' '.join([str(p) for p in center]))
            f.write(' '.join([str(p) for p in center]) + '\n')
        f.close()

if __name__ == '__main__':
    FarthestFirstTraversal()