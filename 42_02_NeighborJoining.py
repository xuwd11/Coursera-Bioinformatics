# python3

import sys
import numpy as np

'''
In 1987, Naruya Saitou and Masatoshi Nei developed the neighbor-joining algorithm for evolutionary tree reconstruction. Given 
an additive distance matrix, this algorithm, which we call NeighborJoining, finds a pair of neighboring leaves and substitutes 
them by a single leaf, thus reducing the size of the tree. NeighborJoining can thus recursively construct a tree fitting the 
additive matrix. This algorithm also provides a heuristic for non-additive distance matrices that performs well in practice.

The central idea of NeighborJoining is that although finding a minimum element in a distance matrix D is not guaranteed to 
yield a pair of neighbors in Tree(D), we can convert D into a different matrix whose minimum element does yield a pair of 
neighbors. First, given an n × n distance matrix D, we define TotalDistanceD(i) as the sum ∑1≤k≤n Di,k of distances from 
leaf i to all other leaves. The neighbor-joining matrix D* (see below) is defined such that for any i and j, D*i,i = 0 
and D*i,j = (n - 2) · Di,j - TotalDistanceD(i) - TotalDistanceD(j).

Implement NeighborJoining.
     Input: An integer n, followed by an n x n distance matrix.
     Output: An adjacency list for the tree resulting from applying the neighbor-joining algorithm. Edge-weights should be 
     accurate to two decimal places (they are provided to three decimal places in the sample output below).

Note on formatting: The adjacency list must have consecutive integer node labels starting from 0. The n leaves must be 
labeled 0, 1, ..., n - 1 in order of their appearance in the distance matrix. Labels for internal nodes may be labeled 
in any order but must start from n and increase consecutively.

Sample Input:
4
0	23	27	20
23	0	30	28
27	30	0	30
20	28	30	0
Sample Output:
0->4:8.000
1->5:13.500
2->5:16.500
3->4:12.000
4->5:2.000
4->0:8.000
4->3:12.000
5->1:13.500
5->2:16.500
5->4:2.000
'''

class NeighborJoining:
    def __init__(self):
        n, disMatrix = self.readFromFile()
        adj = self.runNeighborJoining(disMatrix, n)
        self.printGraph(adj)
        self.saveResult(adj) 
        
    def _input(self):
        data = sys.stdin.read().strip().split('\n')
        n = int(data[0])
        distMatrix = [[0]*n for _ in range(n)]
        for i in range(n):
            d = data[i+1].split()
            for k in range(n):
                distMatrix[i][k] = int(d[k])
        return n, distMatrix
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        n = int(data[0])
        distMatrix = [[0]*n for _ in range(n)]
        for i in range(n):
            d = data[i+1].split()
            for k in range(n):
                distMatrix[i][k] = int(d[k])
        return n, distMatrix

    def saveResult(self, adj):
        f = open('result.txt', 'w')
        for i, nodes in enumerate(adj):
            for d, w in nodes:
                f.write(str(i)+'->'+str(d)+':'+'%0.3f' % w+'\n')

    def printDistMatrix(self, distMatrix):
        for d in distMatrix:
            print(' '.join([str(i) for i in d]))

    def printGraph(self, adj):
        for i, nodes in enumerate(adj):
            for d, w in nodes:
                print(str(i)+'->'+str(d)+':'+'%0.3f' % w)

    def runNeighborJoining(self, disMatrix, n):
        D = np.array(disMatrix, dtype = float)
        clusters = [i for i in range(n)]
        adj = [[] for i in range(n)]
        if len(D) <= 1:
            return adj
        while True:
            if 2 == n:
                adj[len(adj)-1].append((len(adj)-2, D[0][1]))
                adj[len(adj)-2].append((len(adj)-1, D[0][1]))
                break
            totalDist = np.sum(D, axis = 0)
            D1 = (n-2) * D
            D1 = D1 - totalDist
            D1 = D1 - totalDist.reshape((n, 1))
            np.fill_diagonal(D1, 0.)
            #print(D1)
            index = np.argmin(D1)
            i = index // n
            j = index % n
            delta = (totalDist[i] - totalDist[j])/(n-2)
            li = (D[i, j]+delta)/2
            lj = (D[i, j]-delta)/2
            d_new = (D[i, :]+D[j, :]-D[i, j])/2
            D = np.insert(D, n, d_new, axis = 0)
            d_new = np.insert(d_new, n, 0., axis = 0)
            D = np.insert(D, n, d_new, axis = 1)
            D = np.delete(D, [i, j], 0)
            D = np.delete(D, [i, j], 1)
            #print(D)

            m = len(adj)
            adj.append([])
            adj[m].append((clusters[i], li))
            adj[clusters[i]].append((m, li))
            adj[m].append((clusters[j], lj))
            adj[clusters[j]].append((m, lj))
            if i < j:
                del clusters[j]
                del clusters[i]
            else:
                del clusters[i]
                del clusters[j]
            clusters.append(m)
            
            n -= 1
        
        return adj

if __name__ == "__main__":
    NeighborJoining()