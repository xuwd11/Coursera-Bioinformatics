# python3

import sys
import numpy as np
from copy import deepcopy

'''
If we had a molecular clock measuring evolutionary time, then we could assign an age to every node v in a rooted binary tree 
(denoted Age(v)), where all of the leaves of the tree have age 0 because they correspond to present-day species. We could then 
define the weight of an edge (v, w) in the tree as the difference Age(v) - Age(w). Consequently, the length of a path between 
the root and any node would be equal to the difference between their ages. Such a tree, in which the distance from the root to 
any leaf is the same, is called ultrametric.

Our aim is to derive an ultrametric tree that explains a given distance matrix (even if it does so only approximately). UPGMA 
(which stands for Unweighted Pair Group Method with Arithmetic Mean) is a simple clustering heuristic that introduces a hypothetical 
molecular clock for constructing an ultrametric evolutionary tree.

Given an n × n matrix D, UPGMA (which is illustrated in the figure on the next step) first forms n trivial clusters, each containing 
a single leaf. The algorithm then finds a pair of “closest” clusters. To clarify the notion of closest clusters, UPGMA defines the 
distance between clusters C1 and C2 as the average pairwise distance between elements of C1 and C2.

Once UPGMA has identified a pair of closest clusters C1 and C2, it merges them into a cluster C with |C1| + |C2| elements and then 
creates a node for C, which it connects to each of C1 and C2 by a directed edge. The age of C is set to be DC1, C2 /2. UPGMA then 
iterates this process of merging the two closest clusters until only a single cluster remains, which corresponds to the root.

Implement UPGMA.
     Input: An integer n followed by a space separated n x n distance matrix.
     Output: An adjacency list for the ultrametric tree returned by UPGMA. Edge weights should be accurate to two decimal places
     (answers in the sample dataset below are provided to three decimal places).

Note on formatting: The adjacency list must have consecutive integer node labels starting from 0. The n leaves must be labeled 
0, 1, ..., n - 1 in order of their appearance in the distance matrix. Labels for internal nodes may be labeled in any order but 
must start from n and increase consecutively.

Sample Input:
4
0	20	17	11
20	0	20	13
17	20	0	10
11	13	10	0
Sample Output:
0->5:7.000
1->6:8.833
2->4:5.000
3->4:5.000
4->2:5.000
4->3:5.000
4->5:2.000
5->0:7.000
5->4:2.000
5->6:1.833
6->5:1.833
6->1:8.833

UPGMA(D, n) 
        Clusters ← n single-element clusters labeled 1, ... , n
        construct a graph T with n isolated nodes labeled by single elements 1, ... , n
    for every node v in T 
        Age(v) ← 0
    while there is more than one cluster
        find the two closest clusters Ci and Cj
        merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
        add a new node labeled by cluster Cnew to T
        connect node Cnew to Ci and Cj by directed edges 
        Age(Cnew) ← DCi, Cj / 2
        remove the rows and columns of D corresponding to Ci and Cj
        remove Ci and Cj from Clusters
        add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters
        add Cnew to Clusters
    root ← the node in T corresponding to the remaining cluster
    for each edge (v, w) in T
        length of (v, w) ← Age(v) - Age(w)
    return T
'''

class UPGMA:
    def __init__(self):
        #n, disMatrix = self._input()
        n, disMatrix = self.readFromFile()
        adj = self.runUPGMA(disMatrix, n)
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

    def runUPGMA(self, disMatrix, n):
        D = np.array(disMatrix, dtype = float)
        np.fill_diagonal(D, np.inf)        
        clusters = [[i, 1] for i in range(n)]
        adj = [[] for i in range(n)]
        age = [0. for i in range(n)]
        if len(D) <= 1:
            return adj
        while True:
            index = np.argmin(D)
            i = index // len(D)
            j = index % len(D)
            i_new = len(adj)
            adj.append([])
            C_new = [i_new, clusters[i][1] + clusters[j][1]]
            adj[i_new].append(clusters[i][0])
            adj[i_new].append(clusters[j][0])
            adj[clusters[i][0]].append(i_new)
            adj[clusters[j][0]].append(i_new)
            age.append(D[i, j] / 2)

            if 2 == len(D):
                break

            d_new = (D[i,:]*clusters[i][1] + D[j,:]*clusters[j][1]) / (clusters[i][1]+clusters[j][1])
            d_new = np.delete(d_new, [i, j], 0)
            D = np.delete(D, [i, j], 0)
            D = np.delete(D, [i, j], 1)
            D = np.insert(D, len(D), d_new, axis = 0)
            d_new = np.insert(d_new, len(d_new), np.inf, axis = 0)
            D = np.insert(D, len(D)-1, d_new, axis = 1)

            if i < j:
                del clusters[j]
                del clusters[i]
            else:
                del clusters[i]
                del clusters[j]
            
            clusters.append(C_new)

        adjL = deepcopy(adj)
        for i, nodes in enumerate(adj):
            for j, v in enumerate(nodes):
                adjL[i][j] = (v, abs(age[i]-age[v]))
        return adjL

if __name__ == "__main__":
    UPGMA()