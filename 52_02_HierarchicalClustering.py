# python3

import sys
import numpy as np
from copy import deepcopy

'''
Implement HierarchicalClustering.
     Input: An integer n, followed by an n x n distance matrix.
     Output: The result of applying HierarchicalClustering to this distance matrix (using Davg), with each newly created cluster listed
     on each line.
Sample Input:
7
0.00 0.74 0.85 0.54 0.83 0.92 0.89
0.74 0.00 1.59 1.35 1.20 1.48 1.55
0.85 1.59 0.00 0.63 1.13 0.69 0.73
0.54 1.35 0.63 0.00 0.66 0.43 0.88
0.83 1.20 1.13 0.66 0.00 0.72 0.55
0.92 1.48 0.69 0.43 0.72 0.00 0.80
0.89 1.55 0.73 0.88 0.55 0.80 0.00
Sample Output:
4 6
5 7
3 4 6
1 2
5 7 3 4 6
1 2 5 7 3 4 6

HierarchicalClustering, whose pseudocode is shown below, progressively generates n different partitions of the underlying data into clusters, 
all represented by a tree in which each node is labeled by a cluster of genes. The first partition has n single-element clusters represented 
by the leaves of the tree, with each element forming its own cluster. The second partition merges the two “closest” clusters into a single 
cluster consisting of two elements. In general, the i-th partition merges the two closest clusters from the (i - 1)-th partition and has 
n - i + 1 clusters. We hope this algorithm looks familiar — it is UPGMA (from the chapter on evolutionary tree constuction) in disguise.

HierarchicalClustering(D, n)
    Clusters ← n single-element clusters labeled 1, ... , n 
       construct a graph T with n isolated nodes labeled by single elements 1, ... , n 
    while there is more than one cluster 
        find the two closest clusters Ci and Cj
        merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
        add a new node labeled by cluster Cnew to T
        connect node Cnew to Ci and Cj by directed edges
        remove the rows and columns of D corresponding to Ci and Cj
        remove Ci and Cj from Clusters
        add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters 
        add Cnew to Clusters 
    assign root in T as a node with no incoming edges
    return T

The distance function that we encountered with UPGMA uses the average distance between elements in two clusters,

Davg(C1,C2)=∑all points i in cluster C1 ∑all points j in cluster C2Di,j / (|C1|⋅|C2|)
'''

class HierarchicalClustering:
    def __init__(self):
        n, dist = self.readFromFile()
        adj, newClusters = self.clustering(n, dist)
        print('\n'.join([' '.join([str(c) for c in clusters]) for clusters in newClusters]))
        f = open('result.txt', 'w')
        f.write('\n'.join([' '.join([str(c) for c in clusters]) for clusters in newClusters]))
        f.close()
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().strip().split('\n')
        f.close()
        n = int(data[0])
        dist = np.array([[float(v) for v in d.split()] for d in data[1:]])
        np.fill_diagonal(dist, np.inf)
        return n, dist

    def clustering(self, n, dist):
        clusters = [[i, 1] for i in range(n)]
        newClusters = []
        adj = [[] for _ in range(n)]
        while len(dist) > 1:
            node_new = len(adj)
            index = np.argmin(dist)
            i = index // len(dist)
            j = index % len(dist)
            d_new = (dist[i, :] * clusters[i][1] + dist[j, :] * clusters[j][1]) / (clusters[i][1] + clusters[j][1])
            d_new = np.delete(d_new, [i, j], 0)
            dist = np.delete(dist, [i, j], 0)
            dist = np.delete(dist, [i, j], 1)
            dist = np.insert(dist, len(dist), d_new, 0)
            d_new = np.insert(d_new, len(d_new), np.inf, 0)
            dist = np.insert(dist, len(dist)-1, d_new, 1)
            adj.append([clusters[i][0], clusters[j][0]])
            clusters.append([node_new, clusters[i][1] + clusters[j][1]])
            if i < j:
                del clusters[j]
                del clusters[i]
            else:
                del clusters[i]
                del clusters[j]
            newClusters.append(self.findLeafs(adj, node_new))
        return adj, newClusters

    def findLeafs(self, adj, v):
        leafs = []
        visited = [False for _ in range(len(adj))]
        stack = []
        stack.append(v)
        while len(stack) > 0:
            v = stack.pop()
            if 0 == len(adj[v]):
                leafs.append(v + 1)
            if not visited[v]:
                visited[v] = True
                for w in adj[v]:
                    stack.append(w)
        return leafs

if __name__ == "__main__":
    HierarchicalClustering()