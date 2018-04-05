# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Solve the Nearest Neighbors of a Tree Problem.
     Input: Two internal nodes a and b specifying an edge e, followed by an adjacency list of an unrooted binary tree.
     Output: Two adjacency lists representing the nearest neighbors of the tree with respect to e. Separate the
     adjacency lists with a blank line.

Sample Input:
5 4
0->4
4->0
1->4
4->1
2->5
5->2
3->5
5->3
4->5
5->4
Sample Output:
1->4
0->5
3->4
2->5
5->2
5->4
5->0
4->1
4->5
4->3

1->5
0->4
3->4
2->5
5->2
5->4
5->1
4->0
4->5
4->3
'''

class NearestNeighbors0:
    def __init__(self):
        #edge, adj = self._input()
        edge, adj = self.readFromFile()
        adj1, adj2 = self.findNearestNeighbors(edge, adj)
        self.printResults(adj1, adj2)
        self.saveResults(adj1, adj2)       
    
    def _input(self):
        data = sys.stdin.read().strip().split('\n')
        edge = [int(i) for i in data[0].split()]
        adj = []
        for d in data[1:]:
            d = [int(i) for i in d.split('->')]
            if max(d) > len(adj)-1:
                adj.extend([[] for _ in range(max(d)-len(adj)+1)])
            adj[d[0]].append(d[1])
        return edge, adj
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        edge = [int(i) for i in data[0].split()]
        adj = []
        for d in data[1:]:
            d = [int(i) for i in d.split('->')]
            if max(d) > len(adj)-1:
                adj.extend([[] for _ in range(max(d)-len(adj)+1)])
            adj[d[0]].append(d[1])
        f.close()
        return edge, adj

    def printResults(self, adj1, adj2):
        for u, e in enumerate(adj1):
            for v in e:
                print(str(u)+'->'+str(v))
        print('')
        for u, e in enumerate(adj2):
            for v in e:
                print(str(u)+'->'+str(v))

    def saveResults(self, adj1, adj2):
        f = open('result.txt', 'w')
        for u, e in enumerate(adj1):
            for v in e:
                f.write(str(u)+'->'+str(v)+'\n')
        f.write('\n')
        for u, e in enumerate(adj2):
            for v in e:
                f.write(str(u)+'->'+str(v)+'\n')
        f.close()
    
    def findNearestNeighbors(self, edge, adj):
        adj1 = deepcopy(adj)
        adj2 = deepcopy(adj)

        adj1[edge[0]].remove(edge[1])
        adj1[edge[1]].remove(edge[0])
        adj1[edge[0]].append(adj1[edge[1]][0])
        adj1[edge[1]].append(adj1[edge[0]][0])
        adj1[adj1[edge[1]][0]].append(edge[0])
        adj1[adj1[edge[0]][0]].append(edge[1])
        adj1[adj1[edge[1]][0]].remove(edge[1])
        adj1[adj1[edge[0]][0]].remove(edge[0])
        del adj1[edge[0]][0]
        del adj1[edge[1]][0]
        adj1[edge[0]].append(edge[1])
        adj1[edge[1]].append(edge[0])

        adj2[edge[0]].remove(edge[1])
        adj2[edge[1]].remove(edge[0])
        adj2[edge[0]].append(adj2[edge[1]][1])
        adj2[edge[1]].append(adj2[edge[0]][0])
        adj2[adj2[edge[1]][1]].append(edge[0])
        adj2[adj2[edge[0]][0]].append(edge[1])
        adj2[adj2[edge[1]][1]].remove(edge[1])
        adj2[adj2[edge[0]][0]].remove(edge[0])
        del adj2[edge[0]][0]
        del adj2[edge[1]][1]
        adj2[edge[0]].append(edge[1])
        adj2[edge[1]].append(edge[0])
        return adj1, adj2

class NearestNeighbors:
    def __init__(self):
        #edge, adj = self._input()
        edge, adj = self.readFromFile()
        adj1, adj2 = self.findNearestNeighbors(edge, adj)
        self.printResults(adj1, adj2)
        self.saveResults(adj1, adj2)       
    
    def _input(self):
        data = sys.stdin.read().strip().split('\n')
        edge = [int(i) for i in data[0].split()]
        adj = []
        for d in data[1:]:
            d = [int(i) for i in d.split('->')]
            if max(d) > len(adj)-1:
                adj.extend([dict() for _ in range(max(d)-len(adj)+1)])
            adj[d[0]][d[1]] = 0
        return edge, adj
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        edge = [int(i) for i in data[0].split()]
        adj = []
        for d in data[1:]:
            d = [int(i) for i in d.split('->')]
            if max(d) > len(adj)-1:
                adj.extend([dict() for _ in range(max(d)-len(adj)+1)])
            adj[d[0]][d[1]] = 0
        f.close()
        return edge, adj

    def printResults(self, adj1, adj2):
        for u, e in enumerate(adj1):
            for v in e.keys():
                print(str(u)+'->'+str(v))
        print('')
        for u, e in enumerate(adj2):
            for v in e.keys():
                print(str(u)+'->'+str(v))

    def saveResults(self, adj1, adj2):
        f = open('result.txt', 'w')
        for u, e in enumerate(adj1):
            for v in e.keys():
                f.write(str(u)+'->'+str(v)+'\n')
        f.write('\n')
        for u, e in enumerate(adj2):
            for v in e.keys():
                f.write(str(u)+'->'+str(v)+'\n')
        f.close()
    
    def findNearestNeighbors(self, edge, adj):
        adj1 = deepcopy(adj)
        adj2 = deepcopy(adj)

        del adj1[edge[0]][edge[1]]
        del adj1[edge[1]][edge[0]]
        e0 = list(adj1[edge[0]].keys())
        e1 = list(adj1[edge[1]].keys())
        adj1[edge[0]][e1[0]] = 0
        adj1[edge[1]][e0[0]] = 0
        adj1[e1[0]][edge[0]] = 0
        adj1[e0[0]][edge[1]] = 0
        del adj1[e1[0]][edge[1]]
        del adj1[e0[0]][edge[0]]
        del adj1[edge[0]][e0[0]]
        del adj1[edge[1]][e1[0]]
        adj1[edge[0]][edge[1]] = 0
        adj1[edge[1]][edge[0]] = 0

        adj2[edge[0]][e1[1]] = 0
        adj2[edge[1]][e0[0]] = 0
        adj2[e1[1]][edge[0]] = 0
        adj2[e0[0]][edge[1]] = 0
        del adj2[e1[1]][edge[1]]
        del adj2[e0[0]][edge[0]]
        del adj2[edge[0]][e0[0]]
        del adj2[edge[1]][e1[1]]
        return adj1, adj2

if __name__ == "__main__":
    NearestNeighbors()