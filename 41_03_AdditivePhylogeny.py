# python3

import sys
import queue

'''
Implement AdditivePhylogeny to solve the Distance-Based Phylogeny Problem.
     Input: An integer n followed by a space-separated n x n distance matrix.
     Output: A weighted adjacency list for the simple tree fitting this matrix.

Note on formatting: The adjacency list must have consecutive integer node labels starting from 0. The n leaves must be labeled 0, 1, 
..., n - 1 in order of their appearance in the distance matrix. Labels for internal nodes may be labeled in any order but must start 
from n and increase consecutively.

Sample Input:
4
0	13	21	22
13	0	12	13
21	12	0	13
22	13	13	0
Sample Output:
0->4:11
1->4:2
2->5:6
3->5:7
4->0:11
4->1:2
4->5:4
5->4:4
5->3:7
5->2:6

AdditivePhylogeny(D, n)
    if n = 2
        return the tree consisting of a single edge of length D1,2
    limbLength ← Limb(D, n)
    for j ← 1 to n - 1
        Dj,n ← Dj,n - limbLength
        Dn,j ← Dj,n
    (i, k) ← two leaves such that Di,k = Di,n + Dn,k
    x ← Di,n
    remove row n and column n from D
    T ← AdditivePhylogeny(D, n - 1)
    v ← the (potentially new) node in T at distance x from i on the path between i and k
    add leaf n back to T by creating a limb (v, n) of length limbLength
    return T
'''

class AdditivePhylogeny:
    def __init__(self):
        n, disMatrix = self._input()
        adj = self.reconstructPhylogeny(disMatrix, n)
        self.printGraph(adj)        

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

    def saveResult(self, result):
        f = open('result.txt', 'w')

    def printDistMatrix(self, distMatrix):
        for d in distMatrix:
            print(' '.join([str(i) for i in d]))

    def printGraph(self, adj):
        for i, dicts in enumerate(adj):
            for d, w in dicts.items():
                print(str(i)+'->'+str(d)+':'+str(w))

    def calculateLimbLength(self, distMatrix, n, j):
        limbLength = float('inf')
        if j > 0:
            i = j - 1
        else:
            i = j + 1
        for k in range(n):
            if i != k and k != j:
                currLength = (distMatrix[i][j] + distMatrix[j][k] - distMatrix[i][k])//2
                if currLength < limbLength:
                    limbLength = currLength
                    currIndex = (i, k)
        return limbLength, currIndex[0], currIndex[1]

    def reconstructPhylogeny(self, D, n):
        def addNode(adj, j, limbLength, i, k, x):
            l = len(adj)
            dist = [float('inf')] * l
            parent = [-1] * l
            q = queue.Queue()
            dist[i] = 0
            q.put(i)
            while not q.empty():
                currNode = q.get()
                for node, weight in adj[currNode].items():
                    if float('inf') == dist[node]:
                        dist[node] = dist[currNode] + weight
                        parent[node] = currNode
                        q.put(node)
                        if node == k:
                            prevNode = node
                            while x < dist[prevNode]:
                                currNode = prevNode
                                prevNode = parent[currNode]
                            if x == dist[prevNode]:
                                adj[prevNode][j] = limbLength
                                adj[j][prevNode] = limbLength
                            else:
                                adj.append(dict())
                                newNode = len(adj) - 1
                                adj[j][newNode] = limbLength
                                adj[newNode][j] = limbLength
                                del adj[prevNode][currNode]
                                del adj[currNode][prevNode]
                                adj[prevNode][newNode] = x-dist[prevNode]
                                adj[newNode][prevNode] = x-dist[prevNode]
                                adj[currNode][newNode] = dist[currNode]-x
                                adj[newNode][currNode] = dist[currNode]-x
                            return

        adj = [dict() for _ in range(n)]
        adj[0][1] = D[0][1]
        adj[1][0] = D[1][0]
        for j in range(2, n):
            limbLength, i, k = self.calculateLimbLength(D, j+1, j)
            x = D[i][j] - limbLength
            addNode(adj, j, limbLength, i, k, x)
        return adj

if __name__ == "__main__":
    AdditivePhylogeny()