# python3

import sys
import queue

'''
In this chapter, we define the length of a path in a tree as the sum of the lengths of its edges (rather than the number of edges 
on the path). As a result, the evolutionary distance between two present-day species corresponding to leaves i and j in a tree T is 
equal to the length of the unique path connecting i and j, denoted di,j(T).

Distances Between Leaves Problem: Compute the distances between leaves in a weighted tree.
     Input:  An integer n followed by the adjacency list of a weighted tree with n leaves.
     Output: An n x n matrix (di,j), where di,j is the length of the path between leaves i and j.

Code Challenge: Solve the Distances Between Leaves Problem. The tree is given as an adjacency list of a graph whose leaves are 
integers between 0 and n - 1; the notation a->b:c means that node a is connected to node b by an edge of weight c. The matrix 
you return should be space-separated.

Sample Input:
4
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
Sample Output:
0	13	21	22
13	0	12	13
21	12	0	13
22	13	13	0
'''

class DistanceBetweenLeaves:
    def __init__(self):
        n, adjDict = self._input()
        distMatrix = self.calculateDistMatrix(n, adjDict)
        self.printDistMatrix(distMatrix)

    def _input(self):
        data = sys.stdin.read().strip().split()
        n = int(data[0])
        adjDict = dict()
        for d in data[1:]:
            d = d.split('->')
            d1 = d[1].split(':')
            if not int(d[0]) in adjDict:
                adjDict[int(d[0])] = []
            adjDict[int(d[0])].append((int(d1[0]), int(d1[1])))
        return n, adjDict
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        n = int(data[0])
        adjDict = dict()
        for d in data[1:]:
            d = d.split('->')
            d1 = d[1].split(':')
            if not int(d[0]) in adjDict:
                adjDict[int(d[0])] = []
            adjDict[int(d[0])].append((int(d1[0]), int(d1[1])))
        return n, adjDict

    def saveResult(self, result):
        f = open('result.txt', 'w')

    def printDistMatrix(self, distMatrix):
        for d in distMatrix:
            print(' '.join([str(i) for i in d]))

    def calculateDistMatrix(self, n, adjDict):
        distMatrix = [[0]*n for _ in range(n)]

        def bfs(s):
            dist = dict()
            q = queue.Queue()
            dist[s] = 0
            q.put(s)
            while not q.empty():
                currNode = q.get()
                for node, weight in adjDict[currNode]:
                    if not node in dist:
                        dist[node] = dist[currNode] + weight
                        if node < n:
                            distMatrix[s][node] = dist[node]
                        q.put(node)

        for i in range(n):
            bfs(i)
        
        return distMatrix

if __name__ == "__main__":
    DistanceBetweenLeaves()