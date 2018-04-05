# python3

import sys

'''
Solve the Eulerian Path Problem.
     Input: The adjacency list of a directed graph that has an Eulerian path.
     Output: An Eulerian path in this graph.

Sample Input:
     0 -> 2
     1 -> 3
     2 -> 1
     3 -> 0,4
     6 -> 3,7
     7 -> 8
     8 -> 9
     9 -> 6

Sample Output:
     6->7->8->9->6->3->0->2->1->3->4

Consider a graph that does not have an Eulerian cycle but does have an Eulerian path. If an Eulerian path in this graph 
connects a node v to a different node w, then the graph is nearly balanced, meaning that all its nodes except v and w are 
balanced. In this case, adding an extra edge from w to v transforms the Eulerian path into an Eulerian cycle. Thus, a nearly 
balanced graph has an Eulerian path if and only if adding an edge between its unbalanced nodes makes the graph balanced and 
strongly connected.
'''

class EulerianPath:
    def __init__(self):
        self.adj = []
        self.n = None
        self.nUnEdges = 0 # number of explored edges
        self.nodesWUE = dict() # key: node with unused edges; value: the position of such node in the current path
        self.inDeg = []
        self.outDeg = []
        self.adjCurPos = []
        self.path = []
        self.unbalancedNode = []
        self._input()
        self.calculateEulerianPath()
        self.printPath()
        self.saveResult()

    def _input(self):
        data = list(sys.stdin.read().strip().split())
        curMax = 0
        for i in range(len(data) // 3):
            curMax = max(int(data[i*3]), curMax, max(list(map(int, data[i*3+2].split(',')))))
        self.n = curMax + 1
        self.adj = [[]] * self.n
        self.unusedEdges = [[]] * self.n
        self.inDeg = [0] * self.n
        self.outDeg = [0] * self.n
        self.adjCurPos = [0] * self.n
        for i in range(len(data) // 3):
            curIn = int(data[i*3])
            self.adj[curIn] = list(map(int, data[i*3+2].split(',')))
            for v in self.adj[curIn]:
                self.inDeg[v] += 1
            l = len(self.adj[curIn])
            self.outDeg[curIn] = l
            self.nUnEdges += l
    
    def addEdge(self):
        for v in range(self.n):
            if self.inDeg[v] != self.outDeg[v]:
                if self.inDeg[v] < self.outDeg[v]:
                    self.unbalancedNode.append(v)
                else:
                    self.unbalancedNode.insert(0, v)
        if len(self.unbalancedNode) > 0:
            self.adj[self.unbalancedNode[0]].append(self.unbalancedNode[1])
            self.outDeg[self.unbalancedNode[0]] += 1
            self.inDeg[self.unbalancedNode[1]] += 1
    
    def explore(self, s):
        self.path.append(s)
        curPos = self.adjCurPos[s]
        curMaxPos = self.outDeg[s]
        while curPos < curMaxPos:
            self.adjCurPos[s] = curPos + 1
            if curPos + 1 < curMaxPos:
                self.nodesWUE[s] = len(self.path) - 1
            else:
                if s in self.nodesWUE:
                    del self.nodesWUE[s]
            v = self.adj[s][curPos]
            self.path.append(v)
            s = v
            curPos = self.adjCurPos[s]
            curMaxPos = self.outDeg[s]
            self.nUnEdges -= 1
        return

    def updatePath(self, startPos):
        l = len(self.path) - 1
        self.path = self.path[startPos:l] + self.path[:startPos]
        for node, pos in self.nodesWUE.items():
            if pos < startPos:
                self.nodesWUE[node] = pos + l - startPos
            else:
                self.nodesWUE[node] = pos - startPos
        return

    def calculateEulerianCycle(self):
        self.explore(0)
        while self.nUnEdges > 0:
            node, pos = self.nodesWUE.popitem()
            self.updatePath(pos)
            self.explore(node)
        return self.path
    
    def calculateEulerianPath(self):
        self.addEdge()
        self.calculateEulerianCycle()
        if len(self.unbalancedNode) > 0:
            for i in range(len(self.path)-1):
                if self.path[i] == self.unbalancedNode[0] and self.path[i+1] == self.unbalancedNode[1]:
                    self.updatePath(i+1)
                    break
        return           

    def printPath(self):
        print('->'.join([str(node) for node in self.path]))     

    def saveResult(self):
        f = open('result.txt', 'w')
        f.write('->'.join([str(node) for node in self.path]))

if __name__ == "__main__":
    EulerianPath()