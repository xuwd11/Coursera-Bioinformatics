# python3

import sys

'''
Solve the Eulerian Cycle Problem.
     Input: The adjacency list of an Eulerian directed graph.
     Output: An Eulerian cycle in this graph.

Sample Input:
     0 -> 3
     1 -> 0
     2 -> 1,6
     3 -> 2
     4 -> 2
     5 -> 4
     6 -> 5,8
     7 -> 9
     8 -> 7
     9 -> 6

Sample Output:
     6->8->7->9->6->5->4->2->1->0->3->2->6

Pseudocode:
    EulerianCycle(Graph)
        form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
        while there are unexplored edges in Graph
            select a node newStart in Cycle with still unexplored edges
            form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking 
            Cycle ← Cycle’
        return Cycle

It may not be obvious, but a good implementation of EulerianCycle will work in linear time. To achieve this runtime speedup, 
you would need to use an efficient data structure in order to maintain the current cycle that Leo is building as well the list 
of unused edges incident to each node and the list of nodes on the current cycle that have unused edges.
'''

class EulerianCycle:
    def __init__(self):
        self.adj = []
        self.n = None
        self.nUnEdges = 0 # number of explored edges
        #self.edgesExplored = []
        self.nodesWUE = dict() # key: node with unused edges; value: the position of such node in the current path
        self.outDeg = []
        self.adjCurPos = []
        self.path = []
        self._input()
        self.calculateEulerianCycle()
        self.printPath()

    def _input(self):
        data = list(sys.stdin.read().strip().split())
        self.n = len(data) // 3
        self.adj = [[]] * self.n
        self.unusedEdges = [[]] * self.n
        self.outDeg = [None] * self.n
        self.adjCurPos = [0] * self.n
        for i in range(self.n):
            curIn = int(data[i*3])
            self.adj[curIn] = list(map(int, data[i*3+2].split(',')))
            l = len(self.adj[curIn])
            self.outDeg[curIn] = l
            #self.edgesExplored[curIn] = [False] * l
            self.nUnEdges += l
    
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

    def printPath(self):
        print('->'.join([str(node) for node in self.path]))       

if __name__ == "__main__":
    EulerianCycle()