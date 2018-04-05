# python3

import sys
import numpy as np
import copy

'''
Solve the Longest Path in a DAG Problem.
     Input: An integer representing the source node of a graph, followed by an integer representing the sink node of the graph, followed by
     a list of edges in the graph. The edge notation "0->1:7" indicates that an edge connects node 0 to node 1 with weight 7. 
     Output: The length of a longest path in the graph, followed by a longest path. (If multiple longest paths exist, you may return any one.)

Sample Input:
     0
     4
     0->1:7
     0->2:4
     2->3:2
     1->4:1
     3->4:3

Sample Output:
     9
     0->2->3->4
'''

class LongestPath:
    def __init__(self):
        self._input()
        self.sortedNodes = []
        self.topologicalSort(self.nodes, self.outDict, self.inDict)
        self.length, self.backtrack = self.LPBackTrack(self.source, self.sortedNodes, self.inDict)
        self.path = self.OutputLP(self.backtrack, self.source, self.sink)
        print(self.length[self.sink])
        print('->'.join(self.path))

    def _input(self):
        data = sys.stdin.read().strip().split()
        self.source = data[0]
        self.sink = data[1]
        outDict = dict()
        inDict = dict()
        nodes = set()
        for e in data[2:]:
            edge = e.split('->')
            nFrom = edge[0]
            nTo = edge[1].split(':')[0]
            weight = int(edge[1].split(':')[1])
            nodes.add(nFrom)
            nodes.add(nTo)
            if not nFrom in outDict:
                outDict[nFrom] = dict()
                outDict[nFrom][nTo] = weight
            else:
                outDict[nFrom][nTo] = weight
            if not nTo in inDict:
                inDict[nTo] = dict()
                inDict[nTo][nFrom] = weight
            else:
                inDict[nTo][nFrom] = weight
        self.outDict = outDict
        self.inDict = inDict
        self.nodes = nodes
    
    def topologicalSort(self, nodes, outDict, inDict):
        self.unvisited = copy.deepcopy(nodes)
        l = len(nodes)
        self.tempVisited = set()
        while len(self.unvisited) > 0:
            n = list(self.unvisited)[0]
            self.visit(n)

    def visit(self, n):
        if n in self.tempVisited:
            print('not a DAG!')
            return
        if n in self.unvisited:
            self.tempVisited.add(n)
            if n in self.outDict:
                for m in self.outDict[n].keys():
                    self.visit(m)
            self.unvisited.discard(n)
            self.tempVisited.discard(n)
            self.sortedNodes.insert(0, n)

    def LPBackTrack(self, source, nodes, inDict):
        backtrack = dict()
        length = dict()
        isSourceFound = False
        for n in nodes:
            if not isSourceFound:
                if n == source:
                    length[n] = 0
                    isSourceFound = True
                else:
                    continue
            if not n in inDict:
                continue
            else:
                maxLength = -float('inf')
                prevNode = None
                for m, weight in inDict[n].items():
                    if m in length:
                        currLength = length[m] + weight
                        if currLength > maxLength:
                            maxLength = currLength
                            prevNode = m
                if prevNode != None:
                    length[n] = maxLength
                    backtrack[n] = prevNode
        return length, backtrack

    def OutputLP(self, backtrack, source, sink):
        path = [sink]
        n = sink
        while n in backtrack:
            if n == source:
                return path
            n = backtrack[n]
            path.insert(0, n)
        return path

if __name__ == "__main__":
    LongestPath()