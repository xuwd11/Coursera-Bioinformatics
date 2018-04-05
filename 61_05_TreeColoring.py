# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Tree Coloring Problem: Color the internal nodes of a tree given the colors of its leaves.
     Input: An adjacency list, followed by color labels for leaf nodes.
     Input: Color labels for all nodes, in any order.

Sample Input:
0 -> {}
1 -> {}
2 -> 0,1
3 -> {}
4 -> {}
5 -> 3,2
6 -> {}
7 -> 4,5,6
-
0: red
1: red
3: blue
4: blue
6: red
Sample Output:
0: red
1: red
2: red
3: blue
4: blue
5: purple
6: red
7: purple

TreeColoring, whose pseudocode is shown below and illustrated on the next step, colors the nodes of a suffix tree from the 
leaves upward. This algorithm assumes that the leaves of the suffix tree have been labeled "blue" or "red" and all other nodes 
have been labeled "gray". A node in a tree is called ripe if it is gray but has no gray children.

    TreeColoring(ColoredTree)
        while ColoredTree has ripe nodes
            for each ripe node v in ColoredTree
                if there exist differently colored children of v
                    Color(v) ← "purple"
                else
                    Color(v) ← color of all children of v
        return ColoredTree
'''

class dfs:
    def __init__(self, G):
        self.G = G
    
    def search(self):
        G = self.G
        forward = 1     # traversing edge (v,w) from v to w
        reverse = -1    # returning backwards on (v,w) from w to v
        nontree = 0     # edge (v,w) is not part of the DFS tree
        """
        Generate sequence of triples (v,w,edgetype) for DFS of graph G.
        The subsequence for each root of each tree in the DFS forest starts
        with (root,root,forward) and ends with (root,root,reverse).
        If the initial vertex is given, it is used as the root and vertices
        not reachable from it are not searched.
        """
        visited = set()
    
        for v in range(len(G)):
            if v not in visited:
                yield v,v,forward
                visited.add(v)
                stack = [(v,iter(G[v]))]
                while stack:
                    parent,children = stack[-1]
                    try:
                        child = next(children)
                        if child in visited:
                            yield parent,child,nontree
                        else:
                            yield parent,child,forward
                            visited.add(child)
                            stack.append((child,iter(G[child])))
                    except StopIteration:
                        stack.pop()
                        if stack:
                            yield stack[-1][0],parent,reverse
                yield v,v,reverse
    
    def postorder(self):
        """Generate all vertices of graph G in depth-first postorder."""
        forward = 1     # traversing edge (v,w) from v to w
        reverse = -1    # returning backwards on (v,w) from w to v
        nontree = 0     # edge (v,w) is not part of the DFS tree
        for v,w,edgetype in self.search():
            if edgetype is reverse:
                yield w

class TreeColoring:
    # 0: red
    # 1: blue
    # 2: purple
    def __init__(self):
        #adj, color = self._input()
        adj, color = self.readFromFile()
        #print(adj, color)
        color = self.getColors(adj, color)
        self.saveColors(color)

    def _input(self):
        data = sys.stdin.read().strip()
        data = data.split('\n')
        i = [i for i, d in enumerate(data) if '-' == d][0]
        adj = dict()
        for d in data[:i]:
            d = d.split(' -> ')
            adj[int(d[0])] = []
            if not '{}' in d[1]:
                adj[int(d[0])] = list(map(int, d[1].split(',')))
            else:
                adj[int(d[0])] = []
        color = [-1] * len(adj)
        for d in data[i+1:]:
            d = d.split(': ')
            if 'red' == d[1]:
                color[int(d[0])] = 0
            else:
                color[int(d[0])] = 1
        return adj, color

    def readFromFile(self, fileName = 'input.txt'):
        data = open(fileName, 'r').read().strip().split('\n')
        i = [i for i, d in enumerate(data) if '-' == d][0]
        adj = dict()
        for d in data[:i]:
            d = d.split(' -> ')
            adj[int(d[0])] = []
            if not '{}' in d[1]:
                adj[int(d[0])] = list(map(int, d[1].split(',')))
            else:
                adj[int(d[0])] = []
        color = [-1] * len(adj)
        for d in data[i+1:]:
            d = d.split(': ')
            if 'red' == d[1]:
                color[int(d[0])] = 0
            else:
                color[int(d[0])] = 1
        return adj, color

    def getColors(self, adj, color):
        post = list(dfs(adj).postorder())
        for v in post:
            if -1 == color[v]:
                if len(set([color[w] for w in adj[v]])) == 1:
                    color[v] = color[adj[v][0]]
                else:
                    color[v] = 2
        return color
    
    def saveColors(self, color):
        colorDict = {0: 'red', 1: 'blue', 2: 'purple'}
        f = open('result.txt', 'w')
        for i, c in enumerate(color):
            print('{}: {}'.format(i, colorDict[c]))
            f.write('{}: {}\n'.format(i, colorDict[c]))


if __name__ == "__main__":
    TreeColoring()