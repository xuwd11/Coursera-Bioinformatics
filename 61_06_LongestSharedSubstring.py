# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Longest Shared Substring Problem: Find the longest substring shared by two strings.
     Input: Strings Text1 and Text2.
     Output: The longest substring that occurs in both Text1 and Text2.

Sample Input:
TCGGTAGATTGCGCCCACTC
AGGGGCTCGCAGTGTAAGAA
Sample Output:
AGA
'''

class Node:
    total = 0
    def __init__(self, isLeaf = True, startIdx = 0, depth = 0):
        Node.total += 1
        self.id = self.total
        self.isLeaf = isLeaf
        self.startIdx = startIdx
        self.depth = depth

class Edge:
    def __init__(self, startIdx, endIdx, text, startNode, endNode, leafLabel = None):
        self.startIdx = startIdx
        self.endIdx = endIdx
        self.text = text
        self.startNode = startNode
        self.endNode = endNode
        self.leafLabel = leafLabel

    def length(self):
        return self.endIdx - self.startIdx + 1
    
    def str(self):
        return self.text[self.startIdx:self.endIdx+1]
    
    def startChar(self):
        return self.text[self.startIdx]    
    
    def __str__(self):
        return self.startChar()

class dfs:
    def __init__(self, G, root):
        self.G = G
        self.root = root
    
    def search(self):
        G = self.G
        root = self.root
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
    
        for v in [root]:
            if v not in visited:
                yield v,v,forward
                visited.add(v)
                stack = [(v,iter([edge.endNode for edge in G[v].values()]))]
                while stack:
                    parent,children = stack[-1]
                    try:
                        child = next(children)
                        if child in visited:
                            yield parent,child,nontree
                        else:
                            yield parent,child,forward
                            visited.add(child)
                            stack.append((child,iter([edge.endNode for edge in G[child].values()])))
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

class SuffixTree: #Naive algorithm
    def __init__(self):
        #text1, text2 = self._input()
        text1, text2 = self.readFromFile()
        self.text = text1+'#'+text2+'$'
        self.l1, self.l2 = len(text1), len(text2)
        self.root = Node()
        self.tree = dict()
        self.build(self.root, self.text)
        lss = self.longestSharedSubstring(self.tree, self.root, self.text, self.l1, self.l2)
        print(lss)
        f = open('result.txt', 'w')
        f.write(lss)
        f.close()
        
    def _input(self):
        data = sys.stdin.read().strip().split()
        return data[0], data[1]

    def readFromFile(self, fileName = 'input.txt'):
        f = open(fileName, 'r')
        data = f.read().strip().split()
        f.close()
        return data[0], data[1]    

    def match(self, i, root, text):
        l = len(text)
        currNode = root
        atNode = True
        for j in range(i, l):
            if atNode:
                currPos = 0
                if not text[j] in self.tree[currNode]:
                    return (currNode, None, j, -1)
                else:
                    currEdge = self.tree[currNode][text[j]]
                    currString = currEdge.str()
                    lrString = len(currString) - 1
                    if lrString == 0:
                        currNode = currEdge.endNode
                        continue
                    else:
                        atNode = False                    
            else:
                currPos += 1
                if text[j] != currString[currPos]:
                    return (currNode, currEdge, j, currEdge.startIdx + currPos)
                else:
                    lrString -= 1
                    if lrString == 0:
                        currNode = currEdge.endNode
                        atNode = True

    def addEdge(self, node, startIdx, endIdx, leafLabel):
        newEndNode = Node(True, leafLabel, endIdx-leafLabel)
        self.tree[newEndNode] = dict()
        newEdge = Edge(startIdx, endIdx, self.text, node, newEndNode, leafLabel)
        self.tree[node][newEdge.startChar()] = newEdge
        node.isLeaf = False


    def splitEdge(self, edge, startIdx, endIdx, cutIdx, leafLabel):
        newNode = Node(False, leafLabel, startIdx-leafLabel)
        newEndNode = Node(True, leafLabel, endIdx-leafLabel)
        newEdge = Edge(startIdx, endIdx, self.text, newNode, newEndNode, leafLabel)
        self.tree[newNode] = dict()
        self.tree[newNode][newEdge.startChar()] = newEdge
        self.tree[newEndNode] = dict()
        edge2 = Edge(cutIdx, edge.endIdx, self.text, newNode, edge.endNode)
        self.tree[newNode][edge2.startChar()] = edge2
        self.tree[edge.startNode][edge.startChar()].endIdx = cutIdx - 1
        self.tree[edge.startNode][edge.startChar()].endNode = newNode

    def build(self, root, text):
        l = len(text)
        node = Node(True, 0, l-1)
        edge1 = Edge(0, l-1, text, root, node)
        self.tree[root] = dict()
        self.tree[root][edge1.startChar()] = edge1
        self.tree[node] = dict()
        for i in range(1, l):
            currNode, currEdge, j, cutIdx = self.match(i, root, text)
            if not currEdge:
                self.addEdge(currNode, j, l-1, i)
            else:
                self.splitEdge(currEdge, j, l-1, cutIdx, i)
    
    def exploreEdges(self, tree):
        results = []
        for node in tree.keys():
            for edge in tree[node].values():
                results.append(edge.str())
        return results
    
    def printTree(self, tree):
        for node in tree.keys():
            print(node.id)
            for edge in tree[node].values():
                if edge.endNode == None:
                    e = None
                else:
                    e = edge.endNode.id
                print(edge.startIdx, edge.endIdx, edge.startNode.id, e, edge.str())
        print('')
    
    def printColors(self, color):
        for v, c in color.items():
            print(v.id, c)
    
    def saveEdges(self, tree):
        f = open('result.txt', 'w')
        f.write('\n'.join(self.exploreEdges(tree)))
        f.close()
    
    def longestSharedSubstring(self, tree, root, text, l1, l2):
        # 0: red; a path ending in it spells out a substring that appears in text2 but not in text1
        # 1: blue; a path ending in it spells out a substring that appears in text1 but not in text2
        # 2: purple; a path ending in a purple node in the suffix tree of Text1 and Text2 spells out a substring shared by 
        #    Text1 and Text2
        post = list(dfs(tree, root).postorder())      
        color = dict()
        for v in post:
            if v.isLeaf:
                color[v] = 0 if v.depth <= l2 else 1
            else:
                ccolors = [color[w.endNode] for w in tree[v].values()]
                if 1 == len(set(ccolors)):
                    color[v] = ccolors[0]
                else:
                    color[v] = 2

        # Find longest shared substring
        best = (0, 0)
        for v in post:
            if 2 == color[v] and v.depth > best[1]:
                best = (v.startIdx, v.depth)
        
        return text[best[0]:best[0]+best[1]]

if __name__ == "__main__":
    SuffixTree()