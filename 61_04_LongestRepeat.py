# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Longest Repeat Problem: Find the longest repeat in a string.
     Input: A string Text.
     Output: A longest substring of Text that appears in Text more than once.
     (Multiple solutions may exist, in which case you may return any one.)

Sample Input:
ATATCGTTTTATCGTT
Sample Output:
TATCGTT
'''

class Node:
    total = 0
    def __init__(self, startIdx = 0, depth = 0):
        Node.total += 1
        self.id = self.total
        self.startIdx = startIdx
        self.depth = depth

class Edge:
    def __init__(self, startIdx, endIdx, text, startNode, endNode = None, leafLabel = None):
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


class SuffixTree: #Naive algorithm
    def __init__(self, text):
        self.root = Node()
        self.text = text
        self.tree = dict()
        self.build(self.root, self.text)
        lr = self.longestRepeat(self.tree, self.text)
        print(lr)
        f = open('result.txt', 'w')
        f.write(lr)
        f.close()

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
        newEdge = Edge(startIdx, endIdx, self.text, node, None, leafLabel)
        self.tree[node][newEdge.startChar()] = newEdge

    def splitEdge(self, edge, startIdx, endIdx, cutIdx, leafLabel):
        newNode = Node(leafLabel, startIdx-leafLabel)
        newEdge = Edge(startIdx, endIdx, self.text, newNode, None, leafLabel)
        self.tree[newNode] = dict()
        self.tree[newNode][newEdge.startChar()] = newEdge
        edge2 = Edge(cutIdx, edge.endIdx, self.text, newNode, edge.endNode)
        self.tree[newNode][edge2.startChar()] = edge2
        self.tree[edge.startNode][edge.startChar()].endIdx = cutIdx - 1
        self.tree[edge.startNode][edge.startChar()].endNode = newNode

    def build(self, root, text):
        l = len(text)
        edge1 = Edge(0, l-1, text, root)
        self.tree[root] = dict()
        self.tree[root][edge1.startChar()] = edge1
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
    
    def saveEdges(self, tree):
        f = open('result.txt', 'w')
        f.write('\n'.join(self.exploreEdges(tree)))
        f.close()
    
    def longestRepeat(self, tree, text):
        best = (0, 0)
        for node in tree.keys():
            if node.depth > best[1]:
                best = (node.startIdx, node.depth)
        return text[best[0]:best[0]+best[1]]

if __name__ == "__main__":
    f = open('input.txt', 'r')
    for line in f:
        data = line
    f.close()
    SuffixTree(data+'$')