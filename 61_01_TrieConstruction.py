# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Trie Construction Problem: Construct a trie from a set of patterns.
     Input: A collection of strings Patterns.
     Output: Trie(Patterns). The adjacency list corresponding to Trie(Patterns), in the following format. If Trie(Patterns) has 
     n nodes, first label the root with 0 and then label the remaining nodes with the integers 1 through n - 1 in any order you 
     like. Each edge of the adjacency list of Trie(Patterns) will be encoded by a triple: the first two members of the triple 
     must be the integers labeling the initial and terminal nodes of the edge, respectively; the third member of the triple must 
     be the symbol labeling the edge.
Sample Input:
ATAGA
ATC
GAT
Sample Output:
0->1:A
1->2:T
2->3:A
3->4:G
4->5:A
2->6:C
0->7:G
7->8:A
8->9:T

TrieConstruction(Patterns)
    Trie ← a graph consisting of a single node root
    for each string Pattern in Patterns
        currentNode ← root
        for i ← 1 to |Pattern|
            if there is an outgoing edge from currentNode with label currentSymbol
                currentNode ← ending node of this edge
            else
                add a new node newNode to Trie
                add a new edge from currentNode to newNode with label currentSymbol
                currentNode ← newNode
    return Trie
'''

class TrieConstruction:
    def __init__(self):
        #patterns = self._input()
        patterns = self.readFromFile()
        trie = self.buildTrie(patterns)
        self.saveResults(trie)
    
    def _input(self):
        return sys.stdin.read().split()
    
    def readFromFile(self, fileName = 'input.txt'):
        f = open(fileName, 'r')
        patterns = []
        for line in f:
            patterns.append(line.strip())
        f.close()
        return patterns
    
    def saveResults(self, trie, fileName = 'result.txt'):
        f = open(fileName, 'w')
        for node in trie:
            for c in trie[node]:
                print('{}->{}:{}'.format(node, trie[node][c], c))
                f.write('{}->{}:{}\n'.format(node, trie[node][c], c))
        f.close()
        return
    
    def buildTrie(self, patterns):
        trie = dict()
        trie[0] = dict()
        newNode = 1
        for pattern in patterns:
            currNode = 0
            for i in range(len(pattern)):
                currSymbol = pattern[i]
                if currSymbol in trie[currNode]:
                    currNode = trie[currNode][currSymbol]
                else:
                    trie[newNode] = dict()
                    trie[currNode][currSymbol] = newNode
                    currNode = newNode
                    newNode += 1
        return trie
                    

if __name__ == "__main__":
    TrieConstruction()