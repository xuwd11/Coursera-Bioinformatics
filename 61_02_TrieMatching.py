# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Implement TrieMatching to solve the Multiple Pattern Matching Problem.
     Input: A string Text and a collection of strings Patterns.
     Output: All starting positions in Text where a string from Patterns appears as a substring.
Sample Input:
AATCGGGTTCAATCGGGGT
ATCG
GGGT
Sample Output:
1 4 11 15

PrefixTrieMatching(Text, Trie)
    symbol ← first letter of Text
    v ← root of Trie
    while forever
        if v is a leaf in Trie
            return the pattern spelled by the path from the root to v
        else if there is an edge (v, w) in Trie labeled by symbol
            symbol ← next letter of Text
            v ← w
        else
            output "no matches found"
            return

PrefixTrieMatching finds whether any strings in Patterns match a prefix of Text. To find whether any strings in Patterns match 
a substring of Text starting at position k, we chop off the first k − 1 symbols from Text and run PrefixTrieMatching on the 
shortened string. As a result, to solve the Multiple Pattern Matching Problem, we simply iterate PrefixTrieMatching |Text| 
times, chopping the first symbol off of Text before each new iteration.

TrieMatching(Text, Trie)
    while Text is nonempty
        PrefixTrieMatching(Text, Trie)
        remove first symbol from Text
'''

class TrieConstruction:
    def __init__(self):
        #patterns = self._input()
        text, patterns = self.readFromFile()
        pos = self.trieMatching(text, patterns)
        self.saveResults(pos)
    
    def _input(self):
        data =  sys.stdin.read().split()
        text = data[0]
        patterns = data[1:]
        return text, patterns
    
    def readFromFile(self, fileName = 'input.txt'):
        f = open(fileName, 'r')
        data = []
        for line in f:
            data.append(line.strip())
        f.close()
        return data[0], data[1:]
    
    def saveResults(self, pos, fileName = 'result.txt'):
        f = open(fileName, 'w')
        print(' '.join(map(str, pos)))
        f.write(' '.join(map(str, pos)))
        f.close()
        return
    
    def buildTrie(self, patterns):
        trie = dict()
        trie[0] = dict()
        newNode = 1
        for pattern in patterns:
            pattern += '$'
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

    def prefixTrieMatching(self, text, trie):
        i = 0
        symbol = text[i]
        l = len(text)
        v = 0
        while True:
            if '$' in trie[v]:
                return True
            elif symbol in trie[v]:
                v = trie[v][symbol]
                if i < l-1:
                    i += 1
                    symbol = text[i]
                elif not '$' in trie[v]:
                    return False
            else:
                return False                

    def trieMatching(self, text, patterns):
        result = []
        trie = self.buildTrie(patterns)
        for i in range(len(text)):
            if self.prefixTrieMatching(text[i:], trie):
                result.append(i)
        return result        

if __name__ == "__main__":
    TrieConstruction()