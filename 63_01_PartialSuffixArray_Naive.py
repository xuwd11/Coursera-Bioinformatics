# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Construct a partial suffix array.
     Input: A string Text and a positive integer K.
     Output: SuffixArrayK(Text), in the form of a list of ordered pairs (i, SuffixArray(i)) for all nonempty entries in the 
     partial suffix array.

Sample Input:
PANAMABANANAS$
5
Sample Output:
1,5
11,10
12,0
'''

class SuffixArray: # O(|S|log|S|+|Sigma|)
    '''
    Build suffix array of the string text and
    return a list result of the same length as the text
    such that the value result[i] is the index (0-based)
    in text where the i-th lexicographically smallest
    suffix of text starts.
    '''
    def __init__(self, text):
        #text = self._input()
        #text = self.readFromFile()
        self.order = self.buildSuffixArray(text)
        #print(', '.join(map(str, order)))
        #f = open('result.txt', 'w')
        #f.write(', '.join(map(str, order)))
    
    def _input(self):
        return sys.stdin.readline().strip()
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        text = f.read().strip()
        f.close()
        return text

    def sortCharacters(self, S):
        l = len(S)
        order = [0] * l
        count = dict()
        for i in range(l):
            count[S[i]] = count.get(S[i], 0) + 1
        charList = sorted(count.keys())
        prevChar = charList[0]
        for char in charList[1:]:
            count[char] += count[prevChar]
            prevChar = char
        for i in range(l-1, -1, -1):
            c = S[i]
            count[c] = count[c] - 1
            order[count[c]] = i
        return order

    def computeCharClasses(self, S, order):
        l = len(S)
        charClass = [0] * l
        charClass[order[0]] = 0
        for i in range(1, l):
            if S[order[i]] != S[order[i-1]]:
                charClass[order[i]] = charClass[order[i-1]] + 1
            else:
                charClass[order[i]] = charClass[order[i-1]]
        return charClass        
    
    def sortDoubled(self, S, L, order, _class):
        sLen = len(S)
        count = [0] * sLen
        newOrder = [0] * sLen
        for i in range(sLen):
            count[_class[i]] += 1
        for j in range(1, sLen):
            count[j] += count[j-1]
        for i in range(sLen-1, -1, -1):
            start = (order[i]-L+sLen) % sLen
            cl = _class[start]
            count[cl] -= 1
            newOrder[count[cl]] = start
        return newOrder
    
    def updateClasses(self, newOrder, _class, L):
        n = len(newOrder)
        newClass = [0] * n
        newClass[newOrder[0]] = 0
        for i in range(1, n):
            curr = newOrder[i]
            prev = newOrder[i-1]
            mid = curr + L
            midPrev = (prev + L) % n
            if _class[curr] != _class[prev] or _class[mid] != _class[midPrev]:
                newClass[curr] = newClass[prev] + 1
            else:
                newClass[curr] = newClass[prev]
        return newClass
    
    def buildSuffixArray(self, S):
        sLen = len(S)
        order = self.sortCharacters(S)
        _class = self.computeCharClasses(S, order)
        L = 1
        while L < sLen:
            order = self.sortDoubled(S, L, order, _class)
            _class = self.updateClasses(order, _class, L)
            L = 2 * L
        return order

    def getSuffixArray(self):
        return self.order            

class PartialSuffixArray:
    def __init__(self):
        text, K = self.readFromFile()
        partial = self.buildPartialSuffixArray(text, K)
        print('\n'.join([','.join(map(str, partial[i])) for i in range(len(partial))]))
        f = open('result.txt', 'w')
        f.write('\n'.join([','.join(map(str, partial[i])) for i in range(len(partial))]))
        f.close()

    def _input(self):
        data = sys.stdin.readline().strip()
        return data[0], int(data[1])
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().strip().split()
        f.close()
        return data[0], int(data[1])
    
    def buildPartialSuffixArray(self, text, K):
        order = SuffixArray(text).getSuffixArray()
        partial = [(i, order[i]) for i in range(len(order)) if 0 == order[i]%K]
        return partial

if __name__ == '__main__':
    PartialSuffixArray()