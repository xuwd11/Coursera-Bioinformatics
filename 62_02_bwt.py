# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Burrows-Wheeler Transform Construction Problem: Construct the Burrows-Wheeler transform of a string.
     Input: A string Text.
     Output: BWT(Text).

Sample Input:
GCGTGCCTGGTCA$
Sample Output:
ACTGGCT$TGCGGC

Ref:
https://www.quora.com/How-can-I-optimize-burrows-wheeler-transform-and-inverse-transform-to-work-in-O-n-time-O-n-space
'''

class SuffixArray:
    '''
    Build suffix array of the string text and
    return a list result of the same length as the text
    such that the value result[i] is the index (0-based)
    in text where the i-th lexicographically smallest
    suffix of text starts.
    '''
    def __init__(self, text):
        self.sa = self.buildSuffixArray(text)        
    
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
        return self.sa

class BWT:
    def __init__(self):
        #text = self._input()
        text = self.readFromFile()
        #print(self.bwtNaive(text))
        bwt = self.bwtFromSuffixArray(text)
        print(bwt)
        f = open('result.txt', 'w')
        f.write(bwt)
        f.close()

    def _input(self):
        return sys.stdin.readline().strip()

    def readFromFile(self):
        f = open('input.txt', 'r')
        text = f.read().strip()
        f.close()
        return text

    def bwtNaive(self, text):
        return ''.join([m[-1] for m in sorted([text[i:]+text[:i] for i in range(len(text))])])   

    def bwtFromSuffixArray(self, text):
        sa = SuffixArray(text).getSuffixArray()
        return ''.join([text[(sa[i]+len(text)-1)%len(text)] for i in range(len(text))])

if __name__ == '__main__':
    BWT()