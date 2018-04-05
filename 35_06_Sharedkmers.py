# python3

import sys
import numpy as np
import copy

'''
Shared k-mers Problem: Given two strings, find all their shared k-mers.
     Input: An integer k and two strings.
     Output: All k-mers shared by these strings, in the form of ordered pairs (x, y) corresponding to starting positions
     of these k-mers in the respective strings.

Sample Input:
3
AAACTCATC
TTTCAAATC
Sample Output:
(0, 4)
(0, 0)
(4, 2)
(6, 6)
'''

class Sharedkmers:
    def __init__(self):
        k, seqs = self._input()
        #k, seqs = self.readFromFile()
        posList = self.calculateSharedkmers(k, seqs[0], seqs[1])
        print(len(posList))
        for p in posList:
            print(p)
        self.saveResult(posList)

    def _input(self):
        data = sys.stdin.read().strip().split()
        k = int(data[0])
        seqs = data[1:]
        return k, seqs
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        k = int(data[0])
        seqs = data[1:]
        return k, seqs

    def saveResult(self, result):
        f = open('result.txt', 'w')
        for r in result:
            f.write(str(r)+'\n')

    def ReverseComplement(self, seq):
        for base in seq:
            if base not in 'ATCGatcg':
                print("Error: NOT a DNA sequence")
                return None
        seq1 = 'ATCGTAGCatcgtagc'
        seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
        return "".join([seq_dict[base] for base in reversed(seq)])
    
    def calculateSharedkmersNaive(self, k, seq1, seq2):
        posList = []
        for i in range(len(seq1)+1-k):
            for j in range(len(seq2)+1-k):
                if seq1[i:i+k] == seq2[j:j+k] or seq1[i:i+k] == self.ReverseComplement(seq2[j:j+k]):
                    posList.append((i, j))
        return posList
    
    def calculateSharedkmers(self, k, seq1, seq2):
        posList = []
        seq1Dict = dict()
        l = len(seq1)
        m = len(seq2)
        seq1RC = self.ReverseComplement(seq1)
        for i in range(l+1-k):
            s = seq1[i:i+k]
            if not s in seq1Dict:
                seq1Dict[s] = set()
            seq1Dict[s].add(i)
            s = seq1RC[i:i+k]
            if not s in seq1Dict:
                seq1Dict[s] = set()
            seq1Dict[s].add(l-i-k)
        for j in range(m+1-k):
            s = seq2[j:j+k]
            if s in seq1Dict:
                for i in seq1Dict[s]:
                    posList.append((i, j))
        return posList

if __name__ == "__main__":
    Sharedkmers()