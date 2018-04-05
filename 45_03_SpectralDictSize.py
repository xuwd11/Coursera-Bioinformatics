# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Size of Spectral Dictionary Problem: Find the size of a spectral dictionary for a given spectrum and score threshold.
     Input: A spectral vector Spectrum' and an integer threshold.
     Output: The number of peptides in Dictionarythreshold(Spectrum').

We will use dynamic programming to solve the Size of Spectral Dictionary Problem. Given a spectral vector Spectrum' = 
(s1, . . . , sm), we define its i-prefix (for i between 1 and m) as Spectrum'i = (s1, . . . , si ) and introduce a variable 
Size(i, t) as the number of peptides Peptide of mass i such that Score(Peptide, Spectrum'i) is equal to t.

The key to establishing a recurrence relation for computing Size(i, t) is to realize that the set of peptides contributing 
to Size(i, t) can be split into 20 subsets depending on their final amino acid a. Each peptide ending in a specific amino 
acid a results in a shorter peptide with mass i − |a| and score t − si   if we remove a from the peptide (here, |a| denotes 
the mass of a). Thus,

Size(i,t)=∑all amino acids aSize(i−|a|,t−si).

Since there is a single “empty” peptide of length zero, we initialize Size(0, 0) = 1. We also define Size(0, t) = 0 for all 
possible scores t, and set Size(i, t) = 0 for negative values of i. Using the above recurrence, we can compute the size of a 
spectral dictionary of Spectrum' = (s1, . . . , sm) as

|Dictionarythreshold(Spectrum)|=∑t≥thresholdSize(m,t).

Solve the Size of Spectral Dictionary Problem.
     Given: A spectral vector Spectrum', an integer threshold, and an integer max_score.
     Return: The size of the dictionary Dictionarythreshold(Spectrum').

Note: Use the 20 amino acid alphabet as well as the provided max_score for the height of your table. Your answer should be the 
number of peptides whose score is at least T and at most max_score.

Sample Input:
4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 3
1
8
Sample Output:
3

Input
14 -4 -3 -3 5 9 0 14 2 1 -4 6 -1 13 2 -5 13 -8 -8 3 0 -10 14 4 14 14 8 -8 1 3 -10 -2 2 -9 3 6 13 -10 6 -8 12 2 8 -1 -5 -6 -6 10 3 -3 12 -4 14 3 11 14 15 12 -7 -5 -2 11 13 -9 15 -8 -10 5 -8 5 6 -9 -2 7 -6 -1 -2 12 12 -3 9 0 3 0 5 -6 3 3 -7 6 0 -6 5 8 -7 5 3 13 13 -2 -9 0 2 13 13 12 7 -2 -10 -5 -7 7 13 11 14 -4 -9 15 -10 5 -7 -6 -7 -6 11 5 9 8 -4 7 1 -9 12 2 8 12 -6 0 2 -5 10 11 14 15 -1 3 -3 3 -3 12 15 4 -2 14 13 8 -10 2 -3 0 -6 8 3 10 0 9 10 13 15 6 9 -10 -9 1 -3 -10 8 1 -10 2 -1 14 -3 15 -1 0 1 6 -7 5 12 6 -9 2 1 -2 14 -5 1 -8 -6 11 -5 2 -3 -8 7 -6 -10 8 6 13 -8 -5 -10 12 -5 -9 8 9 0 10 15 -1 4 2 -8 9 1 -9 -6 -8 -1 -1 5 10 -4 7 3 11 4 12 6 6 13 -3 12 -3 1 7 11 6 13 8 3 -6 5 11 4 -1 15 10 -8 -7 0 4 7 5 -4 8 -3 -4 -8 9 -2 -3 13 1 12 4 -1 13 -1 -5 -5 7 7 -7 -5 6 6 -2 -5 7 10 14 11 12 -9 6 -3 4 15 -8 11 -3 -7 5 -4 7 9 15 -9 8 13 6 -2 -3 9 6 5 14 10 -7 -9 -8 10 2 -3 -1 2 3 12 13 6 -2 8 -5 5 -3 -8 10 3 0 12 -7 10 6 15 8 7 -2 8 14 -2 13 -1 8 15 -7 -7 -7 7 -3 -2 5 -4 -3 15 11 -4 9 11 13 15 8 4 -6 7 12 14 6 -10 -5 -9 4 -9 13 -3 0 12 3 12 -5 11 1 15 -8 5 3 -5 7 15 -2 -9 0 0 1 1 -1 -4 -1 5 12 12 -5 8 5 14 12 5 -9 2 -10 -9 4 -2 6 5 -3 -7 7 5 -8 -10 8 -1 7 3 6 -6 14 -8 6 -5 -8 -10 14 -2 12 4 5 -2 9 -4 1 -5 -3 -6 -8 -9 -10 10 4 9 11 -7 6 -4 5 13 -8 -7 -3 10 14 4 10 6 4 0 13 -3 11 -9 2 -8 6 -8 4 -1
37
200
Output
330
'''

class SpectralDictSize:
    def __init__(self):
        massList = self.AminoAcidMassList()
        sVector, threshold, maxScore = self.readFromFile()
        s = self.dictSize(sVector, threshold, maxScore, massList)
        print(s)
        
    def AminoAcidMassList(self):
        massTable = '''
G 57
A 71
S 87
P 97
V 99
T 101
C 103
I 113
L 113
N 114
D 115
K 128
Q 128
E 129
M 131
H 137
F 147
R 156
Y 163
W 186'''
        mass = massTable.split()
        return [int(mass[i+1]) for i in range(0, len(mass), 2)]

    def readFromFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip().split())
        sVector = [0] + list(map(int, data[0]))
        threshold = int(data[1][0])
        maxScore = int(data[2][0])
        return sVector, threshold, maxScore
    
    def dictSize(self, sVector, threshold, maxScore, massList):
        size = dict()
        size[(0, 0)] = 1
        s = sum([self.getSize(len(sVector)-1, t, sVector, massList, size) for t in range(threshold, maxScore+1)])
        return s
    
    def getSize(self, i, t, sVector, massList, size):
        if (i, t) in size:
            return size[(i, t)]
        if i < 0 or t < 0:
            size[(i, t)] = 0
            return 0
        s = sum([self.getSize(i-m, t-sVector[i], sVector, massList, size) for m in massList])
        size[(i, t)] = s
        return s

if __name__ == "__main__":
    SpectralDictSize()