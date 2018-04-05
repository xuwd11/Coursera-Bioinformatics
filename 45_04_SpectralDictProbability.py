# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Define Pr(i, t) as the sum of probabilities of all peptides with mass i for which Score(Peptide, Spectrum'i ) is equal to t. 
The set of peptides contributing to Pr(i, t) can be split into 20 subsets depending on their final amino acid. Each peptide 
Peptide ending in a specific amino acid a results in a shorter peptide Peptidea if we remove a; Peptidea has mass i − |a| and 
score t − si. Since the probability of Peptide is 20 times smaller than the probability of Peptidea, the contribution of Peptide 
to Pr(i, t) is 20 times smaller than contribution of Peptidea to Pr(i − |a|, t − si ). Therefore, Pr(i, t) can be computed as

Pr(i,t)=∑all amino acids a120⋅Pr(i−|a|,t−si),

which differs from the recurrence for computing Size(i, t) only in the presence of the factor 1/20. We can now compute the 
probability of a spectral dictionary as 

Pr(Dictionarythreshold(Spectrum'))=∑t≥thresholdPr(m,t).

Solve the Probability of Spectral Dictionary Problem.
     Given: A spectral vector Spectrum', an integer threshold, and an integer max_score.
     Return: The probability of the dictionary Dictionarythreshold(Spectrum').

Note: Use the provided max_score for the height of your table.

Sample Input:
4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 3
1
8
Sample Output:
0.375

Input
-10 11 3 10 11 12 -6 -5 4 4 -2 9 6 -8 9 -6 -1 10 -6 14 4 13 1 -6 5 -7 13 0 -1 12 -2 11 7 -10 9 13 14 -7 7 -9 -6 4 14 2 -9 1 12 13 15 6 15 13 -6 -10 -10 -8 -8 -7 -10 -7 -6 -4 6 9 -6 7 11 -1 -8 1 9 -5 6 7 -3 -10 -9 -1 4 7 7 -6 14 -6 12 15 7 8 11 -5 8 -8 12 -3 -1 -7 -6 9 13 12 -3 7 7 6 3 1 2 4 10 11 -10 -3 14 9 6 8 -9 1 5 -6 -8 5 -7 6 -6 -7 4 1 -3 7 5 10 11 12 0 -10 12 13 11 3 9 8 -10 9 -8 0 15 4 1 1 -4 12 2 4 0 15 -10 4 -10 -10 6 -5 -5 0 10 -5 8 1 14 6 -3 12 9 -7 -4 -9 -9 7 2 6 4 -10 -9 8 -4 -5 0 7 -4 -3 5 12 -10 3 -6 -10 6 10 -6 3 -5 15 4 14 -1 10 -9 13 11 -7 -5 -3 14 15 6 -3 -8 -5 0 12 0 12 2 8 -1 6 2 4 -6 3 11 -4 -10 1 -5 0 14 -5 -6 -1 15 13 12 -10 6 4 0 14 -1 5 15 13 4 -6 13 12 7 14 6 15 10 -9 1 -8 10 9 6 6 2 9 -2 5 11 -4 -6 -10 -7 10 9 8 -6 1 -8 2 -1 -1 -4 -2 0 9 11 -6 9 11 5 5 14 7 -10 14 -4 7 4 14 14 14 8 2 5 14 -4 13 7 10 14 -7 -6 11 -7 -2 -6 -3 1 -7 7 10 15 -6 -2 0 14 1 9 -7 5 -3 -5 5 -5 0 -4 1 3 11 9 -4 -3 -4 0 1 -4 15 -8 -3 0 0 11 -9 11 5 -9 1 -1 -7 -3 8 -9 11 5 4 4 -7 11 -1 -4 -5 7 -7 7 3 6 13 -1 11 -3 13 11 4 3 2 3 0 12 -6 3 12 -10 -8 -9 12 -2 12 5 -3 5 11 5 1 -2 3 5 1 11 6 -6 -2 0 -7 15 14 15 -10 0 6 13 9 10 -2 10 2 8 6 -6 5 -2 1 13 8 14 1 -4 11 11 -8 0 8 5 5 9 -1 -7 3 15 -7 -8 -3 11 9 0 10 2 1 13 4 0 -6 15 15 -1 10 3 1 2
30
200
Output
0.00132187890625
'''

class SpectralDictProb:
    def __init__(self):
        massList = self.AminoAcidMassList()
        sVector, threshold, maxScore = self.readFromFile()
        p = self.dictProb(sVector, threshold, maxScore, massList)
        print(p)
        f = open('result.txt', 'w')
        f.write(str(p))
        f.close()
        
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
    
    def dictProb(self, sVector, threshold, maxScore, massList):
        prob = dict()
        prob[(0, 0)] = 1
        p = sum([self.getProb(len(sVector)-1, t, sVector, massList, prob) for t in range(threshold, maxScore+1)])
        return p
    
    def getProb(self, i, t, sVector, massList, prob):
        if (i, t) in prob:
            return prob[(i, t)]
        if i < 0 or t < 0:
            prob[(i, t)] = 0
            return 0
        p = sum([self.getProb(i-m, t-sVector[i], sVector, massList, prob)/20 for m in massList])
        prob[(i, t)] = p
        return p

if __name__ == "__main__":
    SpectralDictProb()