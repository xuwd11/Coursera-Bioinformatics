# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
We will use the term block indel to refer to the addition or removal of a block of consecutive zeroes from a binary vector. 
Thus, applying k modifications to an amino ￼acid string Peptide corresponds to applying k block indels to its peptide vector 
Peptide'. We define Variantsk(Peptide) as the set of all modified variants of Peptide with up to k modifications.

Given a peptide Peptide and a spectral vector Spectrum', our goal is to find a modified peptide from Variantsk(Peptide) with 
maximum score against Spectrum'.

Spectral Alignment Problem: Given a peptide and a spectral vector, find a modified variant of this peptide that maximizes the 
peptide-spectrum score, among all variants of the peptide with up to k modifications.
     Input: An amino acid string Peptide, a spectral vector Spectrum', and an integer k.
     Output: A peptide of maximum score against Spectrum' among all peptides in Variantsk(Peptide).

Sample Input:
XXZ
4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 -1
2
Sample Output:
XX(-1)Z(+2)

Input
LVWSTE
-3 1 12 1 5 12 -4 10 0 9 -4 -9 -6 14 -10 8 11 -3 3 2 11 12 -10 -4 15 -5 2 1 4 -10 -10 7 9 10 8 -9 15 13 -9 7 7 1 13 15 0 -3 -2 7 6 6 7 15 -5 10 7 -4 3 11 14 10 -10 2 -8 12 0 0 -6 15 1 14 3 12 4 11 -4 -7 5 3 8 -10 3 7 11 8 5 15 7 13 3 1 -1 -1 9 5 7 3 9 10 3 -9 10 6 13 8 5 2 -5 14 3 14 12 -7 11 -3 15 11 -8 10 -4 14 -4 5 -2 8 -3 2 -7 12 -10 0 -3 14 0 6 14 0 14 7 -5 3 -3 -1 -3 7 8 -4 2 -6 10 5 13 -9 -3 10 1 4 0 11 12 5 -10 -6 1 -9 1 13 -3 7 -5 6 0 5 3 -9 3 -4 1 -6 -4 11 -10 8 -2 15 6 5 -7 12 3 11 5 15 -8 6 -4 8 10 -7 12 -3 -4 10 11 -4 -3 9 0 -5 1 -6 11 7 -3 -7 -9 13 10 -7 -6 14 15 -10 14 15 -8 5 1 10 -5 2 -7 4 14 10 7 -1 1 12 -7 7 -8 -1 12 3 -9 4 -7 6 8 10 8 5 6 12 4 4 4 12 -4 6 7 -10 7 -7 0 -1 -2 6 -10 -5 -4 5 -3 -1 1 2 10 -9 -6 14 9 -3 -8 -2 -4 -1 0 -10 6 10 -6 -3 4 -1 -7 0 6 3 -3 5 4 -4 0 3 0 -1 -1 -4 -10 2 -3 -9 -10 5 7 10 5 2 -8 2 -6 13 6 -7 5 5 -5 15 8 11 12 4 8 5 11 3 -7 14 3 8 0 2 -4 6 -7 1 8 13 -1 14 3 -9 12 -3 12 14 0 -2 11 11 13 3 15 -9 -5 2 1 -10 7 -8 -10 -6 -7 1 0 7 9 5 -9 5 0 -3 -2 13 14 -10 4 3 5 13 -3 -4 -5 11 14 5 1 -10 -1 4 -6 13 -6 1 -5 2 -2 3 12 0 -7 15 10 -10 11 7 0 9 -7 -5 9 2 -8 -9 -5 -6 13 7 11 -3 10 1 -2 0 6 -2 -4 5 5 -4 15 11 4 -7 10 -1 9 -2 5 13 -10 14 7 6 4 -8 1 12 12 2 5 13 7 6 -4 -2 7 4 -8 -7 8 14 -8 14 -2 -5 -9 -5 11 10 7 3 1 11 1 14 4 -9 14 11 12 -1 -3 -7 10 13 -7 4 6 8 2 1 -3 -3 -9 11 7 -1 5 -5 -6 -9 2 11 8 13 1 -6 7 5 3 2 2 6 12 -9 12 3 -3 -1 -10 -4 -4 10 10 14 12 -2 7 14 9 3 10 -5 5 -7 4 9 1 13 10 13 13 -8 6 4 10 -4 15 -9 -1 14 9 -10 3 15 0 -3 2 4 -5 1 8 3 15 12 3 0 9 5 1 3 5 14 -6 13 0 13 -7 4 -9 -2 -9 -3 8 8 13 13 8 4 -2 0 5 -9 -2 0 11 7 0 -2 -6 -6 11 2 -3 7 1 5 -10 7 -10 0 12 12 -2 7 -3 2 -6 13 8 -4 11 -10 -9 4 9 -8 14 -3 7 -2 -4 8 4 -3 -1 -2 -5 12 -9 8 15 10 10 13 10 6 4 -8 7 -9 -1 7 -2 14 -10 2 6 -9 6 9 0 6 15 10 2 6 12 10 9 -2 5 -7 4 -3 9 -6 6 1 6 14 13 14 3 4 9 -8 0 10 7 -6 3 11 7 9 -8 8 5 -8 8 4 10 0
3
Output
L(-61)VW(-9)STE(+69)
'''

class SpectralAlignment:
    def __init__(self):
        massDict, aaDict = self.AminoAcidMassDict()
        peptide, sVector, k = self.readFromFile()
        mPeptide = self.constructAlignment(peptide, sVector, k, aaDict)
        print(mPeptide)
        f = open('result.txt', 'w')
        f.write(mPeptide)
        f.close()
    
    def AminoAcidMassDict(self):
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
        return {int(mass[i+1]):mass[i] for i in range(0, len(mass), 2)}, {mass[i]:int(mass[i+1]) for i in range(0, len(mass), 2)}
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip().split())
        peptide = data[0][0]
        sVector = [0] + list(map(int, data[1]))
        k = int(data[2][0])
        return peptide, sVector, k
    
    def getScore(self, node, score):
        if node in score:
            return score[node]
        else:
            score[node] = -np.inf
            return -np.inf

    def constructAlignment(self, peptide, sVector, k, aaDict):
        prefixMasses = [0] + [sum([aaDict[aa] for aa in peptide[:i+1]]) for i in range(len(peptide))]
        diff = {prefixMasses[i+1]:(prefixMasses[i+1]-prefixMasses[i]) for i in range(len(peptide))}
        score = dict()
        score[(0,0,0)] = 0
        backtrack = dict()
        for i in prefixMasses[1:]:
            if i < len(sVector):
                score[(i,i,0)] = sVector[i]+score[(i-diff[i],i-diff[i],0)]
                backtrack[(i,i,0)] = (i-diff[i],i-diff[i],0)
            else:
                break
        for t in range(1, k+1):
            for i in prefixMasses[1:]:
                for j in range(1, len(sVector)):
                    prevList = [(i-diff[i],j-diff[i],t)] + [(i-diff[i],j1,t-1) for j1 in range(j)]
                    prevIndex = np.argmax([sVector[j] + self.getScore(node, score) for node in prevList])
                    score[(i, j, t)] = sVector[j] + self.getScore(prevList[prevIndex], score)
                    backtrack[(i, j, t)] = prevList[prevIndex]
        lastNodes = [(prefixMasses[-1], len(sVector)-1, t) for t in range(k+1)]
        t = np.argmax([self.getScore(node, score) for node in lastNodes])
        prevNode = lastNodes[t]
        mPeptide = ''
        pList = [peptide[i] for i in range(len(peptide))]
        while (0,0,0) != prevNode:
            node = backtrack[prevNode]
            if node[2] == prevNode[2]:
                mPeptide = pList.pop() + mPeptide
                prevNode = node
            else:
                indel = prevNode[1]-node[1]-diff[prevNode[0]]
                if indel > 0:
                    mPeptide = '(+'+str(indel)+')'+mPeptide
                else:
                    mPeptide = '('+str(indel)+')'+mPeptide 
                mPeptide = pList.pop() + mPeptide
                prevNode = node
        return mPeptide        

if __name__ == "__main__":
    SpectralAlignment()