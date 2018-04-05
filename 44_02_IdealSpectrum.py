# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
Given an amino acid string Peptide, its ideal spectrum, denoted IdealSpectrum(Peptide), is the collection of integer masses of 
all its prefixes and suffixes. Note that an ideal spectrum may have repeated masses; for example, IdealSpectrum(GPG) = {0, 57, 
57, 154, 154, 211}. We say that an amino acid string Peptide explains a collection of integers Spectrum if 
IdealSpectrum(Peptide) = Spectrum.

Decoding an Ideal Spectrum Problem: Reconstruct a peptide from its ideal spectrum.
     Input: A collection of integers Spectrum.
     Output: An amino acid string Peptide that explains Spectrum.

Sample Input:
57 71 154 185 301 332 415 429 486
Sample Output:
GPFNA
'''

class IdealSpectrum:
    def __init__(self):
        spectrum = self.readFromFile()
        peptide = self.decodeSpectrum(spectrum)
        print(peptide)
        f = open('result.txt', 'w')
        f.write(peptide)
        f.close()
    
    def _input(self):
        data = sys.stdin.read().strip()
        spectrum = sorted([int(i) for i in data.split()])
        return spectrum
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        for line in f:
            data = line.strip()
        spectrum = sorted([int(i) for i in data.split()])
        f.close()
        return spectrum

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

    def printGraph(self, adj, spectrum):
        for i, aaList in enumerate(adj):
            for j, aa in aaList:
                print(str(spectrum[i])+'->'+str(spectrum[j])+':'+aa)
    
    def saveGraph(self, adj, spectrum):
        f = open('result.txt', 'w')
        for i, aaList in enumerate(adj):
            for j, aa in aaList:
                f.write(str(spectrum[i])+'->'+str(spectrum[j])+':'+aa+'\n')
        f.close()

    def constructGraph(self, spectrum, massDict):
        adj = [[] for _ in range(len(spectrum))]
        spectrum.insert(0, 0)
        for i in range(len(spectrum)-1):
            for j in range(i+1, len(spectrum)):
                mass = spectrum[j]-spectrum[i]
                if mass in massDict:
                    adj[i].append((j, massDict[mass]))
        return adj

    def LinearSpectrum(self, peptide, aaDict):
        n = len(peptide)
        PrefixMass = [0]
        for i in range(n):
            PrefixMass.append(PrefixMass[i] + aaDict[peptide[i]])
        lSpectrum = [0]
        for i in range(n):
            for j in range(i+1, n+1):
                lSpectrum.append(PrefixMass[j]-PrefixMass[i])
        return sorted(lSpectrum)

    def getIdealSpectrum(self, peptide, aaDict):
        n = len(peptide)
        ispectrum = []
        for i in range(n):
            ispectrum.append(sum([aaDict[aa] for aa in peptide[:i]]))
            ispectrum.append(sum([aaDict[aa] for aa in peptide[i:]]))
        return sorted(ispectrum)

    def decodeSpectrum(self, spectrum):
        massDict, aaDict = self.AminoAcidMassDict()
        adj = self.constructGraph(spectrum, massDict)
        s = 0
        d = adj[-1][-1][0]
        paths = self.findAllPaths(adj, s, d)
        for path in paths:
            ispectrum = self.getIdealSpectrum(path, aaDict)
            if ispectrum == spectrum:
                return ''.join(path)

    def findAllPaths(self, adj, s, d):
        paths = []
        path = []
        self.findAllPathsUtil(adj, '', s, d, path, paths)
        return paths
    
    def findAllPathsUtil(self, adj, char, u, d, path, paths):
        path.append(char)
        if u == d:
            paths.append(path[1:])
        else:
            for v, char in adj[u]:
                self.findAllPathsUtil(adj, char, v, d, path, paths)
        del path[-1]

if __name__ == "__main__":
    IdealSpectrum()