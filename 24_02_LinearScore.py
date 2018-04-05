# python3

import sys

'''
To implement the Trim function in ﻿LeaderboardCyclopeptideSequencing, we first will generate the theoretical spectra of all linear peptides from Leaderboard. Then, we will compute the scores of each theoretical spectrum against an experimental spectrum Spectrum. This requires implementing LinearScore(Peptide, Spectrum).

Code Challenge: Compute the score of a linear peptide with respect to a spectrum.
     Input: An amino acid string Peptide and a collection of integers Spectrum.
     Output: The linear score of Peptide with respect to Spectrum, LinearScore(Peptide, Spectrum).

Sample Input:
NQEL
0 99 113 114 128 227 257 299 355 356 370 371 484
Sample Output:
8
'''

class LinearScore:
    def __init__(self):
        #self._inputFromFile('dataset_4913_1.txt')
        self._input()
        print(self.calculateScore(self.peptide, self.spectrum))

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
        return {mass[i]:int(mass[i+1]) for i in range(0, len(mass), 2)}

    def CyclicSpectrum(self, peptide):
        massDict = self.AminoAcidMassDict()
        n = len(peptide)
        PrefixMass = [0]
        for i in range(n):
            PrefixMass.append(PrefixMass[i] + massDict[peptide[i]])
        peptideMass = PrefixMass[n]
        cSpectrum = [0]
        for i in range(n):
            for j in range(i+1, n+1):
                cSpectrum.append(PrefixMass[j]-PrefixMass[i])
                if i > 0 and j < n:
                    cSpectrum.append(peptideMass-(PrefixMass[j]-PrefixMass[i]))
        return sorted(cSpectrum)
    
    def CyclicSpectrumDict(self, peptide):
        massDict = self.AminoAcidMassDict()
        n = len(peptide)
        PrefixMass = [0]
        for i in range(n):
            PrefixMass.append(PrefixMass[i] + massDict[peptide[i]])
        peptideMass = PrefixMass[n]
        cSpectrumDict = {0:1}
        for i in range(n):
            for j in range(i+1, n+1):
                s = PrefixMass[j]-PrefixMass[i]
                cSpectrumDict[s] = cSpectrumDict.get(s, 0) + 1
                if i > 0 and j < n:
                    s = peptideMass-(PrefixMass[j]-PrefixMass[i])
                    cSpectrumDict[s] = cSpectrumDict.get(s, 0) + 1
        return(cSpectrumDict)
    
    def LinearSpectrum(self, peptide):
        massDict = self.AminoAcidMassDict()
        n = len(peptide)
        PrefixMass = [0]
        for i in range(n):
            PrefixMass.append(PrefixMass[i] + massDict[peptide[i]])
        lSpectrum = [0]
        for i in range(n):
            for j in range(i+1, n+1):
                lSpectrum.append(PrefixMass[j]-PrefixMass[i])
        return sorted(lSpectrum)

    def LinearSpectrumDict(self, peptide):
        massDict = self.AminoAcidMassDict()
        n = len(peptide)
        PrefixMass = [0]
        for i in range(n):
            PrefixMass.append(PrefixMass[i] + massDict[peptide[i]])
        lSpectrumDict = {0:1}
        for i in range(n):
            for j in range(i+1, n+1):
                s = PrefixMass[j]-PrefixMass[i]
                lSpectrumDict[s] = lSpectrumDict.get(s, 0) + 1
        return lSpectrumDict
    
    def _input(self):
        data = list(sys.stdin.read().strip().split())
        self.peptide = data[0]
        self.spectrum = [int(s) for s in data[1:]]
    
    def _inputFromFile(self, fileName):
        data = open(fileName, 'r').read().strip().split()
        self.peptide = data[0]
        self.spectrum = [int(s) for s in data[1:]]       

    def calculateScore(self, peptide, spectrum):
        theoSpectrumDict = self.LinearSpectrumDict(peptide)
        score = 0
        spectrumDict = dict()
        for s in spectrum:
            spectrumDict[s] = spectrumDict.get(s, 0) + 1
        for s, v in theoSpectrumDict.items():
            v0 = spectrumDict.get(s, 0)
            if v0 >= v:
                score += v
            else:
                score += v0
        return score            
    
if __name__ == "__main__":
    LinearScore()