# python3

import sys

'''
Code Challenge: Implement Trim (reproduced below).
     Input: A collection of peptides Leaderboard, a collection of integers Spectrum, and an integer N.
     Output: The N highest-scoring linear peptides on Leaderboard with respect to Spectrum.


Trim(Leaderboard, Spectrum, N, AminoAcid, AminoAcidMass)
    for j ← 1 to |Leaderboard|
        Peptide ← j-th peptide in Leaderboard
        LinearScores(j) ← LinearScore(Peptide, Spectrum)
    sort Leaderboard according to the decreasing order of scores in LinearScores
    sort LinearScores in decreasing order
    for j ← N + 1 to |Leaderboard|
        if LinearScores(j) < LinearScores(N)
            remove all peptides starting from the j-th peptide from Leaderboard
            return Leaderboard
    return Leaderboard

Sample Input:
LAST ALST TLLT TQAS
0 71 87 101 113 158 184 188 259 271 372
2
Sample Output:
LAST ALST
'''

class TrimLeaderboard:
    def __init__(self):
        #self._inputFromFile('dataset_4913_3.txt')
        self._input()
        result = self.Trim(self.leaderboard, self.spectrumDict, self.N)
        print(' '.join(result))
        

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
        self.leaderboard = input().strip().split()
        self.spectrum = [int(s) for s in input().strip().split()]
        self.N = int(input().strip())
        self.spectrumDict = self.getSpectrumDict(self.spectrum)
    
    def _inputFromFile(self, fileName):
        f = open(fileName, 'r')
        data = []
        for d in f:
            data.append(d.strip())
        self.leaderboard = data[0].split()
        self.spectrum = [int(s) for s in data[1].split()]
        self.N = int(data[2])
        self.spectrumDict = self.getSpectrumDict(self.spectrum)

    def getSpectrumDict(self, spectrum):
        spectrumDict = dict()
        for s in spectrum:
            spectrumDict[s] = spectrumDict.get(s, 0) + 1
        return spectrumDict     

    def calculateLinearScore(self, peptide, spectrumDict):
        theoSpectrumDict = self.LinearSpectrumDict(peptide)
        score = 0
        for s, v in theoSpectrumDict.items():
            v0 = spectrumDict.get(s, 0)
            if v0 >= v:
                score += v
            else:
                score += v0
        return score  

    def Trim(self, leaderboard, spectrumDict, N):
        l = len(leaderboard) 
        linearScoreDict = dict()
        for peptide in leaderboard:
            linearScoreDict[peptide] = self.calculateLinearScore(peptide, spectrumDict)
        lbScore = sorted(linearScoreDict.items(), key = lambda a:a[1], reverse = True)
        leaderboard = [p[0] for p in lbScore]
        linearScores = [p[1] for p in lbScore]
        for j in range(N, l):
            if linearScores[j] < linearScores[N-1]:
                return leaderboard[:j]
        return leaderboard
    
if __name__ == "__main__":
    TrimLeaderboard()