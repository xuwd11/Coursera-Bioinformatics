# python3

import sys

'''
Cyclopeptide Sequencing Problem (for spectra with errors): Find a cyclic peptide having maximum score against an 
experimental spectrum.
     Input: A collection of integers Spectrum.
     Output: A cyclic peptide Peptide maximizing Score(Peptide, Spectrum) over all peptides Peptide with mass equal 
     to ParentMass(Spectrum).

Code Challenge: Implement LeaderboardCyclopeptideSequencing.
     Input: An integer N and a collection of integers Spectrum.
     Output: LeaderPeptide after running LeaderboardCyclopeptideSequencing(Spectrum, N).
Sample Input:
10
0 71 113 129 147 200 218 260 313 331 347 389 460
Sample Output:
113-147-71-129

Pseudocode:
    CyclopeptideSequencing(Spectrum)
        Peptides ← a set containing only the empty peptide
        while Peptides is nonempty
            Peptides ← Expand(Peptides)
            for each peptide Peptide in Peptides
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Cyclospectrum(Peptide) = Spectrum
                        output Peptide
                    remove Peptide from Peptides
                else if Peptide is not consistent with Spectrum
                    remove Peptide from Peptides
'''

class LeaderboardCyclopeptideSequencing:
    def __init__(self):
        #self.readSpectrum()
        self.readSpectrumFromFile('dataset_102_8.txt')
        result = self.calculateSeqList(self.spectrumDict, self.N, version = 0)
        print(' '.join(result))
    
    def readSpectrum(self):
        self.N = int(input().strip())
        data = [int(s) for s in input().strip().split()]
        self.parentMass = max(data)
        self.spectrumDict = dict()
        for s in data:
            self.spectrumDict[s] = self.spectrumDict.get(s, 0) + 1
        return
    
    def readSpectrumFromFile(self, fileName):
        f = open(fileName, 'r')
        data = []
        for d in f:
            data.append(d.strip())
        self.N = int(data[0])
        data = [int(s) for s in data[1].strip().split()]
        self.parentMass = max(data)
        self.spectrumDict = dict()
        for s in data:
            self.spectrumDict[s] = self.spectrumDict.get(s, 0) + 1
        return       

    def AminoAcidMass(self, version = 0):
        #G	A	S	P	V	T	C	I/L	N	D	K/Q	E	M	H	F	R	Y	W
        #57	71	87	97	99	101	103	113	114	115	128	129	131	137	147	156	163	186
        if 0 == version:
            mass = '57 71 87 97 99 101 103 113 114 115 128 129 131 137 147 156 163 186'
            return mass.split() #[int(m) for m in mass.split()]
        if 1 == version:
            return [str(aa) for aa in range(57, 201)]

    def expand(self, peptides, version = 0):
        mass = self.AminoAcidMass(version)
        expandedPeptides = set()
        for p in peptides:
            if '' == p:
                for m in mass:
                    expandedPeptides.add(m)
            else:
                for m in mass:
                    expandedPeptides.add(p + '-' + m)
        return expandedPeptides
    
    def calculateMass(self, peptide):
        return sum([int(aa) for aa in peptide.split('-')])
    
    def linearSpectrum(self, aaList):
        n = len(aaList)
        PrefixMass = [0]
        for i in range(n):
            PrefixMass.append(PrefixMass[i] + aaList[i])
        lSpectrum = [0]
        for i in range(n):
            for j in range(i+1, n+1):
                lSpectrum.append(PrefixMass[j] - PrefixMass[i])
        currSpectrumDict = dict()
        for s in lSpectrum:
            currSpectrumDict[s] = currSpectrumDict.get(s, 0) + 1
        return currSpectrumDict
    
    def cycloSpectrum(self, aaList):
        n = len(aaList)
        PrefixMass = [0]
        for i in range(n):
            PrefixMass.append(PrefixMass[i] + aaList[i])
        peptideMass = PrefixMass[n]
        cSpectrum = [0]
        for i in range(n):
            for j in range(i+1, n+1):
                cSpectrum.append(PrefixMass[j] - PrefixMass[i])
                if i > 0 and j < n:
                    cSpectrum.append(peptideMass-(PrefixMass[j]-PrefixMass[i]))
        currSpectrumDict = dict()
        for s in cSpectrum:
            currSpectrumDict[s] = currSpectrumDict.get(s, 0) + 1
        return currSpectrumDict
    
    def consistentCycloSpectrum(self, peptide):
        aaList = [int(aa) for aa in peptide.split('-')]
        if self.cycloSpectrum(aaList) == self.spectrumDict:
            return True
        else:
            return False

    def isConsistent(self, peptide):
        aaList = [int(aa) for aa in peptide.split('-')]
        currSpectrumDict = self.linearSpectrum(aaList)
        for key, value in currSpectrumDict.items():
            if value > self.spectrumDict.get(key, 0):
                return False
        return True
    
    def linearScore(self, peptide, spectrumDict):
        if 0 == len(peptide):
            return 0
        aaList = [int(aa) for aa in peptide.split('-')]
        theoSpectrumDict = self.linearSpectrum(aaList)
        score = 0
        for s, v in theoSpectrumDict.items():
            v0 = spectrumDict.get(s, 0)
            if v0 >= v:
                score += v
            else:
                score += v0
        return score
    
    def cycloScore(self, peptide, spectrumDict):
        if 0 == len(peptide):
            return 0
        aaList = [int(aa) for aa in peptide.split('-')]
        theoSpectrumDict = self.cycloSpectrum(aaList)
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
            linearScoreDict[peptide] = self.linearScore(peptide, spectrumDict)
        lbScore = sorted(linearScoreDict.items(), key = lambda a:a[1], reverse = True)
        leaderboard = [p[0] for p in lbScore]
        linearScores = [p[1] for p in lbScore]
        for j in range(N, l):
            if linearScores[j] < linearScores[N-1]:
                return leaderboard[:j]
        return leaderboard
    
    def calculateSequence(self):
        peptides = {''}
        result = []
        while len(peptides) > 0:
            peptides = self.expand(peptides)
            deletions = []
            for peptide in peptides:
                if self.calculateMass(peptide) == self.parentMass:
                    if self.consistentCycloSpectrum(peptide):
                        result.append(peptide)
                    deletions.append(peptide)
                elif not self.isConsistent(peptide):
                    deletions.append(peptide)
            for p in deletions:
                peptides.remove(p)
        return result                    
    
    def calculateSequence2(self, spectrumDict, N):
        leaderboard = {''}
        leaderpeptide = ''
        bestScore = 0
        while len(leaderboard) > 0:
            leaderboard = self.expand(leaderboard)
            deletions = []
            for peptide in leaderboard:
                if self.calculateMass(peptide) == self.parentMass:
                    currScore = self.cycloScore(peptide, spectrumDict)
                    if currScore > bestScore:
                        leaderpeptide = peptide
                        bestScore = currScore
                elif self.calculateMass(peptide) > self.parentMass:
                    deletions.append(peptide)
            for p in deletions:
                leaderboard.remove(p)
            leaderboard = self.Trim(leaderboard, spectrumDict, N)
        #print(bestScore)
        return leaderpeptide

    def calculateSeqList(self, spectrumDict, N, version = 0):
        leaderboard = {''}
        leaderpeptide = ['']
        bestScore = 0
        while len(leaderboard) > 0:
            leaderboard = self.expand(leaderboard, version)
            deletions = []
            for peptide in leaderboard:
                if self.calculateMass(peptide) == self.parentMass:
                    currScore = self.cycloScore(peptide, spectrumDict)
                    if currScore > bestScore:
                        leaderpeptide = [peptide]
                        bestScore = currScore
                    elif currScore == bestScore:
                        leaderpeptide.append(peptide)
                elif self.calculateMass(peptide) > self.parentMass:
                    deletions.append(peptide)
            for p in deletions:
                leaderboard.remove(p)
            leaderboard = self.Trim(leaderboard, spectrumDict, N)
        #print(bestScore)
        return leaderpeptide

if __name__ == "__main__":
    LeaderboardCyclopeptideSequencing()