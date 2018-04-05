# python3

import sys

'''
Implement ConvolutionCyclopeptideSequencing.
     Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.
     Output: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements (and ties) of the convolution of
     Spectrum that fall between 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).

Sample Input:
     20
     60
     57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493

Sample Output:
     99-71-137-57-72-57
'''

class ConvolutionCyclopeptideSequencing:
    def __init__(self):
        self.readSpectrum()
        #self.readSpectrumFromFile('dataset_104_7.txt')
        result = self.calculateSeqList(self.spectrumDict, self.N, version = 0)
        print(' '.join(result))
    
    def readSpectrum(self):
        self.M = int(input().strip())
        self.N = int(input().strip())
        data = [int(s) for s in input().strip().split()]
        self.spectrum = sorted(data)
        self.parentMass = max(data)
        self.spectrumDict = dict()
        for s in data:
            self.spectrumDict[s] = self.spectrumDict.get(s, 0) + 1
        self.mass = self.SpectralConvolution(self.spectrum, self.M)
        return
    
    def readSpectrumFromFile(self, fileName):
        f = open(fileName, 'r')
        data = []
        for d in f:
            data.append(d.strip())
        self.M = int(data[0])
        self.N = int(data[1])
        data = [int(s) for s in data[2].strip().split()]
        self.spectrum = sorted(data)
        self.parentMass = max(data)
        self.spectrumDict = dict()
        for s in data:
            self.spectrumDict[s] = self.spectrumDict.get(s, 0) + 1
        self.mass = self.SpectralConvolution(self.spectrum, self.M)
        return       

    def AminoAcidMass(self, version = 0):
        #G	A	S	P	V	T	C	I/L	N	D	K/Q	E	M	H	F	R	Y	W
        #57	71	87	97	99	101	103	113	114	115	128	129	131	137	147	156	163	186
        if 0 == version:
            mass = '57 71 87 97 99 101 103 113 114 115 128 129 131 137 147 156 163 186'
            return mass.split() #[int(m) for m in mass.split()]
        if 1 == version:
            return [str(aa) for aa in range(57, 201)]
    
    def SpectralConvolution(self, spectrum, M):
        l = len(spectrum)
        cDict = dict()
        for i in range(l):
            for j in range(i+1, l):
                a = spectrum[j] - spectrum[i]
                if 57 <= a and a <= 200:
                    cDict[a] = cDict.get(a, 0) + 1
        sortedMass = sorted(cDict.items(), key = lambda a:a[1], reverse = True)
        mass = [str(a[0]) for a in sortedMass]
        multi = [a[1] for a in sortedMass]
        t = multi[M-1]
        for j in range(M, l):
            if multi[j] < t:
                return mass[:j]
        return mass

    def expand(self, peptides, version = 0):
        #mass = self.AminoAcidMass(version)
        mass = self.mass
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
        if l <= N:
            return leaderboard 
        linearScoreDict = dict()
        for peptide in leaderboard:
            linearScoreDict[peptide] = self.linearScore(peptide, spectrumDict)
        lbScore = sorted(linearScoreDict.items(), key = lambda a:a[1], reverse = True)
        leaderboard = [p[0] for p in lbScore]
        linearScores = [p[1] for p in lbScore]
        t = linearScores[N-1]
        for j in range(N, l):
            if linearScores[j] < t:
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
        print(bestScore)
        return leaderpeptide

if __name__ == "__main__":
    ConvolutionCyclopeptideSequencing()