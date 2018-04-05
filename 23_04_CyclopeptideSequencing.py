# python3

import sys

'''
Cyclopeptide Sequencing Problem: Given an ideal spectrum, find a cyclic peptide whose theoretical spectrum matches 
the experimental spectrum.
     Input: A collection of (possibly repeated) integers Spectrum corresponding to an ideal
     spectrum.
     Output: An amino acid string Peptide such that Cyclospectrum(Peptide) = Spectrum (if such a
     string exists).
Sample Input:
0 113 128 186 241 299 314 427
Sample Output:
186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186

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

class CyclopeptideSequencing:
    def __init__(self):
        self.readSpectrum()
        print(' '.join(self.calculateSequence()))
    
    def readSpectrum(self):
        data = [int(s) for s in list(sys.stdin.read().strip().split())]
        self.parentMass = max(data)
        self.spectrumDict = dict()
        for s in data:
            self.spectrumDict[s] = self.spectrumDict.get(s, 0) + 1
        return

    def AminoAcidMass(self):
        #G	A	S	P	V	T	C	I/L	N	D	K/Q	E	M	H	F	R	Y	W
        #57	71	87	97	99	101	103	113	114	115	128	129	131	137	147	156	163	186
        mass = '57 71 87 97 99 101 103 113 114 115 128 129 131 137 147 156 163 186'
        return mass.split() #[int(m) for m in mass.split()]
    
    def expand(self, peptides):
        mass = self.AminoAcidMass()
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
    
if __name__ == "__main__":
    CyclopeptideSequencing()