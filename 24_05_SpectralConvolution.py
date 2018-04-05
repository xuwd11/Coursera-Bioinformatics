# python3

import sys

class SpectralConvolution:
    def __init__(self):
        self._input()
        #self.SpectralConvolution(self.spectrum, 1)
        result = self.convolutionList(self.spectrum)
        print(' '.join([str(s) for s in result]))
    
    def _input(self):
        self.spectrum = sorted([int(s) for s in input().strip().split()])
    
    def convolutionList(self, spectrum):
        l = len(spectrum)
        cList = []
        for i in range(l):
            for j in range(i+1, l):
                a = spectrum[j]-spectrum[i]
                if a != 0:
                    cList.append(a)
        return cList
    
    def SpectralConvolution(self, spectrum, M):
        l = len(spectrum)
        cDict = dict()
        for i in range(l):
            for j in range(i+1, l):
                a = spectrum[j] - spectrum[i]
                if 57 <= a and a <= 200:
                    cDict[a] = cDict.get(a, 0) + 1
        sortedMass = sorted(cDict.items(), key = lambda a:a[1], reverse = True)
        #print(sortedMass)
        mass = [str(a[0]) for a in sortedMass]
        multi = [a[1] for a in sortedMass]
        t = multi[M-1]
        for j in range(M, l):
            if multi[j] < t:
                return mass[:j]
        return mass

if __name__ == "__main__":
    SpectralConvolution()