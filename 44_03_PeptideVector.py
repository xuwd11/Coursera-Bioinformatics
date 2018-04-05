# python3

import sys

'''
Given an amino acid string Peptide = a1 . . . an of length n, we will represent its prefix masses using a binary peptide 
vector Peptide' with mass(Peptide) coordinates. This vector contains a 1 at each of the n prefix coordinates
mass(a1), mass(a1 a2), . . . , mass(a1 a2 . . . an ) , and it contains a 0 in each of the remaining noise coordinates. The 
toy peptide XZZXX, whose prefix masses are 4, 9, 14, 18, and 22, corresponds to the peptide vector 
(0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1) of length 22.

1. Converting a Peptide into a Peptide Vector Problem: Convert a peptide into a peptide vector.
     Input: An amino acid string Peptide.
     Output: The peptide vector Peptide'.

Sample Input:
XZZXX
Sample Output:
0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 0 1

2. Converting a Peptide Vector into a Peptide Problem: Convert a peptide vector into a peptide.
     Input: A binary vector P.
     Output: A peptide whose peptide vector is equal to P (if such a peptide exists).

Sample Input:
0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 0 1
Sample Output:
XZZXX
'''

class PeptideVector:
    def __init__(self):
        massDict, aaDict = self.AminoAcidMassDict()
        '''
        peptide = self.readPeptideFromFile()
        peptideVector = self.peptide2vector(peptide, aaDict)
        print(' '.join(map(str, peptideVector)))
        f = open('result.txt', 'w')
        f.write(' '.join(map(str, peptideVector)))
        f.close()
        '''
        peptideVector = self.readVectorFromFile()
        peptide = self.vector2peptide(peptideVector, massDict)
        print(peptide)
        f = open('result.txt', 'w')
        f.write(peptide)
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

    def _inputPeptide(self):
        data = sys.stdin.read().strip()
        peptide = data
        return peptide
    
    def readPeptideFromFile(self):
        f = open('input.txt', 'r')
        for line in f:
            data = line.strip()
        peptide = data
        f.close()
        return peptide
    
    def readVectorFromFile(self):
        f = open('input.txt', 'r')
        for line in f:
            data = line.strip().split()
        peptideVector = list(map(int, data))
        f.close()
        return peptideVector

    def peptide2vector(self, peptide, aaDict):
        l = len(peptide)
        prefixMasses = []
        for i in range(l):
            prefixMasses.append(sum([aaDict[aa] for aa in peptide[:i+1]]))
        peptideVector = [0] * prefixMasses[-1]
        for m in prefixMasses:
            peptideVector[m-1] = 1
        return peptideVector

    def vector2peptide(self, peptideVector, massDict):
        prefixMasses = [i+1 for i, v in enumerate(peptideVector) if 1 == v]
        prefixMasses.insert(0, 0)
        peptide = ''.join([massDict[prefixMasses[i+1]-prefixMasses[i]] for i in range(len(prefixMasses)-1)])
        return peptide

if __name__ == "__main__":
    PeptideVector()