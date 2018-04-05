# python3

'''
The theoretical spectrum of a cyclic peptide Peptide, denoted Cyclospectrum(Peptide), is the collection of all of the
masses of its subpeptides, in addition to the mass 0 and the mass of the entire peptide, with masses ordered from 
smallest to largest. We will assume that the theoretical spectrum can contain duplicate elements, as is the case for 
NQEL (shown below), where NQ and EL have the same mass.

Generating Theoretical Spectrum Problem: Generate the theoretical spectrum of a cyclic peptide.
     Input: An amino acid string Peptide.
     Output: Cyclospectrum(Peptide).
Sample Input:
LEQN
Sample Output:
0 113 114 128 129 227 242 242 257 355 356 370 371 484
'''

def AminoAcidMassDict():
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

def LinearSpectrum(peptide):
    massDict = AminoAcidMassDict()
    n = len(peptide)
    PrefixMass = [0]
    for i in range(n):
        PrefixMass.append(PrefixMass[i] + massDict[peptide[i]])
    lSpectrum = [0]
    for i in range(n):
        for j in range(i+1, n+1):
            lSpectrum.append(PrefixMass[j]-PrefixMass[i])
    return sorted(lSpectrum)

def CyclicSpectrum(peptide):
    massDict = AminoAcidMassDict()
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
    
if __name__ == "__main__":
    cSpectrum = LinearSpectrum(input().strip())
    print(' '.join([str(s) for s in cSpectrum]))