# python3

import sys
import queue
import numpy as np
from copy import deepcopy

'''
We represent the masses in a spectrum as a sequence Spectrum of integers  s1,…,sms1,…,sm in increasing order, where s1s1 is 
zero and smsm is the total mass of the (unknown) peptide. We define a labeled graph Graph(Spectrum) by forming a node for each 
element of Spectrum, then connecting nodes sisi and sjsj by a directed edge labeled by an amino acid aa if sj−sisj−si is equal 
to the mass of aa. As we assumed when sequencing antibiotics, we do not distinguish between amino acids having the same integer 
masses (i.e., the pairs K/Q and I/L).

Construct the graph of a spectrum.
     Given: A space-delimited list of integers Spectrum.
     Return: Graph(Spectrum).

Note: Throughout this chapter, all dataset problems implicitly use the standard integer-valued mass table for the regular twenty 
amino acids. Examples sometimes use the toy amino acid alphabet {X, Z} whose masses are 4 and 5, respectively.

Sample Input:
57 71 154 185 301 332 415 429 486
Sample Output:
0->57:G
0->71:A
57->154:P
57->185:K
71->185:N
154->301:F
185->332:F
301->415:N
301->429:K
332->429:P
415->486:A
429->486:G
'''

class SpectrumGraph:
    def __init__(self):
        spectrum = self.readFromFile()
        massDict = self.AminoAcidMassDict()
        adj = self.constructGraph(spectrum, massDict)
        self.printResult(adj, spectrum)
        self.saveResult(adj, spectrum)
    
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
        return {int(mass[i+1]):mass[i] for i in range(0, len(mass), 2)}

    def printResult(self, adj, spectrum):
        for i, aaList in enumerate(adj):
            for j, aa in aaList:
                print(str(spectrum[i])+'->'+str(spectrum[j])+':'+aa)
    
    def saveResult(self, adj, spectrum):
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

if __name__ == "__main__":
    SpectrumGraph()