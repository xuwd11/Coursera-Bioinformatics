# python3

import math

'''
Counting Peptides with Given Mass Problem: Compute the number of peptides of given mass.
     Input: An integer m.
     Output: The number of linear peptides having integer mass m.

Exercise Break: Solve the Counting Peptides with Given Mass Problem. Recall that we assume that peptides are formed 
from the following 18 amino acid masses:

G	A	S	P	V	T	C	I/L	N	D	K/Q	E	M	H	F	R	Y	W
57	71	87	97	99	101	103	113	114	115	128	129	131	137	147	156	163	186

Sample Input:
1024
Sample Output:
14712706211
'''

def AminoAcidMass():
    '''
    G	A	S	P	V	T	C	I/L	N	D	K/Q	E	M	H	F	R	Y	W
    57	71	87	97	99	101	103	113	114	115	128	129	131	137	147	156	163	186
    '''
    mass = '57 71 87 97 99 101 103 113 114 115 128 129 131 137 147 156 163 186'
    return [int(m) for m in mass.split()]

def CountingPeptides(n):
    massTalbe = AminoAcidMass()
    m = len(massTalbe)
    table = [0] * (n+1)
    table[0] = 1
    for i in range(n+1):
        currSum = 0
        for j in range(m):
            if i - massTalbe[j] >= 0:
                currSum += table[i-massTalbe[j]]
        table[i] += currSum    
    return table[n]
    
if __name__ == "__main__":
    print(CountingPeptides(int(input().strip())))
    '''
    a1 = CountingPeptides(5000)
    a2 = CountingPeptides(5001)
    c = 2**(math.log2(a2) - math.log2(a1))
    print(c)
    '''