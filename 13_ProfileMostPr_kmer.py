# python3

import sys
import math
import numpy as np

def Profile(motifs):
    k = len(motifs[0])
    n = len(motifs)
    s = 1 / n
    seq1 = 'ACGTacgt01230123'
    seq_dict = { seq1[i]:int(seq1[i+8]) for i in range(8) }
    P = [[0 for _ in range(k)] for __ in range(4)]
    for motif in motifs:
        for i in range(k):
            P[seq_dict[motif[i]]][i] += s
    return P   

def Pr(pattern, profile):
    seq1 = 'ACGTacgt01230123'
    seq_dict = { seq1[i]:int(seq1[i+8]) for i in range(8) }
    p = 1
    k = len(pattern)
    for i in range(k):
        p *= profile[seq_dict[pattern[i]]][i]
    return p

def ProfileMostPr_kmer(seq, k, profile):
    l = len(seq)
    pmax = -1
    imax = -1
    for i in range(l-k+1):
        p = Pr(seq[i:i+k], profile)
        if p > pmax:
            pmax = p
            imax = i
    return seq[imax:imax+k]

if __name__ == "__main__":
    seq = input().strip()
    data = list(sys.stdin.read().split())
    k = int(data[0])
    profile = []
    for i in range(int((len(data)-1)/k)):
        profile.append(list(map(float, data[1+i*k:2+(i+1)*k])))
    print(ProfileMostPr_kmer(seq, k, profile))