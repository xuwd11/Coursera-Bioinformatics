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

def PseudoProfile(motifs):
    k = len(motifs[0])
    n = len(motifs)
    s = 1 / (n + 4)
    seq1 = 'ACGTacgt01230123'
    seq_dict = { seq1[i]:int(seq1[i+8]) for i in range(8) }
    P = [[1 for _ in range(k)] for __ in range(4)]
    for motif in motifs:
        for i in range(k):
            P[seq_dict[motif[i]]][i] += s
    return P


def Consensus(motifs):
    P = Profile(motifs)
    k = len(P[0])
    Pt = [[row[i] for row in P] for i in range(k)]
    seq_dict = ['A', 'C', 'G', 'T']
    return ''.join([seq_dict[np.argmax(Pt[i])] for i in range(k)])

def Score(motifs):
    k = len(motifs[0])
    n = len(motifs)
    seq1 = 'ACGTacgt01230123'
    seq_dict = { seq1[i]:int(seq1[i+8]) for i in range(8) }
    P = [[0 for _ in range(4)] for __ in range(k)]
    for motif in motifs:
        for i in range(k):
            P[i][seq_dict[motif[i]]] += 1
    Sm = 0
    for i in range(k):
        Sm += max(P[i])
    return n * k - Sm

def Count(motifs):
    k = len(motifs[0])
    n = len(motifs)
    seq1 = 'ACGTacgt01230123'
    seq_dict = { seq1[i]:int(seq1[i+8]) for i in range(8) }
    P = [[0 for _ in range(k)] for __ in range(4)]
    for motif in motifs:
        for i in range(k):
            P[seq_dict[motif[i]]][i] += 1
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

def GreedyMotifSearch(dna, k, t):
    BestMotifs = [dna[i][0:k] for i in range(t)]
    BestScore = float('inf')
    dna1 = dna[0]
    l1 = len(dna1)
    for i in range(l1-k+1):
        motifs = []
        motifs.append(dna1[i:i+k])
        for i in range(1, t):
            P = Profile(motifs)
            motifs.append(ProfileMostPr_kmer(dna[i], k, P))
        currScore = Score(motifs)
        if currScore < BestScore:
            BestMotifs = motifs
            BestScore = currScore
    return BestMotifs

def GreedyMotifSearch2(dna, k, t):
    #GreedyMotifSearch with pseudocounts
    BestMotifs = [dna[i][0:k] for i in range(t)]
    BestScore = float('inf')
    dna1 = dna[0]
    l1 = len(dna1)
    for i in range(l1-k+1):
        motifs = []
        motifs.append(dna1[i:i+k])
        for i in range(1, t):
            P = PseudoProfile(motifs)
            motifs.append(ProfileMostPr_kmer(dna[i], k, P))
        currScore = Score(motifs)
        if currScore < BestScore:
            BestMotifs = motifs
            BestScore = currScore
    return BestMotifs 

if __name__ == "__main__":
    data = list(sys.stdin.read().strip().split())
    k = int(data[0])
    t = int(data[1])
    dna = data[2:]
    BestMotifs = GreedyMotifSearch2(dna, k, t)
    for motif in BestMotifs:
        print(motif)