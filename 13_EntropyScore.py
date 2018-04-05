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

def EntropyScore(P):
    k = len(P[0])
    Pt = [[row[i] for row in P] for i in range(k)]
    entropy = [0 for _ in range(k)]
    for i in range(k):
        for s in Pt[i]:
            entropy[i] += 0 if s == 0 else -s * math.log(s, 2)
    return entropy

if __name__ == "__main__":
    data = sys.stdin.read()
    motifs = list(data.split())
    P = Profile(motifs)
    entropy = EntropyScore(P)
    e = sum(entropy)
    print(e)