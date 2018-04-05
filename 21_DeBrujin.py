# python3

import sys

def Composition(k, seq):
    l = len(seq)
    c = []
    for i in range(l-k+1):
        c.append(seq[i:i+k])
    return c

def ReconstructFromPath(patterns):
    n = len(patterns)
    seq = patterns[0]
    for i in range(1, n):
        seq += patterns[i][-1]
    return seq

def OverlapGraph(patterns):
    n = len(patterns)
    k = len(patterns[0])
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if patterns[i][1:] == patterns[j][:k-1]:
                adj[i].append(j)
    return adj

def DeBrujinFromText(k, seq):
    l = len(seq)
    adjdb = dict()
    for i in range(l-k+1):
        if seq[i:i+k-1] in adjdb:
            adjdb[seq[i:i+k-1]].append(seq[i+1:i+k])
        else:
            adjdb[seq[i:i+k-1]] = []
            adjdb[seq[i:i+k-1]].append(seq[i+1:i+k])
    return adjdb

def DeBrujin(patterns):
    k = len(patterns[0])
    adjdb = dict()
    for p in patterns:
        if p[:k-1] in adjdb:
            adjdb[p[:k-1]].append(p[1:])
        else:
            adjdb[p[:k-1]] = []
            adjdb[p[:k-1]].append(p[1:])
    return adjdb

if __name__ == "__main__":
    data = list(sys.stdin.read().strip().split())
    adjdb = DeBrujin(data)
    for seq1, seq2 in adjdb.items():
        print(seq1+' -> '+','.join(seq2))