# python3

import sys
import math

def HammingDistance(seq1, seq2):
    return len([i for i in range(len(seq1)) if seq1[i] != seq2[i]])

def PatternStringDistance(pattern, dna):
    k = len(pattern)
    distance = 0
    for seq in dna:
        l = len(seq)
        hd = float('inf')
        for i in range(l - k +1):
            hdCurr = HammingDistance(pattern, seq[i:i+k])
            if hd > hdCurr:
                hd = hdCurr
        distance += hd
    return distance

def NumberToPattern(n, k):
    p = []
    seq1 = '0123ACGT'
    seq_dict = { int(seq1[i]):seq1[i+4] for i in range(4) }
    for i in range(k):
        p.insert(0, seq_dict[n % 4])
        n //= 4
    return ''.join(p)

def MedianString(dna, k):
    distance = float('inf')
    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        dCurr = PatternStringDistance(pattern, dna)
        if distance > dCurr:
            distance = dCurr
            median = pattern
    return median

def AllMedianString(dna, k):
    distance = float('inf')
    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        dCurr = PatternStringDistance(pattern, dna)
        if distance > dCurr:
            distance = dCurr
            median = []
            median.append(pattern)
        elif distance == dCurr:
            median.append(pattern)
    return median

if __name__ == "__main__":
    data = sys.stdin.read()
    data = list(data.strip().split())
    k = int(data[0].strip())
    dna = data[1:]
    pattern = AllMedianString(dna, k)
    print(pattern)