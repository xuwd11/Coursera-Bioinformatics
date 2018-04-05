# python3

import sys

def ReverseComplement(seq):
    for base in seq:
        if base not in 'ATCGatcg':
            print("Error: NOT a DNA sequence")
            return None
    seq1 = 'ATCGTAGCatcgtagc'
    seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
    return "".join([seq_dict[base] for base in reversed(seq)])
    
def HammingDistance(seq1, seq2):
    return len([i for i in range(len(seq1)) if seq1[i] != seq2[i]])
    
def ApproxPatternMatching(pattern, text, d):
    pos = []
    l = len(pattern)
    for i in range(len(text)-l+1):
        if HammingDistance(pattern, text[i:i+l]) <= d:
            pos.append(i)
    return pos

def ApproxPatternCount(pattern, text, d):
    c = 0
    l = len(pattern)  
    for i in range(len(text)-l+1):
        if HammingDistance(pattern, text[i:i+l]) <= d:
            c += 1
    return c

def ImmediateNeighbors(pattern):
    neighbor = set()
    nset = {'A', 'C', 'G', 'T'}
    for i in range(len(pattern)):
        for n in nset:
            neighbor.add(pattern[:i]+n+pattern[i+1:])
    return neighbor

def Neighbors(pattern, d):
    if d == 0:
        return {pattern}
    ineighbor = ImmediateNeighbors(pattern)
    neighbor = ineighbor
    for j in range(d-1):
        for p in ineighbor:
            neighbor = neighbor.union(ImmediateNeighbors(p))
        ineighbor = neighbor
    return neighbor

def FrequentWordsWithMismatches(text, k, d):
    counts = dict()
    for i in range(len(text)-k+1):
        neighbor = Neighbors(text[i:i+k], d)
        for n in neighbor:
            counts[n] = counts.get(n, 0) + 1
    m = max(counts.values())
    return [t for t, v in counts.items() if v == m]   

def FrequentWordsWMARC(text, k, d):
    #Frequent Words with Mismatches and Reverse Complements
    counts = dict()
    for i in range(len(text)-k+1):
        neighbor = Neighbors(text[i:i+k], d)
        for n in neighbor:
            nrc = ReverseComplement(n)
            counts[n] = counts.get(n, 0) + 1
            counts[nrc] = counts.get(nrc, 0) + 1
    m = max(counts.values())
    return [t for t, v in counts.items() if v == m] 

def DoesPatternAppear(pattern, text, d):
    k = len(pattern)
    l = len(text)
    for i in range(l-k+1):
        if HammingDistance(pattern, text[i:i+k]) <= d:
            return True
    return False               

def MotifEnumeration(Dna, k, d):
    patterns = set()
    neighbor = set()
    text = Dna[0]
    for i in range(len(text)-k+1):
        n = Neighbors(text[i:i+k], d)
        neighbor = neighbor.union(n)
    for n in neighbor:
        ValidPattern = True
        for i in range(1, len(Dna)):
            if not DoesPatternAppear(n, Dna[i], d):
                ValidPattern = False
                break
        if ValidPattern:
            patterns.add(n)
    return patterns        

if __name__ == "__main__":
    '''
    n = int(input().strip())
    k, d = map(int, input().strip().split())
    Dna = []
    for i in range(n):
        Dna.append(input().strip())
    '''
    i = sys.stdin.read()
    k, d, *Dna = list(i.split())
    k = int(k)
    d = int(d)
    p = MotifEnumeration(Dna, k, d)
    print(' '.join(p))