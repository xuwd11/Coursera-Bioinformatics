# python3

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
   
if __name__ == "__main__":
    #print(HammingDistance(input().strip(), input().strip()))
    '''
    f = open('dataset_9_4.txt', 'r')
    seq = []
    for line in f:
        seq.append(line.strip())
    pos = ApproxPatternMatching(seq[0], seq[1], int(seq[2]))
    print(' '.join([str(p) for p in pos]))
    
    pos = ApproxPatternMatching(input().strip(), input().strip(), int(input().strip()))
    print(' '.join([str(p) for p in pos])) 
    '''
    print(ApproxPatternCount(input().strip(), input().strip(), int(input().strip())))  