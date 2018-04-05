# python3
def PatternCount(text, pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if text[i:(i+len(pattern))] == pattern:
            count += 1
    return count

def FrequentWords(text, k):
    counts = dict()
    for i in range(len(text)-k+1):
        counts[text[i:(i+k)]] = counts.get(text[i:(i+k)], 0) + 1
    frequency = max(counts.values())
    return [t for t, v in counts.items() if v == frequency]

def ClumpFinding(genome, k, t, L):
    counts = dict()
    for i in range(L):
        text = genome[i:(i+k)]
        counts[text] = counts.get(text, 0) + 1
    clump = [p for p, v in counts.items() if v >= t]
    for i in range(1, len(genome)-L+1):
        firstPattern = genome[(i-1):(i-1+k)]
        counts[firstPattern] -= 1
        lastPattern = genome[(i+L-k):(i+L)]
        counts[lastPattern] = counts.get(lastPattern, 0) + 1
        if counts[lastPattern] >= t and lastPattern not in clump:
            clump.append(lastPattern)
    return clump
      
if __name__ == "__main__":
    '''
    text = input().strip()
    k, L, t = map(int, input().strip().split())
    clump = ClumpFinding(text, k, t, L)
    print(' '.join(clump))
    
    f = open('dataset_4_5.txt', 'r')
    c = []
    for line in f:
        c.append(line)
    k, L, t = map(int, c[1].strip().split())
    clump = ClumpFinding(c[0], k, t, L)
    print(' '.join(clump))
    '''
    f = open('E-coli.txt', 'r')
    for line in f:
        genome = line
    k, L, t = (9, 500, 3)
    clump = ClumpFinding(genome, k, t, L)
    print(len(clump))