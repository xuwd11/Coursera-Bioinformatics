# python3
def PatternToNumber(pattern):
    n = 0
    seq1 = 'ACGTacgt01230123'
    seq_dict = { seq1[i]:int(seq1[i+8]) for i in range(8) }
    l = len(pattern.strip())
    for i in range(l):
        n += seq_dict[pattern[i]]*4**(l-i-1)
    return n 

def NumberToPattern(n, k):
    p = []
    seq1 = '0123ACGT'
    seq_dict = { int(seq1[i]):seq1[i+4] for i in range(4) }
    for i in range(k):
        p.insert(0, seq_dict[n % 4])
        n //= 4
    return ''.join(p)
    
def ComputingFrequencies(text, k):
    fa = [0] * 4**k
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        j = PatternToNumber(pattern)
        fa[j] += 1
    return fa
          
    
if __name__ == "__main__":
    #print(PatternToNumber(input()))
    #print(NumberToPattern(int(input()), int(input())))    
    text = input().strip()
    k = int(input().strip())
    fa = ComputingFrequencies(text, k)
    print(' '.join([str(n) for n in fa]))