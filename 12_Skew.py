# python3

def SkewDiagram(seq):
    s = [0]
    skew = 0
    seq1 = 'ATGatg001001'
    seq_dict = { seq1[i]:int(seq1[i+6]) for i in range(6) }
    seq_dict['C'] = -1
    seq_dict['c'] = -1
    for nucleotide in seq:
        skew += seq_dict[nucleotide]
        s.append(skew)
    return s
    
def minSkew(seq):
    s = SkewDiagram(seq)
    m = min(s)
    ans = [i for i, v in enumerate(s) if v == m]
    return ans
   
if __name__ == "__main__":
    print(' '.join([str(i) for i in minSkew(input().strip())]))
    '''
    f = open('dataset_7_6.txt', 'r')
    for line in f:
        seq = line.strip()
    mSkew = minSkew(seq)
    print(' '.join([str(i) for i in mSkew])) 
    '''   