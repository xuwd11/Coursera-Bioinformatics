# python3

def ReverseComplement(seq):
    for base in seq:
        if base not in 'ATCGatcg':
            print("Error: NOT a DNA sequence")
            return None
    seq1 = 'ATCGTAGCatcgtagc'
    seq_dict = {seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12}
    return "".join([seq_dict[base] for base in reversed(seq)])

def RNAReverseComplement(seq):
    seq1 = 'AUCGaucgUAGCUAGC'
    seq_dict = {seq1[i]:seq1[i+8] for i in range(8)}
    return ''.join([seq_dict[base] for base in reversed(seq)])

def transcription(seq):
    seq1 = 'ATCGatcgAUCGAUCG'
    seq_dict = {seq1[i]:seq1[i+8] for i in range(8)}
    return ''.join([seq_dict[base] for base in seq])

def translation(seq):
    codonTable = '''
AAA K
AAC N
AAG K
AAU N
ACA T
ACC T
ACG T
ACU T
AGA R
AGC S
AGG R
AGU S
AUA I
AUC I
AUG M
AUU I
CAA Q
CAC H
CAG Q
CAU H
CCA P
CCC P
CCG P
CCU P
CGA R
CGC R
CGG R
CGU R
CUA L
CUC L
CUG L
CUU L
GAA E
GAC D
GAG E
GAU D
GCA A
GCC A
GCG A
GCU A
GGA G
GGC G
GGG G
GGU G
GUA V
GUC V
GUG V
GUU V
UAA 0
UAC Y
UAG 0
UAU Y
UCA S
UCC S
UCG S
UCU S
UGA 0
UGC C
UGG W
UGU C
UUA L
UUC F
UUG L
UUU F'''
    codon = codonTable.split()
    codonDict = {codon[i*2]:codon[i*2+1] for i in range(int(len(codon)/2))}
    protein1 = ''.join([codonDict[seq[3*i:3*i+3]] for i in range(int(len(seq)/3))])
    return protein1.split('0')[0]

def PeptideEncoding(text, pattern):
    l = len(text)
    k = len(pattern)
    if l < k*3:
        return
    rna = transcription(text)
    subStrings = [text[i:i+k*3] for i in range(l-k*3) if translation(rna[i:i+k*3]) == pattern or translation(RNAReverseComplement(rna[i:i+k*3])) == pattern]
    return subStrings
    
if __name__ == "__main__":
    '''
    subStrings = PeptideEncoding(input().strip(), input().strip())
    for s in subStrings:
        print(s)
    #Tyrocidine B1 (Val-Lys-Leu-Phe-Pro-Trp-Phe-Asn-Gln-Tyr)
    #VKLFPWFNQY
    
    f = open('Bacillus_brevis.txt', 'r')
    data = ''
    for line in f:
        data = data + line.strip()
    subStrings = PeptideEncoding(data, 'VKLFPWFNQY')
    print(data)
    print(len(subStrings))
    '''    
    print(translation(input().strip()))