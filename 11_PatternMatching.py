# python3
from sys import stdin

def PatternMatching(text, pattern):
    pos = []
    for i in range(len(text)-len(pattern)+1):
        #print(i)
        if text[i:(i+len(pattern))] == pattern:
            pos.append(i)
    return pos
    
    
if __name__ == "__main__":
    pattern = input().strip()
    text = input().strip()
    pos = PatternMatching(text, pattern)
    print(' '.join([str(i) for i in pos]))
    '''
    pattern = input()
    f = open('week1_Vibrio_cholerae.txt', 'r')
    for line in f:
        text = line
    f = open('dataset_3_5.txt', 'r')
    t = []
    for line in f:
        t.append(line.strip())
    pattern = t[0]
    text = t[1]
    pos = PatternMatching(text, pattern)
    print(' '.join([str(i) for i in pos]))
    '''

    