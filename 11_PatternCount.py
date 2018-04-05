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
    
if __name__ == "__main__":
    text = input()
    ''''
    pattern = input()
    print(PatternCount(text, pattern))
    '''
    k = int(input())
    fw = FrequentWords(text, k)
    print(' '.join(fw))
   