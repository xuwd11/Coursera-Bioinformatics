# put your python code here
# python3

import sys

def PairedComposition(k, d, text):
    l = len(text)
    c = []
    for i in range(l-2*k-d+1):
        c.append((text[i:i+k], text[i+k+d:i+2*k+d]))
    #c.sort()
    return c

class StringReconstruction2: #String reconstruction from paired reads
    def __init__(self):
        self.readData()
        #print(self.patterns)
        self.string = self.StringSpelledByGappedPatterns(self.patterns, self.k, self.d)
        print(self.string)

    def readData(self):
        data = list(sys.stdin.read().strip().split())
        self.k, self.d = int(data[0]), int(data[1])
        self.patterns = [tuple(p.split('|')) for p in data[2:]]
    
    def StringSpelledByGappedPatterns(self, patterns, k, d):
        firstPatterns = patterns[0][0] + ''.join([p[0][-1] for p in patterns[1:]])
        secondPatterns = patterns[0][1] + ''.join([p[1][-1] for p in patterns[1:]])
        l = len(firstPatterns)
        if firstPatterns[k+d:] == secondPatterns[:l-k-d]:
            return firstPatterns + secondPatterns[-(k+d):]
        else:
            return 'there is no string spelled by the gapped patterns'        

if __name__ == "__main__":
    StringReconstruction2()