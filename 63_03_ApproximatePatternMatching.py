# python3
import sys

'''
Solve the Multiple Approximate Pattern Matching Problem.
     Input: A string Text, followed by a collection of strings Patterns, and an integer d.
     Output: All positions where one of the strings in Patterns appears as a substring of Text with at most d mismatches.

Sample Input:
ACATGCTACTTT
ATT GCC GCTA TATT
1
Sample Output:
2 4 4 6 7 8 9

Theorem: If two strings of length n match with at most d mismatches, then they must share a k-mer of length 
k=⌊n/(d+1)⌋.

We now have the outline of an algorithm for matching a string Pattern of length n to Text with at most d mismatches. We first 
divide Pattern into d+1 segments of length k=⌊n/(d+1)⌋, called seeds. After finding which seeds match Text exactly (seed 
detection), we attempt to extend seeds in both directions in order to verify whether Pattern occurs in Text with at most d 
mismatches (seed extension).
'''

class SuffixArray:
    '''
    Build suffix array of the string text and
    return a list result of the same length as the text
    such that the value result[i] is the index (0-based)
    in text where the i-th lexicographically smallest
    suffix of text starts.
    '''
    def __init__(self, text):
        self.order = self.buildSuffixArray(text)
    
    def _input(self):
        self.text = sys.stdin.readline().strip()

    def sortCharacters(self, S):
        l = len(S)
        order = [0] * l
        count = dict()
        for i in range(l):
            count[S[i]] = count.get(S[i], 0) + 1
        charList = sorted(count.keys())
        prevChar = charList[0]
        for char in charList[1:]:
            count[char] += count[prevChar]
            prevChar = char
        for i in range(l-1, -1, -1):
            c = S[i]
            count[c] = count[c] - 1
            order[count[c]] = i
        return order

    def computeCharClasses(self, S, order):
        l = len(S)
        charClass = [0] * l
        charClass[order[0]] = 0
        for i in range(1, l):
            if S[order[i]] != S[order[i-1]]:
                charClass[order[i]] = charClass[order[i-1]] + 1
            else:
                charClass[order[i]] = charClass[order[i-1]]
        return charClass        
    
    def sortDoubled(self, S, L, order, _class):
        sLen = len(S)
        count = [0] * sLen
        newOrder = [0] * sLen
        for i in range(sLen):
            count[_class[i]] += 1
        for j in range(1, sLen):
            count[j] += count[j-1]
        for i in range(sLen-1, -1, -1):
            start = (order[i]-L+sLen) % sLen
            cl = _class[start]
            count[cl] -= 1
            newOrder[count[cl]] = start
        return newOrder
    
    def updateClasses(self, newOrder, _class, L):
        n = len(newOrder)
        newClass = [0] * n
        newClass[newOrder[0]] = 0
        for i in range(1, n):
            curr = newOrder[i]
            prev = newOrder[i-1]
            mid = curr + L
            midPrev = (prev + L) % n
            if _class[curr] != _class[prev] or _class[mid] != _class[midPrev]:
                newClass[curr] = newClass[prev] + 1
            else:
                newClass[curr] = newClass[prev]
        return newClass
    
    def buildSuffixArray(self, S):
        sLen = len(S)
        order = self.sortCharacters(S)
        _class = self.computeCharClasses(S, order)
        L = 1
        while L < sLen:
            order = self.sortDoubled(S, L, order, _class)
            _class = self.updateClasses(order, _class, L)
            L = 2 * L
        return order

class ApproxPatternMatching:
    def __init__(self):
        #text, patterns, d = self._input()
        text, patterns, d = self.readFromFile()
        occs = self.findOccurrences(text, patterns, d)
        print(' '.join(map(str, sorted(occs))))
        f = open('result.txt', 'w')
        f.write(' '.join(map(str, sorted(occs))))
        f.close()
    
    def _input(self):
        text = sys.stdin.readline().strip()
        text = text + '$'
        data = sys.stdin.read().strip().split()
        patterns = data[:-1]
        d = int(data[-1])
        return text, patterns, d

    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().strip().split()
        text = data[0] + '$'
        patterns = data[1:-1]
        d = int(data[-1])
        return text, patterns, d

    def bwtFromSuffixArray(self, text, order, alphabet = ['$', 'A', 'C', 'G', 'T']):
        l = len(text)
        bwt = [''] * l
        for i in range(l):
            bwt[i] = text[(order[i]+l-1)%l]

        counts = dict()
        starts = dict()
        for char in alphabet:
            counts[char] = [0] * (l + 1)
        for i in range(l):
            currChar = bwt[i]
            for char, count in counts.items():
                counts[char][i+1] = counts[char][i]
            counts[currChar][i+1] += 1
        currIndex = 0
        for char in sorted(alphabet):
            starts[char] = currIndex
            currIndex += counts[char][l]
        return bwt, starts, counts
        
    def approxMatching(self, p1, p2, d):
        error = 0
        for i in range(len(p1)):
            if p1[i] != p2[i]:
                error += 1
                if error > d:
                    return False
        return True

    def findOccurrences(self, text, patterns, d):
        order = SuffixArray(text).order
        bwt, starts, counts = self.bwtFromSuffixArray(text, order)        
        
        occs = []
        for pattern in patterns:
            currOccs = set()
            n = len(pattern)
            k = n // (d+1)
            seeds = [(pattern[i*k:(i+1)*k], i*k) for i in range(d)] + [(pattern[d*k:n], d*k)]
            
            for p, offset in seeds:
                top = 0
                bottom = len(bwt) - 1
                currIndex = len(p) - 1
                while top <= bottom:
                    if currIndex >= 0:
                        symbol = p[currIndex]
                        currIndex -= 1
                        if counts[symbol][bottom+1] - counts[symbol][top] > 0:
                            top = starts[symbol] + counts[symbol][top]
                            bottom = starts[symbol] + counts[symbol][bottom+1] - 1
                        else:
                            break
                    else:
                        for i in range(top, bottom + 1):
                            currOccs.add(order[i]-offset)
                        break
            for occ in currOccs:
                if self.approxMatching(text[occ:occ+n], pattern, d):
                    occs.append(occ)
        return occs
    
if __name__ == '__main__':
    ApproxPatternMatching()