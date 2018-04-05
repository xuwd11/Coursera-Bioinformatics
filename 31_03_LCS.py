# python3

import sys
import numpy as np

'''
Use OutputLCS (reproduced below) to solve the Longest Common Subsequence Problem.
     Input: Two strings s and t.
     Output: A longest common subsequence of s and t. (Note: more than one solution may exist, in which case you may output any one.)

Sample Input:
     AACCTTGG
     ACACTGTGA

Sample Output:
     AACTGG

Pseudocode:
    LCSBackTrack(v, w)
        for i ← 0 to |v|
            si, 0 ← 0
        for j ← 0 to |w| 
            s0, j ← 0
        for i ← 1 to |v|
            for j ← 1 to |w|
                si, j ← max{si-1, j , si,j-1 , si-1, j-1 + 1 (if vi = wj)}
                if si,j = si-1,j
                    Backtracki, j ← "↓"
                else if si, j = si, j-1
                    Backtracki, j ← "→"
                else if si, j = si-1, j-1 + 1 and vi = wj
                    Backtracki, j ← "↘"
        return Backtrack

    OutputLCS(backtrack, v, i, j)
        if i = 0 or j = 0
            return
        if backtracki, j = "↓"
            OutputLCS(backtrack, v, i - 1, j)
        else if backtracki, j = "→"
            OutputLCS(backtrack, v, i, j - 1)
        else
            OutputLCS(backtrack, v, i - 1, j - 1)
            output vi
'''

class LCS:
    def __init__(self):
        self._input()
        self.backtrack = self.LCSBackTrack(self.seq1, self.seq2)
        self.lcs = ''
        self.OutputLCS(self.backtrack, self.seq1, self.n, self.m)
        print(self.lcs)
    
    def _input(self):
        self.seq1, self.seq2 = sys.stdin.read().strip().split('\n')
        self.n = len(self.seq1)
        self.m = len(self.seq2)
        
    def LCSBackTrack(self, v, w):
        n = len(v)
        m = len(w)
        s = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        backtrack = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        for i in range(1, n+1):
            for j in range(1, m+1):
                s[i, j] = max([s[i-1, j], s[i, j-1], s[i-1, j-1] + 1 if v[i-1]==w[j-1] else s[i-1, j-1]])
                if s[i, j] == s[i-1, j]:
                    backtrack[i, j] = 1 #down
                elif s[i, j] == s[i, j-1]:
                    backtrack[i, j] = 2 #right
                elif s[i, j] == s[i-1, j-1] + 1 and v[i-1] == w[j-1]:
                    backtrack[i, j] = 3
        return backtrack

    def OutputLCS(self, backtrack, v, i, j):
        while i > 0 and j > 0:
            if 1 == backtrack[i, j]:
                i -= 1
                continue
            elif 2 == backtrack[i, j]:
                j -= 1
                continue
            else:
                i -= 1
                j -= 1
                self.lcs = v[i] + self.lcs

    def OutputLCS0(self, backtrack, v, i, j):
        if 0 == i or 0 == j:
            return
        if 1 == backtrack[i, j]:
            self.OutputLCS(backtrack, v, i-1, j)
        elif 2 == backtrack[i, j]:
            self.OutputLCS(backtrack, v, i, j-1)
        else:
            self.OutputLCS(backtrack, v, i-1, j-1)
            self.lcs += v[i-1]

if __name__ == "__main__":
    LCS()