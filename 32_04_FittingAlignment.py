# python3

import sys
import numpy as np

'''
Code Challenge: Solve the Fitting Alignment Problem.
     Input: Two nucleotide strings v and w, where v has length at most 1000 and w has length at most 100.
     Output: A highest-scoring fitting alignment between v and w. Use the simple scoring method in which matches count +1 and both the
     mismatch and indel penalties are 1.

Sample Input:
     GTAGGCTTAAGGTTA
     TAGATA

Sample Output:
     2
     TAGGCTTA
     TAGA--TA
'''

class FittingAlignment:
    def __init__(self):
        self._input()
        maxScore, s1, s2 = self.computeCLS(self.seq1, self.seq2)
        print(maxScore)
        print(s1)
        print(s2)

    def _input(self):
        self.seq1, self.seq2 = sys.stdin.read().strip().split('\n')
        self.seq1 = self.seq1.strip()
        self.seq2 = self.seq2.strip()
    
    def LCSBackTrack(self, v, w, penalty = 1):
        n = len(v)
        m = len(w)
        s = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        backtrack = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        '''
        for i in range(1, n+1):
            s[i, 0] = s[i-1, 0] - penalty
            backtrack[i, 0] = 1
        '''
        for j in range(1, m+1):
            s[0, j] = s[0, j-1] - penalty
            backtrack[0, j] = 2
        for i in range(1, n+1):
            for j in range(1, m+1):
                score1 = s[i-1, j] - penalty
                score2 = s[i, j-1] - penalty
                score3 = s[i-1, j-1] + 1 if v[i-1]==w[j-1] else s[i-1, j-1] - penalty
                s[i, j] = max([score1, score2, score3])
                if s[i, j] == score1:
                    backtrack[i, j] = 1 # down
                elif s[i, j] == score2:
                    backtrack[i, j] = 2 # right
                elif s[i, j] == score3:
                    backtrack[i, j] = 3
        i = np.argmax(s[:, m])
        j = m
        return s[i, j], backtrack, i, j

    def OutputLCS(self, backtrack, v, w, i, j):
        s1 = ''
        s2 = ''
        while i > 0 or j > 0:
            if 1 == backtrack[i, j]:
                i -= 1
                s1 = v[i] + s1
                s2 = '-' + s2
                continue
            elif 2 == backtrack[i, j]:
                j -= 1
                s1 = '-' + s1
                s2 = w[j] + s2
                continue
            elif 3 == backtrack[i, j]:
                i -= 1
                j -= 1
                s1 = v[i] + s1
                s2 = w[j] + s2
                continue
            else:
                i = 0
                j = 0
        return s1, s2
    
    def computeCLS(self, seq1, seq2):
        maxSocre, backtrack, i, j = self.LCSBackTrack(seq1, seq2)
        s1, s2 = self.OutputLCS(backtrack, seq1, seq2, i, j)
        return maxSocre, s1, s2

if __name__ == "__main__":
    FittingAlignment()