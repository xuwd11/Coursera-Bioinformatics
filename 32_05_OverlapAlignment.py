# python3

import sys
import numpy as np

'''
An overlap alignment of strings v = v1 ... vn and w = w1 ... wm is a global alignment of a suffix of v with a prefix of w. An optimal 
overlap alignment of strings v and w maximizes the global alignment score between an i-suffix of v and a j-prefix of w (i.e., between 
vi ... vn and w1 ... wj) among all i and j.

Overlap Alignment Problem: Construct a highest-scoring overlap alignment between two strings.
     Input: Two strings and a matrix score.
     Output: A highest-scoring overlap alignment between the two strings as defined by the scoring matrix score.

Sample Input:
     PAWHEAE
     HEAGAWGHEE

Sample Output:
     1
     HEAE
     HEAG
'''

class OverlapAlignment:
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
    
    def LCSBackTrack(self, v, w, penalty = 2):
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
        i = n
        j = np.argmax(s[n, :])
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
    OverlapAlignment()