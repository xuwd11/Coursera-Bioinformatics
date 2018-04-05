# python3

import sys
import numpy as np

'''
Solve the Alignment with Affine Gap Penalties Problem.
     Input: Two amino acid strings v and w (each of length at most 100).
     Output: The maximum alignment score between v and w, followed by an alignment of v and w achieving this maximum score. Use the
     BLOSUM62 scoring matrix, a gap opening penalty of 11, and a gap extension penalty of 1.

Sample Input:
     PRTEINS
     PRTWPSEIN

Sample Output:
     8
     PRT---EINS
     PRTWPSEIN-
'''

class GlobalLCSGapPenalties:
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
    
    def scoringMatrix(self): # BLOSUM62 scoring matrix
        sMatrixTxt = '''
   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7
'''
        sMatrixList = sMatrixTxt.strip().split('\n')
        aaList = sMatrixList[0].split()
        sMatrix = dict()
        for aa in aaList:
            sMatrix[aa] = dict()
        for i in range(1, len(aaList) + 1):
            currRow = sMatrixList[i].split()
            for j in range(len(aaList)):
                sMatrix[currRow[0]][aaList[j]] = int(currRow[j + 1])
        return sMatrix

    def LCSBackTrack(self, v, w, sigma = 11, epsilon = 1, inf = 1000000):
        sMatrix = self.scoringMatrix()
        n = len(v)
        m = len(w)
        sLower = np.matrix(-inf*np.ones((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        sMiddle = np.matrix(-inf*np.ones((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        sUpper = np.matrix(-inf*np.ones((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        s = [sLower, sMiddle, sUpper]
        backtrackLower = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        backtrackMiddle = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        backtrackUpper = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        backtrack = [backtrackLower, backtrackMiddle, backtrackUpper]
        s[0][0, 0] = 0
        s[1][0, 0] = 0
        s[2][0, 0] = 0
        s[0][1, 0] = -sigma
        s[1][1, 0] = -sigma
        for i in range(2, n+1):
            s[0][i, 0] = s[0][i-1, 0] - epsilon
            s[1][i, 0] = s[0][i, 0]
        s[2][0, 1] = -sigma
        s[1][0, 1] = -sigma
        for j in range(2, m+1):
            s[2][0, j] = s[2][0, j-1] - epsilon
            s[1][0, j] = s[2][0, j]
        for i in range(1, n+1):
            for j in range(1, m+1):
                score1 = s[0][i-1, j] - epsilon
                score2 = s[1][i-1, j] - sigma
                lowerScore = max(score1, score2)
                s[0][i, j] = lowerScore
                if lowerScore == score2:
                    backtrack[0][i, j] = 1 # From middle
                score1 = s[2][i, j-1] - epsilon
                score2 = s[1][i, j-1] - sigma
                upperScore = max(score1, score2)
                s[2][i, j] = upperScore
                if upperScore == score2:
                    backtrack[2][i, j] = 1 # From middle
                score3 = s[1][i-1, j-1] + sMatrix[v[i-1]][w[j-1]]
                middleScore = max(lowerScore, score3, upperScore)
                s[1][i, j] = middleScore
                if middleScore == score3:
                    backtrack[1][i, j] = 1 # From middle
                elif middleScore == upperScore:
                    backtrack[1][i, j] = 2 # From upper
        return s[1][n, m], backtrack, n, m                       


    def OutputLCS(self, backtrack, v, w, i, j):
        s1 = ''
        s2 = ''
        level = 1 # middle
        while i > 0 or j > 0:
            if 0 == level:
                if 1 == backtrack[0][i, j]:
                    level = 1
                i -= 1
                s1 = v[i] + s1
                s2 = '-' + s2
                continue
            elif 2 == level:
                if 1 == backtrack[2][i, j]:
                    level = 1
                j -= 1
                s1 = '-' + s1
                s2 = w[j] + s2
                continue
            else:
                if 1 == backtrack[1][i, j]:
                    i -= 1
                    j -= 1
                    s1 = v[i] + s1
                    s2 = w[j] + s2
                    continue
                else:
                    level = backtrack[1][i, j]
        return s1, s2
    
    def computeCLS(self, seq1, seq2):
        maxSocre, backtrack, i, j = self.LCSBackTrack(seq1, seq2)
        s1, s2 = self.OutputLCS(backtrack, seq1, seq2, i, j)
        return maxSocre, s1, s2

if __name__ == "__main__":
    GlobalLCSGapPenalties()