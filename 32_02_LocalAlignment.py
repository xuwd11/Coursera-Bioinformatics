# python3

import sys
import numpy as np

'''
Solve the Local Alignment Problem.
     Input: Two protein strings written in the single-letter amino acid alphabet.
     Output: The maximum score of a local alignment of the strings, followed by a local alignment of these strings achieving the maximum
     score. Use the PAM250 scoring matrix and indel penalty σ = 5.

Sample Input:
     MEANLY
     PENALTY

Sample Output:
     15
     EANL-Y
     ENALTY
'''

class LocalLCS:
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
    
    def scoringMatrix(self): # PAM250 scoring matrix
        sMatrixTxt = '''
   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  2 -2  0  0 -3  1 -1 -1 -1 -2 -1  0  1  0 -2  1  1  0 -6 -3
C -2 12 -5 -5 -4 -3 -3 -2 -5 -6 -5 -4 -3 -5 -4  0 -2 -2 -8  0
D  0 -5  4  3 -6  1  1 -2  0 -4 -3  2 -1  2 -1  0  0 -2 -7 -4
E  0 -5  3  4 -5  0  1 -2  0 -3 -2  1 -1  2 -1  0  0 -2 -7 -4
F -3 -4 -6 -5  9 -5 -2  1 -5  2  0 -3 -5 -5 -4 -3 -3 -1  0  7
G  1 -3  1  0 -5  5 -2 -3 -2 -4 -3  0  0 -1 -3  1  0 -1 -7 -5
H -1 -3  1  1 -2 -2  6 -2  0 -2 -2  2  0  3  2 -1 -1 -2 -3  0
I -1 -2 -2 -2  1 -3 -2  5 -2  2  2 -2 -2 -2 -2 -1  0  4 -5 -1
K -1 -5  0  0 -5 -2  0 -2  5 -3  0  1 -1  1  3  0  0 -2 -3 -4
L -2 -6 -4 -3  2 -4 -2  2 -3  6  4 -3 -3 -2 -3 -3 -2  2 -2 -1
M -1 -5 -3 -2  0 -3 -2  2  0  4  6 -2 -2 -1  0 -2 -1  2 -4 -2
N  0 -4  2  1 -3  0  2 -2  1 -3 -2  2  0  1  0  1  0 -2 -4 -2
P  1 -3 -1 -1 -5  0  0 -2 -1 -3 -2  0  6  0  0  1  0 -1 -6 -5
Q  0 -5  2  2 -5 -1  3 -2  1 -2 -1  1  0  4  1 -1 -1 -2 -5 -4
R -2 -4 -1 -1 -4 -3  2 -2  3 -3  0  0  0  1  6  0 -1 -2  2 -4
S  1  0  0  0 -3  1 -1 -1  0 -3 -2  1  1 -1  0  2  1 -1 -2 -3
T  1 -2  0  0 -3  0 -1  0  0 -2 -1  0  0 -1 -1  1  3  0 -5 -3
V  0 -2 -2 -2 -1 -1 -2  4 -2  2  2 -2 -1 -2 -2 -1  0  4 -6 -2
W -6 -8 -7 -7  0 -7 -3 -5 -3 -2 -4 -4 -6 -5  2 -2 -5 -6 17  0
Y -3  0 -4 -4  7 -5  0 -1 -4 -1 -2 -2 -5 -4 -4 -3 -3 -2  0 10
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

    def LCSBackTrack(self, v, w, sigma = 5):
        sMatrix = self.scoringMatrix()
        n = len(v)
        m = len(w)
        s = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        backtrack = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        '''
        for i in range(1, n+1):
            s[i, 0] = s[i-1, 0] - sigma
            backtrack[i, 0] = 1
        for j in range(1, m+1):
            s[0, j] = s[0, j-1] - sigma
            backtrack[0, j] = 2
        '''
        for i in range(1, n+1):
            for j in range(1, m+1):
                score1 = s[i-1, j] - sigma
                score2 = s[i, j-1] - sigma
                score3 = s[i-1, j-1] + sMatrix[v[i-1]][w[j-1]]
                score4 = 0
                s[i, j] = max([score4, score1, score2, score3])
                if s[i, j] == score4:
                    backtrack[i, j] = 0
                elif s[i, j] == score1:
                    backtrack[i, j] = 1 # down
                elif s[i, j] == score2:
                    backtrack[i, j] = 2 # right
                elif s[i, j] == score3:
                    backtrack[i, j] = 3
        ind = np.argmax(s)
        j = ind % (m+1)
        i = (ind-m-1) // (m+1) + 1
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
    LocalLCS()