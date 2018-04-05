# python3

import sys
import numpy as np

'''
Solve the Global Alignment Problem.
     Input: Two protein strings written in the single-letter amino acid alphabet.
     Output: The maximum alignment score of these strings followed by an alignment achieving this maximum score. Use the
     BLOSUM62 scoring matrix and indel penalty σ = 5.

Sample Input:
     PLEASANTLY
     MEANLY

Sample Output:
     8
     PLEASANTLY
     -MEA--N-LY
'''

class GlobalLCS:
    def __init__(self):
        #self._input()
        self.readFile()
        maxScore, s1, s2 = self.computeCLS(self.seq1, self.seq2)
        self.saveData([str(maxScore), s1, s2])
        print(maxScore)
        print(s1)
        print(s2)

    def _input(self):
        self.seq1, self.seq2 = sys.stdin.read().strip().split('\n')
        self.seq1 = self.seq1.strip()
        self.seq2 = self.seq2.strip()
    
    def readFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        self.seq1 = data[0]
        self.seq2 = data[1]
        self.n = len(self.seq1)
        self.m = len(self.seq2)

    def saveData(self, data):
        f = open('result.txt', 'w')
        for d in data:
            f.write(d+'\n')
    
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

    def LCSBackTrack(self, v, w, sigma = 5):
        sMatrix = self.scoringMatrix()
        n = len(v)
        m = len(w)
        s = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        backtrack = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        for i in range(1, n+1):
            s[i, 0] = s[i-1, 0] - sigma
            backtrack[i, 0] = 1
        for j in range(1, m+1):
            s[0, j] = s[0, j-1] - sigma
            backtrack[0, j] = 2
        for i in range(1, n+1):
            for j in range(1, m+1):
                score1 = s[i-1, j] - sigma
                score2 = s[i, j-1] - sigma
                score3 = s[i-1, j-1] + sMatrix[v[i-1]][w[j-1]]
                s[i, j] = max([score1, score2, score3])
                if s[i, j] == score1:
                    backtrack[i, j] = 1 # down
                if s[i, j] == score2:
                    backtrack[i, j] = 2 # right
                if s[i, j] == score3:
                    backtrack[i, j] = 3
        return s[n, m], backtrack

    def OutputLCS(self, backtrack, v, w):
        s1 = ''
        s2 = ''
        i = len(v)
        j = len(w)
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
            else:
                i -= 1
                j -= 1
                s1 = v[i] + s1
                s2 = w[j] + s2
        return s1, s2
    
    def computeCLS(self, seq1, seq2):
        maxSocre, backtrack = self.LCSBackTrack(seq1, seq2)
        s1, s2 = self.OutputLCS(backtrack, seq1, seq2)
        return maxSocre, s1, s2

if __name__ == "__main__":
    GlobalLCS()