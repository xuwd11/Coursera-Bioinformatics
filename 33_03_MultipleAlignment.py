# python3

import sys
import numpy as np

'''
Solve the Multiple Longest Common Subsequence Problem.
     Input: Three DNA strings of length at most 10.
     Output: The length of a longest common subsequence of these three strings, followed by a multiple alignment of the three strings
     corresponding to such an alignment.

Sample Input:
ATATCCG
TCCGA
ATGTACTG
Sample Output:
3
ATATCC-G-
---TCC-GA
ATGTACTG-
'''

class MultipleLCS:
    def __init__(self):
        self._input()
        #self.readFile()
        maxScore, s1, s2, s3 = self.computeCLS(self.seq1, self.seq2, self.seq3)
        self.saveData([str(maxScore), s1, s2, s3])
        print(maxScore)
        print(s1)
        print(s2)
        print(s3)

    def _input(self):
        self.seq1, self.seq2, self.seq3 = sys.stdin.read().strip().split('\n')
        self.seq1 = self.seq1.strip()
        self.seq2 = self.seq2.strip()
        self.seq3 = self.seq3.strip()
    
    def readFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        self.seq1 = data[0]
        self.seq2 = data[1]
        self.seq3 = data[2]

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

    def LCSBackTrack(self, u, v, w, score = 1):
        n = len(u)
        m = len(v)
        l = len(w)
        s = [[[0 for _ in range(l+1)] for _ in range(m+1)] for _ in range(n+1)]
        backtrack = [[[0 for _ in range(l+1)] for _ in range(m+1)] for _ in range(n+1)]
        for i in range(1, n+1):
            backtrack[i][0][0] = 1
        for j in range(1, m+1):
            backtrack[0][j][0] = 2
        for k in range(1, l+1):
            backtrack[0][0][k] = 3
        for i in range(1, n+1):
            for j in range(1, m+1):
                backtrack[i][j][0] = 4
            for k in range(1, l+1):
                backtrack[i][0][k] = 5
        for j in range(1, m+1):
            for k in range(1, l+1):
                backtrack[0][j][k] = 6
                #for i in range(1, n+1):
                    #backtrack[i][j][k] = 0
        for i in range(1, n+1):
            for j in range(1, m+1):
                for k in range(1, l+1):
                    score = [s[i-1][j-1][k-1]+1 if u[i-1]==v[j-1]==w[k-1] else s[i-1][j-1][k-1], \
                    s[i-1][j][k], s[i][j-1][k], s[i][j][k-1], s[i-1][j-1][k], s[i-1][j][k-1], s[i][j-1][k-1]]
                    index = np.argmax(score)
                    backtrack[i][j][k] = index
                    s[i][j][k] = score[index]
        return s[n][m][l], backtrack, n, m, l

    def OutputLCS(self, backtrack, u, v, w, i, j, k):
        s1 = ''
        s2 = ''
        s3 = ''
        while i > 0 or j > 0 or k > 0:
            if 1 == backtrack[i][j][k]:
                i -= 1
                s1 = u[i] + s1
                s2 = '-' + s2
                s3 = '-' + s3
                continue
            elif 2 == backtrack[i][j][k]:
                j -= 1
                s1 = '-' + s1
                s2 = v[j] + s2
                s3 = '-' + s3
                continue
            elif 3 == backtrack[i][j][k]:
                k -= 1
                s1 = '-' + s1
                s2 = '-' + s2
                s3 = w[k] + s3
                continue
            elif 4 == backtrack[i][j][k]:
                i -= 1
                j -= 1
                s1 = u[i] + s1
                s2 = v[j] + s2
                s3 = '-' + s3
                continue
            elif 5 == backtrack[i][j][k]:
                i -= 1
                k -= 1
                s1 = u[i] + s1
                s2 = '-' + s2
                s3 = w[k] + s3
                continue
            elif 6 == backtrack[i][j][k]:
                j -= 1
                k -= 1
                s1 = '-' + s1
                s2 = v[j] + s2
                s3 = w[k] + s3
                continue
            elif 0 == backtrack[i][j][k]:
                i -= 1
                j -= 1
                k -= 1
                s1 = u[i] + s1
                s2 = v[j] + s2
                s3 = w[k] + s3
        return s1, s2, s3
    
    def computeCLS(self, seq1, seq2, seq3):
        maxSocre, backtrack, i, j, k = self.LCSBackTrack(seq1, seq2, seq3)
        s1, s2, s3 = self.OutputLCS(backtrack, seq1, seq2, seq3, i, j, k)
        return maxSocre, s1, s2, s3

if __name__ == "__main__":
    MultipleLCS()