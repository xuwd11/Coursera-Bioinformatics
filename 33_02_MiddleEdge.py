# python3

import sys
import math

'''
Solve the Middle Edge in Linear Space Problem (for protein strings).
     Input: Two amino acid strings.
     Output: A middle edge in the alignment graph in the form "(i, j) (k, l)", where (i, j) connects to (k, l). To compute scores, use the
     BLOSUM62 scoring matrix and a (linear) indel penalty equal to 5.

Sample Input:
     PLEASANTLY
     MEASNLY

Sample Output:
     (4, 3) (5, 4)
'''

class MiddleEdge:
    def __init__(self):
        #self._input()
        self.readFile()
        self.sMatrix = self.scoringMatrix()
        self.sigma = 5
        i1, j1, i2, j2 = self.findMiddleEdge(self.seq1, self.seq2, 0, self.n, 0, self.m, self.sMatrix, self.sigma)
        print('('+str(i1)+', '+str(j1)+') ('+str(i2)+', '+str(j2)+')')

    def _input(self):
        self.seq1, self.seq2 = sys.stdin.read().strip().split('\n')
        self.seq1 = self.seq1.strip()
        self.seq2 = self.seq2.strip()
        self.n = len(self.seq1)
        self.m = len(self.seq2)
    
    def readFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        self.seq1 = data[0]
        self.seq2 = data[1]
        self.n = len(self.seq1)
        self.m = len(self.seq2)
                
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
    
    def findMiddleEdge(self, v, w, top, bottom, left, right, sMatrix, sigma):
        middle = math.floor((right-left)/2)
        n = bottom - top
        m = right - left
        fromSource = [[0]*(n + 1), [0]*(n + 1)]
        currColumn = 0
        for i in range(1, n+1):
            fromSource[currColumn][i] = fromSource[currColumn][i-1] - sigma
        currColumn = 1 - currColumn
        for jChar in range(left, left+middle+2):
            if jChar > bottom-1:
                continue
            fromSource[currColumn][0] = fromSource[1-currColumn][0] - sigma
            for i in range(1, n+1):
                iChar = i + top - 1
                score1 = fromSource[currColumn][i-1] - sigma
                score2 = fromSource[1-currColumn][i] - sigma
                score3 = fromSource[1-currColumn][i-1] + sMatrix[v[iChar]][w[jChar]]
                fromSource[currColumn][i] = max(score1, score2, score3)
            currColumn = 1 - currColumn
        leftColumn1 = currColumn

        toSink = [[0]*(n + 1), [0]*(n + 1)]
        currColumn = 0
        for i in range(n-1, -1, -1):
            toSink[currColumn][i] = toSink[currColumn][i+1] - sigma
        currColumn = 1- currColumn
        for jChar in range(right-1, left+middle-1, -1):
            toSink[currColumn][n] = toSink[1-currColumn][n] - sigma
            for i in range(n-1, -1, -1):
                iChar = i + top -1
                score1 = toSink[currColumn][i+1] - sigma
                score2 = toSink[1-currColumn][i] - sigma
                score3 = toSink[1-currColumn][i+1] + sMatrix[v[iChar]][w[jChar]]
                toSink[currColumn][i] = max(score1, score2, score3)
            currColumn = 1 - currColumn
        leftColumn2 = 1 - currColumn

        length = [0]*(n+1)
        for i in range(n+1):
            length[i] = fromSource[leftColumn1][i] + toSink[leftColumn2][i]
        #print(length)
        iMax = max(range(len(length)), key = lambda x: length[x])
        i1 = top + iMax - 1
        j1 = left + middle
        if iMax == n:
            i2 = i1
            j2 = j1 + 1
        else:
            score1 = fromSource[1-leftColumn1][iMax] + toSink[1-leftColumn2][iMax]
            score2 = fromSource[leftColumn1][iMax+1] + toSink[leftColumn2][iMax+1]
            score3 = fromSource[1-leftColumn1][iMax+1] + toSink[1-leftColumn2][iMax+1]
            sMax = max(score1, score2, score3)
            if sMax == score3:
                i2 = i1 + 1
                j2 = j1 + 1
            elif sMax == score1:
                i2 = i1
                j2 = j1 + 1
            else:
                i2 = i1 + 1
                j2 = j1
        return i1, j1, i2, j2

if __name__ == "__main__":
    MiddleEdge()