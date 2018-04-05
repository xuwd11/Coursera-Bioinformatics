# python3

import sys
import queue

'''
We now have an algorithm for solving the Limb Length Problem. For each j, we can compute LimbLength(j) by finding the minimum value 
of (Di,j+Dj,k−Di,k)/2(Di,j+Dj,k−Di,k)/2 over all pairs of leaves i and k.

Code Challenge: Solve the Limb Length Problem.
     Input: An integer n, followed by an integer j between 0 and n - 1, followed by a space-separated
     additive distance matrix D (whose elements are integers).
     Output: The limb length of the leaf in Tree(D) corresponding to row j of this distance
     matrix (use 0-based indexing).

Sample Input:
4
1
0	13	21	22
13	0	12	13
21	12	0	13
22	13	13	0
Sample Output:
2
'''

class LimbLength:
    def __init__(self):
        n, j, distMatrix = self._input()
        limbLength = self.calculateLimbLength(n, j, distMatrix)
        print(limbLength)

    def _input(self):
        data = sys.stdin.read().strip().split('\n')
        n = int(data[0])
        j = int(data[1])
        distMatrix = [[0]*n for _ in range(n)]
        for i in range(n):
            d = data[i+2].split()
            for k in range(n):
                distMatrix[i][k] = int(d[k])
        return n, j, distMatrix
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        n = int(data[0])
        j = int(data[1])
        distMatrix = [[0]*n for _ in range(n)]
        for i in range(n):
            d = data[i+2].split()
            for k in range(n):
                distMatrix[i][k] = int(d[k])
        return n, j, distMatrix

    def saveResult(self, result):
        f = open('result.txt', 'w')

    def printDistMatrix(self, distMatrix):
        for d in distMatrix:
            print(' '.join([str(i) for i in d]))

    def calculateLimbLength(self, n, j, distMatrix):
        limbLength = float('inf')
        if j > 0:
            i = j - 1
        else:
            i = j + 1
        for k in range(n):
            if i != k and k != j:
                currLength = (distMatrix[i][j] + distMatrix[j][k] - distMatrix[i][k])//2
                if currLength < limbLength:
                    limbLength = currLength
        return limbLength

if __name__ == "__main__":
    LimbLength()