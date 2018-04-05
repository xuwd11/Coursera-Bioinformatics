# python3

import sys
import numpy as np

'''
Edit Distance Problem: Find the edit distance between two strings.
     Input: Two strings.
     Output: The edit distance between these strings.

Sample Input:
PLEASANTLY
MEANLY
Sample Output:
5
'''

class EditDistance:
    def __init__(self):
        self._input()
        distance = self.computeEditDistance(self.seq1, self.seq2)
        print(distance)

    def _input(self):
        self.seq1, self.seq2 = sys.stdin.read().strip().split('\n')
        self.seq1 = self.seq1.strip()
        self.seq2 = self.seq2.strip()
    
    def computeEditDistance(self, v, w):
        n = len(v)
        m = len(w)
        s = np.matrix(np.zeros((n+1)*(m+1), dtype = np.int).reshape((n+1, m+1)))
        for i in range(1, n+1):
            s[i, 0] = s[i-1, 0] + 1
        for j in range(1, m+1):
            s[0, j] = s[0, j-1] + 1
        for i in range(1, n+1):
            for j in range(1, m+1):
                s[i, j] = min([s[i-1, j] + 1, s[i, j-1] + 1, s[i-1, j-1] + 1 if v[i-1]!=w[j-1] else s[i-1, j-1]])
        return s[n, m]

if __name__ == "__main__":
    EditDistance()