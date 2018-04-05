# python3

import sys
import numpy as np

'''
Find the length of a longest path in the Manhattan Tourist Problem.
     Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right. The two matrices are separated
     by the "-" symbol.
     Output: The length of a longest path from source (0, 0) to sink (n, m) in the n × m rectangular grid whose edges are defined by
     the matrices Down and Right.
Sample Input:
4 4
1 0 2 4 3
4 6 5 2 1
4 4 5 2 1
5 6 8 5 3
-
3 2 4 0
3 2 4 2
0 7 3 3
3 3 0 2
1 3 2 2
Sample Output:
34

Pseudocode:
    ManhattanTourist(n, m, Down, Right)
        s0, 0 ← 0
        for i ← 1 to n
            si, 0 ← si-1, 0 + downi, 0
        for j ← 1 to m
            s0, j ← s0, j−1 + right0, j
        for i ← 1 to n
            for j ← 1 to m
                si, j ← max{si - 1, j + downi, j, si, j - 1 + righti, j}
        return sn, m
'''

class ManhattanTourist:
    def __init__(self):
        self._input()
        print(self.longestPath(self.n, self.m, self.down, self.right))

    def _input(self):                    
        data = sys.stdin.read().strip().split('\n')
        n, m = [int(s) for s in data[0].strip().split()]
        down = np.matrix(';'.join(data[1:n+1]))
        right = np.matrix(';'.join(data[2+n:3+2*n]))
        self.n = n
        self.m = m
        self.down = down
        self.right = right
    
    def longestPath(self, n, m, down, right):
        s = np.matrix(np.arange((n+1)*(m+1)).reshape((n+1, m+1)))
        s[0, 0] = 0
        for i in range(1, n+1):
            s[i, 0] = s[i-1, 0] + down[i-1, 0]
        for j in range(1, m+1):
            s[0, j] = s[0, j-1] + right[0, j-1]
        for i in range(1, n+1):
            for j in range(1, m+1):
                s[i, j] = max([s[i-1, j] + down[i-1, j], s[i, j-1] + right[i, j-1]])
        return s[n, m]

if __name__ == "__main__":
    ManhattanTourist()