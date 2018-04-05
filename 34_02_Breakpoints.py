# python3

import sys
import numpy as np
import copy

'''
Find the number of breakpoints in a permutation.
     Input: A permutation.
     Output: The number of breakpoints in this permutation.

Sample Input:
(+3 +4 +5 -12 -8 -7 -6 +1 +2 +10 +9 -11 +13 +14)
Sample Output:
8

Consecutive elements (pi pi+1) in permutation P = (p1 ... pn) form an adjacency if pi+1 − pi is equal to 1. By definition, for any 
positive element k < n, both (k k + 1) and (−(k + 1) −k) are adjacencies. If pi+1 − pi is not equal to 1, then we say that (pi pi+1) 
is a breakpoint.
Because any pair of consecutive elements of a permutation form either a breakpoint or adjacency, we have the following identity for 
any permutation P of length n:

    Adjacencies(P) + Breakpoints(P) = n + 1

This formula implies that a permutation on n elements may have at most n + 1 adjacencies.
'''

class Breakpoints:
    def __init__(self):
        #self._input()
        self.readFile()
        count = self.countBreakpoints(self.P)
        print(count)     

    def _input(self):
        data = sys.stdin.read().strip().split()
        self.P = [int(data[0][1:])] + [int(e) for e in data[1:-1]] + [int(data[-1][:-1])]
    
    def readFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        data = data[0].split()
        self.P = [int(data[0][1:])] + [int(e) for e in data[1:-1]] + [int(data[-1][:-1])] 

    def countBreakpoints(self, P):
        count = 0
        P = copy.deepcopy([0] + P + [len(P)+1])
        for i in range(len(P)-1):
            if 1 != P[i+1] - P[i]:
                count += 1
        return count
        
if __name__ == "__main__":
    Breakpoints()