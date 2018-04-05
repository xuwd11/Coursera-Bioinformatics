# python3
import sys
from math import *
import numpy as np

'''
Decoding Problem: Find an optimal hidden path in an HMM given a string of its emitted symbols.
     Input: A string x emitted by an HMM (Σ, States, Transition, Emission).
     Output: A path π that maximizes the (unconditional) probability Pr(x, π) that x will be
     emitted given that the HMM follow over all possible paths through this HMM.
Implement the Viterbi algorithm solving the Decoding Problem.
     Input: A string x, followed by the alphabet from which x was constructed,
     followed by the states States, transition matrix Transition, and emission matrix
     Emission of an HMM (Σ, States, Transition, Emission).
     Output: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.

Note: You may assume that transitions from the initial state occur with equal probability.

Sample Input:
xyxzzxyxyy
--------
x y z
--------
A B
--------
	A	B
A	0.641	0.359
B	0.729	0.271
--------
	x	y	z
A	0.117	0.691	0.192	
B	0.097	0.42	0.483
Sample Output:
AAABBAAAAA

We will apply a dynamic programming algorithm to solve the Decoding Problem. First, define sk,i as the product weight of 
an optimal path (i.e., a path with maximum product weight) from source to node (k, i). The Viterbi algorithm relies on the 
fact that the first i−1 edges of an optimal path from source to (k, i) must form an optimal path from source to (l, i−1) for 
some (unknown) state l. This observation yields the following recurrence:

sk,i = maxall states l{sl,i−1⋅(weight of edge between nodes(l,i−1) and (k,i))}
= maxall states l{sl,i−1⋅Weighti(l,k)}
= maxall states l{sl,i−1⋅transitionπi−1,πi⋅emissionπi(xi)}

For i = 1, every node in the leftmost column is connected to source, so that the above recurrence becomes

sk,1 = ssource⋅(weight of edge between source and (k,1))
= ssource⋅Weight0(source,k)
= ssource⋅transitionsource,k⋅emissionk(x1)

In order to initialize this recurrence, we set ssourcessource equal to 1. We can now compute the maximum product weight over 
all paths from source to sink as

ssink = maxall states lsl,n.
'''

class Decoding:
    def __init__(self):
        x, transionLog, emissionLog, stateDict = self.readFromFile()
        path = self.viterbi(x, transionLog, emissionLog, stateDict)
        print(path)
        f = open('result.txt', 'w')
        f.write(path)
        f.close()

    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().split()
        x = data[0]
        ind = [i for i in range(len(data)) if '--------' == data[i]]
        alphabet = data[ind[0]+1:ind[1]]
        states = data[ind[1]+1:ind[2]]
        stateDict = {i:states[i] for i in range(len(states))}
        transitionLog = {i:{k:log(float(data[ind[2]+len(states)+2+i*(len(states)+1)+k])) for k in range(len(states))} for i in range(len(states))}
        emissionLog = {i:{alphabet[k]:log(float(data[ind[3]+len(alphabet)+2+i*(len(alphabet)+1)+k])) for k in range(len(alphabet))} for i in range(len(states))}
        f.close()
        return x, transitionLog, emissionLog, stateDict

    def viterbi(self, x, transionLog, emissionLog, stateDict):
        n = len(x)
        l = len(transionLog)
        s = [[0 for _ in range(l)] for __ in range(n)]
        backTrack = [[0 for _ in range(l)] for __ in range(n)]
        for k in range(l):
            s[0][k] = log(1/l) + emissionLog[k][x[0]]
        for i in range(1, n):
            for k in range(l):
                currS = [s[i-1][kpre] + transionLog[kpre][k] + emissionLog[k][x[i]] for kpre in range(l)]
                ind = np.argmax(currS)
                backTrack[i][k] = ind
                s[i][k] = currS[ind]
        
        currState = np.argmax(s[n-1])
        stateList = [currState]
        for i in range(n-1, 0, -1):
            currState = backTrack[i][currState]
            stateList.insert(0, currState)
        path = ''.join([stateDict[state] for state in stateList])
        return path

if __name__ == '__main__':
    Decoding()