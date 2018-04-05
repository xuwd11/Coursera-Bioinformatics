# python3
import sys
from math import *
import numpy as np

'''
Solve the Soft Decoding Problem.
     Input: A string x, followed by the alphabet Σ from which x was constructed,
     followed by the states States, transition matrix Transition, and emission matrix
     Emission of an HMM (Σ, States, Transition, Emission).
     Output: An |x| x |States| matrix whose (i, k)-th element holds the conditional probability Pr(πi = k|x).

Sample Input:
zyxxxxyxzz
--------
x y z
--------
A B
--------
	A	B
A	0.911	0.089
B	0.228	0.772
--------
	x	y	z
A	0.356	0.191	0.453 
B	0.040	0.467	0.493
Sample Output:
A	B 
0.5438	0.4562 
0.6492	0.3508 
0.9647	0.0353 
0.9936	0.0064 
0.9957	0.0043 
0.9891	0.0109 
0.9154	0.0846 
0.964	0.036 
0.8737	0.1263 
0.8167	0.1833

Soft Decoding Problem: Find the probability that an HMM was in a particular state at a particular moment given its emitted 
string.
     Input: A string x = x1 ... xn emitted by an HMM.
     Output: The conditional probability Pr(πi = k|x) that the HMM was in state k at step i
     given that it emitted x.

Recall that we computed forwardk,i when solving the Outcome Likelihood Problem.  If we reverse the edges of the Viterbi 
graph, we can compute 

backwardk,i=∑all states lbackwardl,i+1⋅Weighti(k,l).

Here we have that the reversed edge connecting (l, i+1) to (k, i) has weight Weighti(k,l)=transitionk,l⋅emissionl(xi+1).

Combining the forward-backward algorithm described in the lecture with our solution to the Outcome Likelihood Problem for 
computing Pr(x) yields that

Pr(πi=k|x)=Pr(πi=k,x)/Pr(x)=forwardk,i⋅backwardk,i/forward(sink)
'''

class SoftDecoding:
    def __init__(self):
        x, transition, emission, alphabet, states = self.readFromFile()
        Pr = self.softDecode(x, transition, emission, alphabet, states)
        self.savePr(Pr, states)

    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().strip().split()
        x = data[0]
        ind = [i for i in range(len(data)) if '--------' == data[i]]
        alphabet = data[ind[0]+1:ind[1]]
        states = data[ind[1]+1:ind[2]]
        transition = np.zeros((len(states), len(states)), dtype = float)
        emission = np.zeros((len(states), len(alphabet)), dtype = float)
        for i in range(len(states)):
            transition[i, :] = [float(d) for d in data[ind[2]+len(states)+2+i*(len(states)+1):ind[2]+len(states)+1+(i+1)*(len(states)+1)]]
            emission[i, :] = [float(d) for d in data[ind[3]+len(alphabet)+2+i*(len(alphabet)+1):ind[3]+len(alphabet)+1+(i+1)*(len(alphabet)+1)]]
        return x, transition, emission, alphabet, states

    def softDecode(self, x, transition, emission, alphabet, states):
        n = len(x)
        l = transition.shape[0]
        x2ind = {alphabet[i]:i for i in range(len(alphabet))}
        xList = [x2ind[x[i]] for i in range(len(x))]
        forward = [[0 for _ in range(l)] for __ in range(n)]
        backward = [[0 for _ in range(l)] for __ in range(n)]
        for k in range(l):
            forward[0][k] = emission[k, xList[0]]/l
        for i in range(1, n):
            for k in range(l):
                forward[i][k] = sum([forward[i-1][kpre]*transition[kpre, k]*emission[k, xList[i]] for kpre in range(l)])
        fsink = sum(forward[n-1])

        for k in range(l):
            backward[n-1][k] = 1
        for i in range(n-2, -1, -1):
            for k in range(l):
                backward[i][k] = sum([backward[i+1][kpre]*transition[k, kpre]*emission[kpre, xList[i+1]] for kpre in range(l)])
        
        Pr = np.zeros((n, l), dtype = float)
        for i in range(n):
            for k in range(l):
                Pr[i, k] = forward[i][k]*backward[i][k]/fsink

        return Pr

    def savePr(self, Pr, states):
        f = open('result.txt', 'w')
        print(' '.join(states))
        f.write('\t'.join(states)+'\n')
        for i in range(Pr.shape[0]):
            print(' '.join(list(map(str, Pr[i, :]))))
            f.write('\t'.join(list(map(str, Pr[i, :])))+'\n')
        f.close()

if __name__ == '__main__':
    SoftDecoding()