# python3
import sys

'''
Outcome Likelihood Problem: Find the probability that an HMM emits a given string.
     Input: A string x emitted by an HMM (Σ, States, Transition, Emission).
     Output: The probability Pr(x) that the HMM emits x.
Solve the Outcome Likelihood Problem.
     Input: A string x, followed by the alphabet from which x was constructed,
     followed by the states States, transition matrix Transition, and emission matrix
     Emission of an HMM (Σ, States, Transition, Emission).
     Output: The probability Pr(x) that the HMM emits x.

Note: You may assume that transitions from the initial state occur with equal probability.

Sample Input:
xzyyzzyzyy
--------
x y z
--------
A B
--------
	A	B
A	0.303	0.697 
B	0.831	0.169 
--------
	x	y	z
A	0.533	0.065	0.402 
B	0.342	0.334	0.324
Sample Output:
1.1005510319694847e-06

We denote the total product weight of all paths from source to node (k,i)(k,i) in the Viterbi graph as forwardk,i; 
note that forwardsink is equal to Pr(x). To compute forwardk,iforwardk,i, we will divide all paths connecting 
source to (k, i) into |States| subsets, where each subset contains those paths that pass through node (l, i - 1) 
(with product weight forwardl,i−1) before reaching (k, i) for some l between 1 and |States|. Therefore, 
forwardk,i is the sum of |States| terms:

forwardk,i = ∑all states lforwardl,i−1⋅(weight of edge connecting (l,i−1) and (k,i))
= ∑all states lforwardl,i−1⋅Weighti(l,k)

Note that the only difference between the above recurrence and the Viterbi recurrence is that the maximization in the Viterbi 
algorithm has changed into a summation symbol. We can now solve the Outcome Likelihood Problem by computing forwardsink,
which is equal to

∑all states kforwardk,n.
'''

class OutcomeLikelihood:
    def __init__(self):
        x, transition, emission = self.readFromFile()
        prob = self.calculatePrx(x, transition, emission)
        print(prob)
    
    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().split()
        x = data[0]
        ind = [i for i in range(len(data)) if '--------' == data[i]]
        alphabet = data[ind[0]+1:ind[1]]
        states = data[ind[1]+1:ind[2]]
        stateDict = {i:states[i] for i in range(len(states))}
        transition = {i:{k:float(data[ind[2]+len(states)+2+i*(len(states)+1)+k]) for k in range(len(states))} for i in range(len(states))}
        emission = {i:{alphabet[k]:float(data[ind[3]+len(alphabet)+2+i*(len(alphabet)+1)+k]) for k in range(len(alphabet))} for i in range(len(states))}
        f.close()
        return x, transition, emission
    
    def calculatePrx(self, x, transition, emission):
        n = len(x)
        l = len(transition)
        forward = [[0 for _ in range(l)] for __ in range(n)]
        for k in range(l):
            forward[0][k] = 1/l*emission[k][x[0]]
        for i in range(1, n):
            for k in range(l):
                forward[i][k] = sum([forward[i-1][kpre]*transition[kpre][k]*emission[k][x[i]] for kpre in range(l)])
        return sum(forward[n-1])

if __name__ == '__main__':
    OutcomeLikelihood()