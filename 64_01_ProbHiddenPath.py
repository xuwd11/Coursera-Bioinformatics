# python3
import sys

'''
Solve the Probability of a Hidden Path Problem.
     Given: A hidden path π followed by the states States and transition matrix Transition of an HMM
     (Σ, States, Transition, Emission).
     Return: The probability of this path, Pr(π).

Note: You may assume that transitions from the initial state occur with equal probability.

Sample Input:
ABABBBAAAA
--------
A B
--------
	A	B
A	0.377	0.623
B	0.26	0.74
Sample Output:
0.000384928691755

For simplicity, we assume that in the beginning, the crooked dealer is equally likely to use the fair or biased coin, an 
assumption that is modeled by setting Pr(π0 → π1) = 1/2, where π0 is a “silent” initial state that does not emit any symbols. 
The probability of π is equal to the product of its transition probabilities,

Pr(π)=∏i=1nPr(πi−1→πi)=∏i=1ntransitionπi−1,πi.
'''

class ProbHiddenPath:
    def __init__(self):
        path, transition = self.readFromFile()
        P = self.calculateProb(path, transition)
        print(P)
        f = open('result.txt', 'w')
        f.write(str(P))
        f.close()

    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().split()
        path = data[0]
        transition = {'A':{'A':float(data[-5]), 'B':float(data[-4])}, 'B':{'A':float(data[-2]), 'B':float(data[-1])}}
        f.close()
        return path, transition

    def calculateProb(self, path, transition):
        P = 0.5
        for i in range(len(path)-1):
            P *= transition[path[i]][path[i+1]]
        return P        
    
if __name__ == '__main__':
    ProbHiddenPath()