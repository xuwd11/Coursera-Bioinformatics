# python3
import sys

'''
Solve the Probability of an Outcome Given a Hidden Path Problem.
     Input: A string x, followed by the alphabet from which x was constructed, followed by
     a hidden path π, followed by the states States and emission matrix Emission of an HMM
     (Σ, States, Transition, Emission).
     Output: The conditional probability Pr(x|π) that x will be emitted given that the HMM
     follows the hidden path π.

Note: You may assume that transitions from the initial state occur with equal probability.

Sample Input:
zzzyxyyzzx
--------
x y z
--------
BAAAAAAAAA
--------
A B
--------
	x	y	z
A	0.176	0.596	0.228
B	0.225	0.572	0.203
Sample Output:
3.59748954746e-06

To compute Pr(x|π) for a general HMM, we will write Pr(xi|πi) to denote the emission probability emissionπi(xi) that symbol 
xi was emitted given that the HMM was in state πi. As a result, for a given path π, the HMM emits a string x with probability 
equal to the product of emission probabilities along that path,

Pr(x|π)=∏i=1nPr(xi|πi)=∏i=1nemissionπi(xi).
'''

class ProbOutcomeGivenPath:
    def __init__(self):
        x, path, emission = self.readFromFile()
        P = self.calculatePrxpi(x, path, emission)
        print(P)
        f = open('result.txt', 'w')
        f.write(str(P))
        f.close()

    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().split()
        x = data[0]
        path = data[6]
        emission = {'A':{'x':float(data[-7]), 'y':float(data[-6]), 'z':float(data[-5])}, 'B':{'x':float(data[-3]), 'y':float(data[-2]), 'z':float(data[-1])}}
        f.close()
        return x, path, emission

    def calculatePrxpi(self, x, path, emission):
        P = 1
        for i in range(len(x)):
            P *= emission[path[i]][x[i]]
        return P
    
if __name__ == '__main__':
    ProbOutcomeGivenPath()