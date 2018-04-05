# python3
import sys
import numpy as np

'''
Solve the HMM Parameter Estimation Problem.
     Input: A string x of symbols emitted from an HMM, followed by the HMM's alphabet Σ,
     followed by a path π, followed by the collection of states of the HMM.
     Output: A transition matrix Transition followed by an emission matrix Emission that maximize
     Pr(x, π) over all possible transition and emission matrices.

Sample Input:
yzzzyxzxxx
--------
x y z
--------
BBABABABAB
--------
A B C
Sample Output:
	A	B	C
A	0.0	1.0	0.0
B	0.8	0.2	0.0
C	0.333	0.333	0.333
--------
	x	y	z
A	0.25	0.25	0.5
B	0.5 0.167	0.333
C	0.333	0.333	0.333

HMM Parameter Estimation Problem: Find optimal parameters explaining the emitted string and hidden path of an HMM.
     Input: A string x = x1 ... xn emitted by a k-state HMM with unknown transition and
     emission probabilities following a known hidden path π = π1 . . . − πn .
     Output: A transition matrix Transition and an emission matrix Emission that maximize
     Pr(x, π) over all possible transition and emission matrices.

If we know both xx and ππ, then we can compute empirical estimates for the transition and emission probabilities using a 
method similar to one we used for estimating parameters for profile HMMs. If Tl,kTl,k denotes the number of transitions 
from state ll to state kk in the hidden path ππ, then we can estimate the probability transitionl,ktransitionl,k by computing 
the ratio of Tl,k to the total number of transitions leaving state l,

transitionl,k=Tl,k/∑all states jTl,j.

Likewise, if Ek(b) denotes the number of times symbol b is emitted when the hidden path π is in state k, then we can estimate 
the probability emissionk(b) as the ratio of Ek(b) to the total number of emitted symbols from state k,

emissionk(b)=Ek(b)/∑all symbols c in the alphabetEk(c).

It turns out that the above two formulas for computing Transition and Emission result in parameters solving the HMM Parameter 
Estimation Problem.
'''

class HMMParameterEstimation:
    def __init__(self):
        x, path, alphabet, states = self.readFromFile()
        transition, emission = self.estimateParameters(x, path, alphabet, states)
        #print(transition, emission)
        self.saveParameters(transition, emission, alphabet, states)

    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().strip().split()
        ind = [i for i in range(len(data)) if '--------' == data[i]]
        x = data[0]
        alphabet = data[ind[0]+1:ind[1]]
        states = data[ind[2]+1:]
        path = data[ind[1]+1]
        f.close()
        return x, path, alphabet, states

    def estimateParameters(self, x, path, alphabet, states):
        alphabet2ind = {alphabet[i]:i for i in range(len(alphabet))}
        states2ind = {states[i]:i for i in range(len(states))}
        transition = np.zeros((len(states), len(states)), dtype = float)
        emission = np.zeros((len(states), len(alphabet)), dtype = float)
        
        for i in range(len(path)-1):
            transition[states2ind[path[i]], states2ind[path[i+1]]] += 1
        
        for i in range(len(x)):
            emission[states2ind[path[i]], alphabet2ind[x[i]]] += 1
        
        for i in range(len(states)):
            sum1 = sum(transition[i,:])
            if 0 == sum1:
                transition[i,:] += 1/len(states)
            else:
                transition[i,:] /= sum1
            sum2 = sum(emission[i,:])
            if 0 == sum2:
                emission[i,:] += 1/len(alphabet)
            else:
                emission[i,:] /= sum2
        
        return transition, emission

    def saveParameters(self, transition, emission, alphabet, states):
        f = open('result.txt', 'w')
        print(' '.join([' ']+states))
        f.write('\t'+'\t'.join(states)+'\n')
        for i in range(len(states)):
            print(' '.join([states[i]]+list(map(str, transition[i, :]))))
            f.write('\t'.join([states[i]]+list(map(str, transition[i, :])))+'\n')
        print('--------')
        f.write('--------\n')
        print(' '.join([' ']+alphabet))
        f.write('\t'+'\t'.join(alphabet)+'\n')
        for i in range(len(states)):
            print(' '.join([states[i]]+list(map(str, emission[i, :]))))
            f.write('\t'.join([states[i]]+list(map(str, emission[i, :])))+'\n')
        f.close()

if __name__ == '__main__':
    HMMParameterEstimation()