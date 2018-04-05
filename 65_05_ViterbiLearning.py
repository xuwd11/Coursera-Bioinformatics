# python3
import sys
from math import *
import numpy as np

'''
Implement Viterbi learning for estimating the parameters of an HMM.
     Input: A number of iterations j, followed by a string x of symbols emitted by an HMM,
     followed by the HMM's alphabet Σ, followed by the HMM's states, followed by initial transition
     and emission matrices for the HMM.
     Output: Emission and transition matrices resulting from applying Viterbi learning for j iterations.

Sample Input:
100
--------
zyzxzxxxzz
--------
x y z
--------
A B
--------
	A	B
A	0.599	0.401	
B	0.294	0.706	
--------
	x	y	z
A	0.424	0.367	0.209	
B	0.262	0.449	0.289
Sample Output:
	A	B
A	0.5	0.5	
B	0.0	1.0	
--------
	x	y	z
A	0.333	0.333	0.333	
B	0.4	0.1	0.5

HMM Parameter Learning Problem: Estimate the parameters of an HMM explaining an emitted string.
     Input: A string x = x1 ... xn emitted by an HMM with unknown transition and emission probabilities.
     Output: A transition matrix Transition and an emission matrix Emission that maximize
     Pr(x, π) over all possible transition and emission matrices and over all hidden paths π.

Unfortunately, the HMM Parameter Learning Problem is intractable, and so we will instead develop a heuristic that is 
analogous to the Lloyd algorithm for k-means clustering. In that algorithm, we iterated two steps, “From Centers to 
Clusters”,

(Data, ?, Centers) → HiddenVector,

and “From Clusters to Centers”,

(Data, HiddenVector, ?) → Centers.

As for HMM parameter estimation, we begin with an initial random guess for Parameters. Then, we use the Viterbi algorithm 
to find the optimal hidden path π:

(x, ?, Parameters) → π

Once we know π, we will question our original choice of Parameters and apply our solution to the HMM Parameter Estimation 
Problem to update Parameters based on x and π:

(x, π, ?) → Parameters′

We then iterate over these two steps, hoping that the estimated parameters are getting closer and closer to the parameters 
solving the HMM Parameter Learning Problem:

(x, ?, Parameters) → (x, π, Parameters) → (x, π, ?) → (x, π, Parameters′) → (x, ?, Parameters′) → (x, π′, Parameters′) → 
(x, π′, ?) → (x, π′, Parameters′′) → . . .

This approach to learning the HMM’s parameters is called Viterbi learning.
'''

class ViterbiLearning:
    def __init__(self):
        x, transition, emission, alphabet, states, iterNo = self.readFromFile()
        transition, emission = self.viterbiLearning(x, transition, emission, alphabet, states, iterNo)
        self.saveTransitionAndEmission(alphabet, states, transition, emission)

    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().strip().split()
        iterNo = int(data[0])
        x = data[2]
        ind = [i for i in range(len(data)) if '--------' == data[i]]
        alphabet = data[ind[1]+1:ind[2]]
        states = data[ind[2]+1:ind[3]]
        transition = np.zeros((len(states), len(states)), dtype = float)
        emission = np.zeros((len(states), len(alphabet)), dtype = float)
        for i in range(len(states)):
            transition[i, :] = [float(d) for d in data[ind[3]+len(states)+2+i*(len(states)+1):ind[3]+len(states)+1+(i+1)*(len(states)+1)]]
            emission[i, :] = [float(d) for d in data[ind[4]+len(alphabet)+2+i*(len(alphabet)+1):ind[4]+len(alphabet)+1+(i+1)*(len(alphabet)+1)]]
        return x, transition, emission, alphabet, states, iterNo
    
    def decode(self, xList, transition, emission):
        n, l = len(xList), transition.shape[0]
        score = [[1. for _ in range(l)] for __ in range(n)]
        backtrack = [[None for _ in range(l)] for __ in range(n)]
        for k in range(l):
            score[0][k] = emission[k, xList[0]]/l
        for i in range(1, n):
            for k in range(l):
                currScore = [score[i-1][kpre] * transition[kpre, k] * emission[k, xList[i]] for kpre in range(l)]
                ind = np.argmax(currScore)
                backtrack[i][k] = ind
                score[i][k] = currScore[ind]
        currState = np.argmax(score[n-1])
        pathList = [currState]
        for i in range(n-1, 0, -1):
            currState = backtrack[i][currState]
            pathList.insert(0, currState)
        return pathList

    def estimateParameters(self, xList, pathList, transition, emission):
        transition = np.zeros_like(transition, dtype = float)
        emission = np.zeros_like(emission, dtype = float)
        n, l = emission.shape
        for i in range(len(pathList)-1):
            transition[pathList[i], pathList[i+1]] += 1
        for i in range(len(xList)):
            emission[pathList[i], xList[i]] += 1
        for i in range(n):
            sum1 = sum(transition[i,:])
            if 0 == sum1:
                transition[i,:] += 1/n
            else:
                transition[i,:] /= sum1
            sum2 = sum(emission[i,:])
            if 0 == sum2:
                emission[i,:] += 1/l
            else:
                emission[i,:] /= sum2
        return transition, emission

    def viterbiLearning(self, x, transition, emission, alphabet, states, iterNo):
        x2ind, ind2x = {alphabet[i]:i for i in range(len(alphabet))}, {i:alphabet[i] for i in range(len(alphabet))}
        path2ind, ind2path = {states[i]:i for i in range(len(states))}, {i:states[i] for i in range(len(states))}
        xList = [x2ind[x[i]] for i in range(len(x))]
        for _ in range(iterNo):
            pathList = self.decode(xList, transition, emission)
            transition, emission = self.estimateParameters(xList, pathList, transition, emission)
        return transition, emission
    
    def saveTransitionAndEmission(self, alphabet, states, fullTransition, emission):
        f = open('result.txt', 'w')
        print(' '.join([' '] + states))
        f.write('\t'+'\t'.join(states) + '\n')
        for i in range(fullTransition.shape[0]):
            print(' '.join([states[i]] + list(map(str, fullTransition[i, :]))))
            f.write('\t'.join([states[i]] + list(map(str, fullTransition[i, :])))+'\n')
        print('--------')
        f.write('--------'+'\n')
        print(' '.join([' '] + alphabet))
        f.write('\t'+'\t'.join(alphabet)+'\n')
        for i in range(emission.shape[0]):
            print(' '.join([states[i]] + list(map(str, emission[i, :]))))
            f.write('\t'.join([states[i]] + list(map(str, emission[i, :])))+'\n')
        f.close()

if __name__ == '__main__':
    ViterbiLearning()