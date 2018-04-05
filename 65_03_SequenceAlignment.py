# python3
import sys
import numpy as np

'''
Solve the Sequence Alignment with Profile HMM Problem.
     Input: A string x followed by a threshold θ and a pseudocount σ, followed by an
     alphabet Σ, followed by a multiple alignment Alignment whose strings are formed from Σ. 
     Output: An optimal hidden path emitting x in HMM(Alignment, θ, σ).

Sample Input:
AEFDFDC
--------
0.4 0.01
--------
A B C D E F
--------
ACDEFACADF
AFDA---CCF
A--EFD-FDC
ACAEF--A-C
ADDEFAAADF
Sample Output:
M1 D2 D3 M4 M5 I5 M6 M7 M8
'''

class SequenceAlignment:
    def __init__(self):
        x, theta, pseudocount, alphabet, alignment = self.readFromFile()
        transition, emission = self.profile(theta, pseudocount, alphabet, alignment)
        states = self.getAllStates(transition.shape[0])
        transitionLog, emissionLog = np.log(transition), np.log(emission)
        path = self.findHiddenPath(x, transitionLog, emissionLog, alphabet)
        self.savePath(path, states)

    def readFromFile(self):
        f = open('input.txt', 'r')
        data = f.read().split()
        x = data[0]
        theta, pseudocount = float(data[2]), float(data[3])
        ind = [i for i in range(len(data)) if '--------' == data[i]]
        alphabet = data[ind[1]+1:ind[2]]
        alignment = np.array([[*s] for s in data[ind[2]+1:]])
        f.close()
        return x, theta, pseudocount, alphabet, alignment

    def profile(self, theta, pseudocount, alphabet, alignment): # with pseudocounts
        alphabetDict = ({alphabet[i]:i for i in range(len(alphabet))}, {i:alphabet[i] for i in range(len(alphabet))})
        n, l = alignment.shape
        threshold = theta * n
        kept = [True] * l
        for i in range(l):
            if sum('-' == alignment[:, i]) >= threshold:
                kept[i] = False
        levels = [0 for _ in range(l)]
        for i in range(l):
            levels[i] = levels[i-1]
            if kept[i]:
                levels[i] += 1

        def getIndex(level, state):
            if 0 == level:
                if 0 == state:
                    return 1
                else:
                    return 0
            if 0 != state:
                return 3*level-2+state
            else:
                return 3*level+1

        transition = np.zeros((sum(kept)*3+3, 3), dtype = float)
        # 0: I; 1: M; 2: D/E
        emission = np.zeros((sum(kept)*3+3, len(alphabet)), dtype = float)

        for i in range(n):
            lastLevel = 0
            lastState = -1
            lastInd = getIndex(lastLevel, lastState)
            for j in range(l):
                currLevel = levels[j]
                if kept[j]:
                    currState = 2 if '-' == alignment[i, j] else 1
                    currInd = getIndex(currLevel, currState)
                    transition[lastInd, currState] += 1
                    if 1 == currState:
                        emission[currInd, alphabetDict[0][alignment[i, j]]] += 1
                    lastInd = currInd
                else:
                    if '-' != alignment[i, j]:
                        currState = 0
                        currInd = getIndex(currLevel, currState)
                        transition[lastInd, currState] += 1
                        emission[currInd, alphabetDict[0][alignment[i, j]]] += 1
                        lastInd = currInd
            transition[lastInd, 2] += 1
        
        for i in range(transition.shape[0]):
            sum1 = sum(transition[i, :])
            if 0 != sum1:
                transition[i, :] /= sum1
            sum2 = sum(emission[i, :])
            if 0 != sum2:
                emission[i, :] /= sum2
        
        # Add pseudocounts
        for i in range(transition.shape[0]-4):
            transition[i, :] += pseudocount
            transition[i, :] /= sum(transition[i, :])
        for i in range(-4, -1):
            transition[i, 0] += pseudocount
            transition[i, 2] += pseudocount
            transition[i, :] /= sum(transition[i, :])

        for i in range(1, emission.shape[0]-1):
            if 0 != i%3:
                emission[i, :] += pseudocount
                emission[i, :] /= sum(emission[i, :])
            
        return transition, emission

    def getAllStates(self, n):
        # n: number of states in total
        states = ['' for _ in range(n)]
        states[0] = 'S'
        states[-1] = 'E'
        states[1] = 'I0'
        s = ('M', 'D', 'I')
        for i in range(2, n-1):
            states[i] = s[(i+1)%3] + str((i+1)//3)
        return states

    def getFullTransition(self, transition):
        fullTransition = np.zeros((transition.shape[0], transition.shape[0]), dtype = float)
        fullTransition[0, 1:4] = transition[0, :]

        for i in range(1, transition.shape[0]-4):
            if 1 == i%3:
                fullTransition[i, i:i+3] = transition[i, :]
            elif 2 == i%3:
                fullTransition[i, i+2:i+5] = transition[i, :]
            else:
                fullTransition[i, i+1:i+4] = transition[i, :]
        for i in range(-4, -1):
            fullTransition[i, -2] = transition[i, 0]
            fullTransition[i, -1] = transition[i, -1]

        return fullTransition

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

    def findHiddenPath(self, x, transitionLog, emissionLog, alphabet):
        alphabet2ind = {alphabet[i]:i for i in range(len(alphabet))}
        l = len(x)
        n = transitionLog.shape[0]
        score = [[-np.inf for _ in range(l+1)] for __ in range(n)]
        backtrack = [[None for _ in range(l+1)] for __ in range(n)]
        score[3][0] = transitionLog[0, 2]
        backtrack[3][0] = (0, 0)
        for i in range(6, n, 3):
            score[i][0] = score[i-3][0] + transitionLog[i-3, 2]
            backtrack[i][0] = (i-3, 0)
        score[1][1] = transitionLog[0, 0] + emissionLog[1, alphabet2ind[x[0]]]
        backtrack[1][1] = (0, 0)
        score[2][1] = transitionLog[0, 1] + emissionLog[2, alphabet2ind[x[0]]]
        backtrack[2][1] = (0, 0)
        score[3][1] = score[1][1] + transitionLog[1, 2]
        backtrack[3][1] = (1, 1)
        for j in range(2, l+1):
            score[1][j] = score[1][j-1] + transitionLog[1, 0] + emissionLog[1, alphabet2ind[x[j-1]]]
            backtrack[1][j] = (1, j-1)
            score[2][j] = score[1][j-1] + transitionLog[1, 1] + emissionLog[2, alphabet2ind[x[j-1]]]
            backtrack[2][j] = (1, j-1)
            score[3][j] = score[1][j] + transitionLog[1, 2]

        for i in range(4, n-1):
            if 1 == i%3: # I
                score[i][1] = score[i-1][0] + transitionLog[i-1, 0] + emissionLog[i, alphabet2ind[x[0]]]
                backtrack[i][1] = (i-1, 0)
                for j in range(2, l+1):                                
                    candBack = ((i, j-1), (i-2, j-1), (i-1, j-1))
                    candScore = (score[i][j-1] + transitionLog[i, 0] + emissionLog[i, alphabet2ind[x[j-1]]], \
                    score[i-2][j-1] + transitionLog[i-2, 0] + emissionLog[i, alphabet2ind[x[j-1]]], \
                    score[i-1][j-1] + transitionLog[i-1, 0] + emissionLog[i, alphabet2ind[x[j-1]]])
                    ind = np.argmax(candScore)
                    score[i][j] = candScore[ind]
                    backtrack[i][j] = candBack[ind]
            elif 2 == i%3: # M
                score[i][1] = score[i-2][0] + transitionLog[i-2, 1] + emissionLog[i, alphabet2ind[x[0]]]
                backtrack[i][1] = (i-2, 0)
                for j in range(2, l+1):
                    candBack = ((i-1, j-1), (i-3, j-1), (i-2, j-1))
                    candScore = (score[i-1][j-1] + transitionLog[i-1, 1] + emissionLog[i, alphabet2ind[x[j-1]]], \
                    score[i-3][j-1] + transitionLog[i-3, 1] + emissionLog[i, alphabet2ind[x[j-1]]], \
                    score[i-2][j-1] + transitionLog[i-2, 1] + emissionLog[i, alphabet2ind[x[j-1]]])
                    ind = np.argmax(candScore)
                    score[i][j] = candScore[ind]
                    backtrack[i][j] = candBack[ind]
            else: # D
                for j in range(1, l+1):
                    candBack = ((i-2, j), (i-4, j), (i-3, j))
                    candScore = (score[i-2][j] + transitionLog[i-2, 2], \
                    score[i-4][j] + transitionLog[i-4, 2], \
                    score[i-3][j] + transitionLog[i-3, 2])
                    ind = np.argmax(candScore)
                    score[i][j] = candScore[ind]
                    backtrack[i][j] = candBack[ind]
        candBack = ((n-2, l), (n-4, l), (n-3, l))
        candScore = (score[n-2][l] + transitionLog[n-2, 2], \
        score[n-4][l] + transitionLog[n-4, 2], \
        score[n-3][l] + transitionLog[n-3, 2])
        ind = np.argmax(candScore)
        score[n-1][l] = candScore[ind]
        backtrack[n-1][l] = candBack[ind]

        # backtrack
        path = []
        currPos = backtrack[n-1][l]
        while 0 != currPos[0]:
            path.insert(0, currPos)
            currPos = backtrack[currPos[0]][currPos[1]]
        return path
    
    def savePath(self, path, states):
        print(' '.join([states[path[i][0]] for i in range(len(path))]))
        f = open('result.txt', 'w')
        f.write(str(' '.join([states[path[i][0]] for i in range(len(path))])))
        f.close()

if __name__ == '__main__':
    SequenceAlignment()