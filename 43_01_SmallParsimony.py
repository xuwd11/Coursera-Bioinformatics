# python3

import sys
import queue
import numpy as np

'''
Implement SmallParsimony to solve the Small Parsimony Problem.
     Input: An integer n followed by an adjacency list for a rooted binary tree with n leaves labeled by DNA strings.
     Output: The minimum parsimony score of this tree, followed by the adjacency list of the tree corresponding to labeling
     internal nodes by DNA strings in order to minimize the parsimony score of the tree.

Note: Remember to run SmallParsimony on each individual index of the strings at the leaves of the tree.

Sample Input:
4
4->CAAATCCC
4->ATTGCGAC
5->CTGCGCTG
5->ATGGACGA
6->4
6->5
Sample Output:
16
ATTGCGAC->ATAGCCAC:2
ATAGACAA->ATAGCCAC:2
ATAGACAA->ATGGACTA:2
ATGGACGA->ATGGACTA:1
CTGCGCTG->ATGGACTA:4
ATGGACTA->CTGCGCTG:4
ATGGACTA->ATGGACGA:1
ATGGACTA->ATAGACAA:2
ATAGCCAC->CAAATCCC:5
ATAGCCAC->ATTGCGAC:2
ATAGCCAC->ATAGACAA:2
CAAATCCC->ATAGCCAC:5

SmallParsimony(T, Character)
    for each node v in tree T
        Tag(v) ← 0
        if v is a leaf
            Tag(v) ← 1
            for each symbol k in the alphabet
                if Character(v) = k
                    sk(v) ← 0
                else
                    sk(v) ← ∞
    while there exist ripe nodes in T
        v ← a ripe node in T
        Tag(v) ← 1
        for each symbol k in the alphabet
            sk(v) ← minimumall symbols i {si(Daughter(v))+δi,k} + minimumall symbols j {sj(Son(v))+δj,k}
    return minimum over all symbols k {sk(v)}
'''

class SmallParsimony:
    def __init__(self):
        n, adjC, adjP, adj, nodes = self._input()
        s = self.runSmallParsimony(n, adjC, adjP, adj, nodes)
        self.printResults(s, adj, nodes)
        

    def _input(self):
        data = sys.stdin.read().strip().split('\n')
        n = int(data[0])
        adjC = [[] for _ in range(n)]
        adjP = [[] for _ in range(n)]
        adj = [dict() for _ in range(n)]
        nodes = ['' for _ in range(n)]
        currNode = 0
        for d in data[1:]:
            d = d.split('->')
            p = int(d[0])
            try:
                c = int(d[1])
            except:
                c = currNode
                nodes[c] = d[1]
                currNode += 1
            if p > len(adjC)-1 or c > len(adjC)-1:
                adjC.extend([[] for _ in range(max([p,c])-len(adjC)+1)])
                adjP.extend([[] for _ in range(max([p,c])-len(adjP)+1)])
                adj.extend([dict() for _ in range(max([p,c])-len(adj)+1)])
            adjC[p].append(c)
            adjP[c].append(p)
            adj[p][c] = 0
            adj[c][p] = 0
        nodes.extend(['' for _ in range(len(adjC)-n)])
        #print(adjC, adjP, adj, nodes)
        return n, adjC, adjP, adj, nodes
    
    def printResults(self, s, adj, nodes):
        print(s)
        for i, d in enumerate(adj):
            for j, w in d.items():
                print(nodes[i]+'->'+nodes[j]+':'+str(w))
    
    def charIndConversion(self):
        char2ind = {'A':0, 'C':1, 'G':2, 'T':3}
        ind2char = {0:'A', 1:'C', 2:'G', 3:'T'}
        return char2ind, ind2char
    
    def singleSmallParsimony(self, n, adjC, adjP, adj, nodes, char2ind, ind2char, charInd):
        s = [[np.inf]*4 for _ in range(len(adjC))]
        backtrack = [[(-1, -1) for _ in range(4)] for __ in range(len(adjC))]
        processed = [0 for _ in range(len(adjC))]
        ripe = set()
        for i in range(n):
            s[i][char2ind[nodes[i][charInd]]] = 0
            processed[i] = 1
            if len(adjP[i]) > 0:
                ripe.add(adjP[i][0])
        
        while len(ripe) > 0:
            v = ripe.pop()
            for k in range(4):
                l = [s[adjC[v][0]][i] + (0 if k == i else 1) for i in range(4)]
                r = [s[adjC[v][1]][i] + (0 if k == i else 1) for i in range(4)]
                largmin = np.argmin(l)
                rargmin = np.argmin(r)
                backtrack[v][k] = (largmin, rargmin)
                s[v][k] = l[largmin] + r[rargmin]
            processed[v] = 1
            if len(adjP[v]) > 0 and all([processed[u] for u in adjC[adjP[v][0]]]):
                ripe.add(adjP[v][0])
        
        ind = np.argmin(s[v])
        nodes[v] += ind2char[ind]
        smin = s[v][ind]

        #dist = [np.inf] * len(adj)
        q = queue.Queue()
        #dist[v] = 0
        q.put((v, ind))
        while not q.empty():
            v, k = q.get()
            if len(adjC[v]) > 0:
                u, w = adjC[v]
                l, r = backtrack[v][k]
                
                if k != l:
                    adj[v][u] += 1
                    adj[u][v] += 1
                if k != r:
                    adj[v][w] += 1
                    adj[w][v] += 1
                if len(adjC[u]) > 0:
                    nodes[u] += ind2char[l]
                    nodes[w] += ind2char[r]
                    q.put((u, l))
                    q.put((w, r))        
        return smin
    
    def runSmallParsimony(self, n ,adjC, adjP, adj, nodes):
        char2ind, ind2char = self.charIndConversion()
        s = 0
        for i in range(len(nodes[0])):
            s += self.singleSmallParsimony(n, adjC, adjP, adj, nodes, char2ind, ind2char, i)
        return s

if __name__ == "__main__":
    SmallParsimony()