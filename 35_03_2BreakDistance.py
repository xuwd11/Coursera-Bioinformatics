# python3

import sys
import numpy as np
import copy

'''
Solve the 2-Break Distance Problem.
     Input: Genomes P and Q.
     Output: The 2-break distance d(P, Q).

Sample Input:
(+1 +2 +3 +4 +5 +6)
(+1 -3 -6 -5)(+2 -4)
Sample Output:
3

Cycle Theorem: Given genomes P and Q, any 2-break applied to P can increase Cycles(P, Q) by at most 1.
2-Break Distance Theorem: The 2-break distance between genomes P and Q is equal to Blocks(P, Q)− Cycles(P, Q).
'''

class TwoBreakDistance:
    def __init__(self):
        genomes = self._inputGenomes()
        #genomes = self.readGenomesFromFile()
        dist = self.calculate2BreakDistance(genomes[0], genomes[1])
        print(dist)
    
    def _inputGenomes(self):
        data = sys.stdin.read().strip().split('\n')
        genomes = []
        for g in data:
            g = g.split(')(')
            genome = []
            for d in g:
                d = d.split()
                if d[0][-1] != ')':
                    genome.append([int(d[0][1:] if '('==d[0][0] else d[0])] + [int(e) for e in d[1:-1]] +\
                    [int(d[-1][:-1] if ')'==d[-1][-1] else d[-1])])
                else:
                    genome.append([int(d[0][:-1])])
            genomes.append(genome)
        return genomes
    
    def readGenomesFromFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        genomes = []
        for g in data:
            g = g.split(')(')
            genome = []
            for d in g:
                d = d.split()
                genome.append([int(d[0][1:] if '('==d[0][0] else d[0])] + [int(e) for e in d[1:-1]] +\
                [int(d[-1][:-1] if ')'==d[-1][-1] else d[-1])])
            genomes.append(genome)
        return genomes

    def chromosomeToCycle(self, chromosome):
        l = len(chromosome)
        nodes = [0]*(2*l)
        for j in range(l):
            i = chromosome[j]
            if i > 0:
                nodes[2*j] = 2*i-1
                nodes[2*j+1] = 2*i
            else:
                nodes[2*j] = -2*i
                nodes[2*j+1] = -2*i-1
        return nodes
        
    def cycleToChromosome(self, nodes):
        l = len(nodes) // 2
        chromosome = [0]*l
        for j in range(l):
            if nodes[2*j] < nodes[2*j+1]:
                chromosome[j] = nodes[2*j+1]//2
            else:
                chromosome[j] = -nodes[2*j]//2
        return chromosome
    
    def printChromosome(self, chromosome):
        print('('+' '.join(['+'+str(e) if e>0 else str(e) for e in chromosome])+')')

    def coloredEdges(self, genome):
        edges = set()
        for chromosome in genome:
            nodes = self.chromosomeToCycle(chromosome)
            nodes.append(nodes[0])
            for j in range(len(chromosome)):
                edges.add((nodes[2*j+1], nodes[2*j+2]))
        return edges
        
    def calculate2BreakDistance(self, P, Q):
        blocks = sum([len(a) for a in P])
        edges = self.coloredEdges(P).union(self.coloredEdges(Q))
        parent = dict()
        rank = dict()
        for e in edges:
            parent[e[0]] = e[0]
            parent[e[1]] = e[1]
            rank[e[0]] = 0
            rank[e[1]] = 0

        def findParent(i):
            if i != parent[i]:
                parent[i] = findParent(parent[i])
            return parent[i]
        
        def union(i, j):
            i_id = findParent(i)
            j_id = findParent(j)
            if i_id == j_id:
                return
            if rank[i_id] > rank[j_id]:
                parent[j_id] = i_id
            else:
                parent[i_id] = j_id
                if rank[i_id] == rank[j_id]:
                    rank[j_id] += 1
        
        for e in edges:
            union(e[0], e[1])

        nodesSets = set()

        for e in edges:
            id = findParent(e[0])
            nodesSets.add(id)
        
        cycles = len(nodesSets)
        dist = blocks - cycles
        return dist

    def printGenome(self, genome):
        result = ''
        for chromosome in genome:
            result += '('+' '.join(['+'+str(e) if e>0 else str(e) for e in chromosome])+')'
        print(result)
        
if __name__ == "__main__":
    TwoBreakDistance()