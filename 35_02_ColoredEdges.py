# python3

import sys
import numpy as np
import copy

'''
1. Implement ColoredEdges.
     Input: A genome P.
     Output: The collection of colored edges in the genome graph of P in the form (x, y).

Sample Input:
(+1 -2 -3)(+4 +5 -6)
Sample Output:
(2, 4), (3, 6), (5, 1), (8, 9), (10, 12), (11, 7)

ColoredEdges(P)
     Edges ← an empty set
     for each chromosome Chromosome in P
          Nodes ← ChromosomeToCycle(Chromosome)
          for j ← 1 to |Chromosome|
               add the edge (Nodes2j, Nodes2j +1) to Edges
     return Edges

2. Implement GraphToGenome.
     Input: The colored edges ColoredEdges of a genome graph.
     Output: The genome P corresponding to this genome graph.

Sample Input:
(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)
Sample Output:
(+1 -2 -3)(-4 +5 -6)

GraphToGenome(GenomeGraph)
     P ← an empty set of chromosomes
     for each cycle Nodes in GenomeGraph
          Nodes﻿ ← sequence of nodes in this cycle (starting from node 1)
          Chromosome ← CycleToChromosome(Nodes)
          add Chromosome to P
     return P
'''

class ColoredEdges:
    def __init__(self):
        '''
        self._inputGenome()
        edges = self.coloredEdges(self.genome)
        print(str(edges)[1:-1])
        '''
        #self._inputGraph()
        self.readGraphFromFile()
        genome = self.graphToGenome(self.edges)
        self.printGenome(genome)


    def _inputChromosome(self):
        data = sys.stdin.read().strip().split()
        self.chromosome = [int(data[0][1:])] + [int(e) for e in data[1:-1]] + [int(data[-1][:-1])]

    def _inputNodes(self):
        data = sys.stdin.read().strip().split()
        self.nodes = [int(data[0][1:])] + [int(e) for e in data[1:-1]] + [int(data[-1][:-1])]
    
    def _inputGenome(self):
        data = sys.stdin.read().strip().split(')(')
        genome = []
        for d in data:
            d = d.split()
            genome.append([int(d[0][1:] if '('==d[0][0] else d[0])] + [int(e) for e in d[1:-1]] +\
            [int(d[-1][:-1] if ')'==d[-1][-1] else d[-1])])
        self.genome = genome

    def _inputGraph(self):
        data = sys.stdin.read().strip().split('), (')
        edges = set()
        for d in data:
            d = d.split(',')
            edges.add((int(d[0][1:] if '('==d[0][0] else d[0]), int(d[1][:-1] if ')'==d[1][-1] else d[1])))
        self.edges = edges
    
    def readGraphFromFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        data = data[0].split('), (')
        edges = set()
        for d in data:
            d = d.split(',')
            edges.add((int(d[0][1:] if '('==d[0][0] else d[0]), int(d[1][:-1] if ')'==d[1][-1] else d[1])))
        self.edges = edges
    
    def readFile(self):
        f = open('input.txt', 'r')
        data = []
        for line in f:
            data.append(line.strip())
        data = data[0].split()
        self.chromosome = [int(data[0][1:])] + [int(e) for e in data[1:-1]] + [int(data[-1][:-1])] 

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
    
    def graphToGenome(self, edges):
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
        
        def unionEdges(edge):
            union(edge[0], edge[1])
            if 1 == edge[0] % 2:
                union(edge[0], edge[0]+1)
            else:
                union(edge[0], edge[0]-1)

            if 1 == edge[1] % 2:
                union(edge[1], edge[1]+1)
            else:
                union(edge[1], edge[1]-1)

        for e in edges:
            unionEdges(e)

        nodesID = dict()
        nodesSets = set()

        for e in edges:
            id = findParent(e[0])
            nodesID[e[0]] = id
            nodesID[e[1]] = id
            nodesSets.add(id)
        
        nodes = dict()
        sortedEdges = sorted(edges)
        for e in sortedEdges:
            id = nodesID[e[0]]
            if not id in nodes:
                nodes[id] = []
            nodes[id].append(e[0])
            nodes[id].append(e[1])
        
        genome = dict()
        for id, n in nodes.items():
            genome[id] = self.cycleToChromosome([n[-1]]+n[:-1])       
        
        genome = sorted(genome.values(), key = lambda x:abs(x[0]))
        return genome

    def printGenome(self, genome):
        result = ''
        for chromosome in genome:
            result += '('+' '.join(['+'+str(e) if e>0 else str(e) for e in chromosome])+')'
        print(result)
        
if __name__ == "__main__":
    ColoredEdges()