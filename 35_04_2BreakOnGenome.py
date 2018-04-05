# python3

import sys
import numpy as np
import copy

'''
1. Implement 2-BreakOnGenomeGraph.
     Input: The colored edges of a genome graph GenomeGraph, followed by indices i, i', j, and j'.
     Output: The colored edges of the genome graph resulting from applying the 2-break operation
     2-BreakOnGenomeGraph(GenomeGraph, i, i′, j, j′).

Sample Input:
(2, 4), (3, 8), (7, 5), (6, 1)
1, 6, 3, 8
Sample Output:
(2, 4), (3, 1), (7, 5), (6, 8)

2-BreakOnGenomeGraph(GenomeGraph, i, i′, j, j′)
     remove colored edges (i, i') and (j, j′) from GenomeGraph
     add colored edges (i, j) and (i′, j') to GenomeGraph
     return GenomeGraph

2. Implement 2-BreakOnGenome.
     Input: A genome P, followed by indices i, i', j, and j'.
     Output: The genome P' resulting from applying the 2-break operation 2-BreakOnGenome(GenomeGraph, i, i′, j, j′).

Sample Input:
(+1 -2 -4 +3)
1, 6, 3, 8
Sample Output:
(+1 -2)(-3 +4)

2-BreakOnGenome(P, i, i′, j, j′)
     GenomeGraph ← BlackEdges(P) and ColoredEdges(P)
     GenomeGraph ← 2-BreakOnGenome(GenomeGraph, i, i′, j, j′)
     P ← GraphToGenome(GenomeGraph)
     return P
'''

class TwoBreakOnGenome:
    def __init__(self):
        '''
        edges, i0, i1, j0, j1 = self._inputGraphAndChange()
        editedEdges = self.twoBreakOnGraph(edges, i0, i1, j0, j1)        
        print(str(editedEdges)[1:-1])
        '''
        genome, i0, i1, j0, j1 = self._inputGenomeAndChange()
        genome = self.twoBreakOnGenome(genome, i0, i1, j0, j1)
        self.printGenome(genome)
    
    def _inputGenomeAndChange(self):
        data = sys.stdin.read().strip().split('\n')
        d0 = data[0].split(')(')
        genome = []
        for d in d0:
            d = d.split()
            genome.append([int(d[0][1:] if '('==d[0][0] else d[0])] + [int(e) for e in d[1:-1]] +\
            [int(d[-1][:-1] if ')'==d[-1][-1] else d[-1])])
        i0, i1, j0, j1 = [int(d) for d in data[1].split(',')]
        return genome, i0, i1, j0, j1 

    def _inputGraphAndChange(self):
        data = sys.stdin.read().strip().split('\n')
        d0 = data[0].split('), (')
        edges = set()
        for d in d0:
            d = d.split(',')
            edges.add((int(d[0][1:] if '('==d[0][0] else d[0]), int(d[1][:-1] if ')'==d[1][-1] else d[1])))
        i0, i1, j0, j1 = [int(d) for d in data[1].split(',')]
        return edges, i0, i1, j0, j1

    def _inputGenomes(self):
        data = sys.stdin.read().strip().split('\n')
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
        
    def printGenome(self, genome):
        result = ''
        for chromosome in genome:
            result += '('+' '.join(['+'+str(e) if e>0 else str(e) for e in chromosome])+')'
        print(result)

    def twoBreakOnGraph(self, edges, i0, i1, j0, j1):
        edges.discard((i0, i1))
        edges.discard((i1, i0))
        edges.discard((j0, j1))
        edges.discard((j1, j0))
        edges.add((i0, j0))
        edges.add((i1, j1))
        return edges

    def groupNodes(self, edges):
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
        
        return nodesSets, nodesID
    
    def buildEdgeDict(self, edges, nodesSet, nodesID):
        edgeDict = dict()
        for e in edges:
            id = nodesID[e[0]]
            if not id in edgeDict:
                edgeDict[id] = dict()
            edgeDict[id][e[0]] = e[1]
            edgeDict[id][e[1]] = e[0]
        return edgeDict
            
    def twoBreakOnGenome(self, genome, i0, i1, j0, j1):
        edges = self.twoBreakOnGraph(self.coloredEdges(genome), i0, i1, j0, j1)
        nodesSet, nodesID = self.groupNodes(edges)
        edgeDict = self.buildEdgeDict(edges, nodesSet, nodesID)
        nodesDict = dict()
        for id, eDict in edgeDict.items():
            nodesDict[id] = []
            currNode0 = list(eDict)[0]
            while len(eDict) > 0:
                nodesDict[id].append(currNode0)
                if 1 == currNode0 % 2:
                    currNode1 = currNode0+1
                else:
                    currNode1 = currNode0-1
                nodesDict[id].append(currNode1)
                newNode = eDict[currNode1]
                del eDict[currNode0]
                del eDict[currNode1]
                currNode0 = newNode
        newGenome = dict()
        for id, nodes in nodesDict.items():
            newGenome[id] = self.cycleToChromosome(nodes)
        newGenome = sorted(newGenome.values(), key = lambda x:abs(x[0]))
        return newGenome

if __name__ == "__main__":
    TwoBreakOnGenome()