# python3

import sys
import numpy as np
import copy

'''
2-Break Sorting Problem: Find a shortest transformation of one genome into another by 2-breaks.
     Input: Two genomes with circular chromosomes on the same set of synteny blocks.
     Output: The sequence of genomes resulting from applying a shortest sequence of 2-breaks
     transforming one genome into the other.

Sample Input:
(+1 -2 -3 +4)
(+1 +2 -4 -3)
Sample Output:
(+1 -2 -3 +4)
(+1 -2 -3)(+4)
(+1 -2 -4 -3)
(-3 +1 +2 -4)

ShortestRearrangementScenario(P, Q)
     output P
     RedEdges ← ColoredEdges(P)
     BlueEdges ← ColoredEdges(Q)
     BreakpointGraph ← the graph formed by RedEdges and BlueEdges
     while BreakpointGraph has a non-trivial cycle Cycle
          (j, i′) ← an arbitrary edge from BlueEdges in a nontrivial red-blue cycle
          (i, j) ← an edge from RedEdges incident to node j
          (i', j') ← an edge from RedEdges incident to node i'
          RedEdges ← RedEdges with edges (i, j) and (i′, j′) removed
          RedEdges ← RedEdges with edges (j, i′) and (j′, i) added
          BreakpointGraph ← the graph formed by RedEdges and BlueEdges
          P ← 2-BreakOnGenome(P, i, j, i', j′)
          output P
'''

class TwoBreakSorting:
    def __init__(self):
        genomes = self._inputGenomes()
        result = self.shortestRearrangement(genomes[0], genomes[1])
        self.saveResult(result)

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

    def saveResult(self, result):
        f = open('result.txt', 'w')
        for r in result:
            d = self.printGenome(r)
            f.write(d+'\n')

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
        return result

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
    
    def edgeFromNontrivialCycle(self, edges, redEdges, blueEdges, blocks):
        # Output whether BreakpointGraph has a non-trivial cycle Cycle
        # Get an arbitrary edge from BlueEdges in a nontrivial red-blue cycle
        # Output removedRedEdges
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

        nodesID = dict()
        nodesSets = set()

        for e in edges:
            id = findParent(e[0])
            nodesID[e[0]] = id
            nodesID[e[1]] = id
            nodesSets.add(id)
        
        cycles = len(nodesSets)
        hasNontrivialCycle = False
        edge = None
        removedRedEdges = []
        if cycles != blocks:
            hasNontrivialCycle = True
            edgeDict = dict()
            redEdgeDict = dict()
            for e in edges:
                id = nodesID[e[0]]
                if not id in edgeDict:
                    edgeDict[id] = dict()
                edgeDict[id][e[0]] = e[1]
                edgeDict[id][e[1]] = e[0]
                if edge == None and len(edgeDict[id]) > 2 and e in blueEdges:
                    edge = (e[0], e[1])
                    edgeID = id
                if e in redEdges:
                    if not id in redEdgeDict:
                        redEdgeDict[id] = dict()
                    redEdgeDict[id][e[0]] = e[1]
                    redEdgeDict[id][e[1]] = e[0]
            removedRedEdges.append((edge[0], redEdgeDict[edgeID][edge[0]]))
            removedRedEdges.append((edge[1], redEdgeDict[edgeID][edge[1]]))
        return hasNontrivialCycle, removedRedEdges        

    def shortestRearrangement(self, P, Q):
        blocks = sum([len(a) for a in P])
        result = [P]
        redEdges = self.coloredEdges(P)
        blueEdges = self.coloredEdges(Q)
        breakpointGraph = redEdges.union(blueEdges)
        hasNontrivialCycle, removedRedEdges = self.edgeFromNontrivialCycle(breakpointGraph, redEdges, blueEdges, blocks)
        while hasNontrivialCycle:
            redEdges = self.twoBreakOnGraph(redEdges, removedRedEdges[0][0], removedRedEdges[0][1], removedRedEdges[1][0], removedRedEdges[1][1])
            breakpointGraph = redEdges.union(blueEdges)
            P = self.twoBreakOnGenome(P, removedRedEdges[0][0], removedRedEdges[0][1], removedRedEdges[1][0], removedRedEdges[1][1])
            hasNontrivialCycle, removedRedEdges = self.edgeFromNontrivialCycle(breakpointGraph, redEdges, blueEdges, blocks)
            result.append(P)
        return result

if __name__ == "__main__":
    TwoBreakSorting()