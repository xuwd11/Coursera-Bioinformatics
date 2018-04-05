# python3

import sys
import numpy as np
import copy

'''
1. Implement ChromosomeToCycle.
     Input: A chromosome Chromosome containing n synteny blocks.
     Output: The sequence Nodes of integers between 1 and 2n resulting from applying ChromosomeToCycle to Chromosome.

Sample Input:
(+1 -2 -3 +4)
Sample Output:
(1 2 4 3 6 5 7 8)

ChromosomeToCycle(Chromosome)
     for j ← 1 to |Chromosome|
          i ← Chromosomej
          if i > 0
               Nodes2j−1 ←2i−1
               Nodes2j ← 2i
          else
               Nodes2j−1 ← -2i
               Nodes2j ←-2i−1
     return Nodes

2. Implement CycleToChromosome.
     Input: A sequence Nodes of integers between 1 and 2n.
     Output: The chromosome Chromosome containing n synteny blocks resulting from applying CycleToChromosome to Nodes.

Sample Input:
(1 2 4 3 6 5 7 8)
Sample Output:
(+1 -2 -3 +4)

CycleToChromosome(Nodes)
     for j ← 1 to |Nodes|/2
          if Nodes2j−1 < Nodes2j
               Chromosomej ← Nodes2j /2
          else
               Chromosomej ← −Nodes2j−1/2
     return Chromosome
'''

class ChromosomeCycle:
    def __init__(self):
        #self._inputChromosome()
        #nodes = self.chromosomeToCycle(self.chromosome)
        #print('('+' '.join([str(n) for n in nodes])+')')
        self._inputNodes()
        chromosome = self.cycleToChromosome(self.nodes)
        self.printChromosome(chromosome)

    def _inputChromosome(self):
        data = sys.stdin.read().strip().split()
        self.chromosome = [int(data[0][1:])] + [int(e) for e in data[1:-1]] + [int(data[-1][:-1])]

    def _inputNodes(self):
        data = sys.stdin.read().strip().split()
        self.nodes = [int(data[0][1:])] + [int(e) for e in data[1:-1]] + [int(data[-1][:-1])]
    
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

if __name__ == "__main__":
    ChromosomeCycle()