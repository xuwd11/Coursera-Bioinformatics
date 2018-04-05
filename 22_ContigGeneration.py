# python3

import sys

'''
Solve the Contig Generation Problem.
     Input: A collection of k-mers Patterns. 
     Output: All contigs in DeBruijn(Patterns).

Sample Input:
ATG
ATG
TGT
TGG
CAT
GGA
GAT
AGA

Sample Output:
AGA ATG ATG CAT GAT TGGA TGT
'''
'''
Implement MaximalNonBranchingPaths.
     Input: The adjacency list of a graph whose nodes are integers.
     Output: The collection of all maximal nonbranching paths in this graph.

Sample Input:
1 -> 2
2 -> 3
3 -> 4,5
6 -> 7
7 -> 6

Sample Output:
1 -> 2 -> 3
3 -> 4
3 -> 5
7 -> 6 -> 7
'''
'''
A node v in a directed graph Graph is called a 1-in-1-out node if its indegree and outdegree are both equal to 1, 
i.e., in(v) = out(v) = 1.  We can rephrase the definition of a "maximal non-branching path" from the main text as
a path whose internal nodes are 1-in-1-out nodes and whose initial and final nodes are not 1-in-1-out nodes. Also, 
note that the definition from the main text does not handle the special case when Graph has a connected component 
that is an isolated cycle, in which all nodes are 1-in-1-out nodes.

The MaximalNonBranchingPaths pseudocode below generates all non-branching paths in a graph. It iterates through all
nodes of the graph that are not 1-in-1-out nodes and generates all non-branching paths starting at each such node. 
In a final step, MaximalNonBranchingPaths finds all isolated cycles in the graph.

    MaximalNonBranchingPaths(Graph)
        Paths ← empty list
        for each node v in Graph
            if v is not a 1-in-1-out node
                if out(v) > 0
                    for each outgoing edge (v, w) from v
                        NonBranchingPath ← the path consisting of the single edge (v, w)
                        while w is a 1-in-1-out node
                            extend NonBranchingPath by the outgoing edge (w, u) from w 
                            w ← u
                        add NonBranchingPath to the set Paths
        for each isolated cycle Cycle in Graph
            add Cycle to Paths
        return Paths
'''

class MaximalNonBrachingPath:
    def __init__(self, adj):
        self.adj = adj
        self.updateAdj(self.adj)
        self.paths = self.findMaximalNonBranchingPaths()
    
    def _input(self):
        data = list(sys.stdin.read().strip().split())
        adj = dict()
        for i in range(len(data) // 3):
            nl = data[i*3]
            nrList = data[i*3+2].split(',')
            adj[nl] = adj.get(nl, []) + nrList
            for nr in nrList:
                if nr not in adj:
                    adj[nr] = []
        return adj
    
    def updateAdj(self, adj):
        self.n = len(adj)
        self.inDeg = dict()
        self.outDeg = dict()
        for w, vList in adj.items():
            self.inDeg[w] = self.inDeg.get(w, 0)
            for v in vList:
                self.inDeg[v] = self.inDeg.get(v, 0) + 1
            self.outDeg[w] = len(vList)
        return
    
    def findMaximalNonBranchingPaths(self):
        paths = []
        nodes1in1out = set() #1-in-1-out nodes
        nExplored = set() #1-in-1-out nodes which were explored
        for v in self.adj.keys():
            if not (1 == self.inDeg[v] and 1 == self.outDeg[v]):
                if self.outDeg[v] > 0:
                    for w in self.adj[v]:
                        nbPath = [v, w] #NonBrachingPath
                        while 1 == self.inDeg[w] and 1 == self.outDeg[w]:
                            nExplored.add(w)
                            u = self.adj[w][0]
                            nbPath.append(u)
                            w = u
                        paths.append(nbPath)
            else:
                nodes1in1out.add(v)            
        for v in nodes1in1out:
            if v not in nExplored:
                w = v
                nbPath = []
                while w in nodes1in1out:
                    nbPath.append(w)
                    if w == v and len(nbPath) > 1:
                        paths.append(nbPath)
                        for node in nbPath:
                            nExplored.add(node)
                        break
                    w = self.adj[w][0]
        return paths
    
    def printPaths(self, paths):
        for p in paths:
            print(' -> '.join(p))

class ContigGeneration:
    def __init__(self):
        self.adj = self.readData()
        self.paths = MaximalNonBrachingPath(self.adj).findMaximalNonBranchingPaths()
        self.contigs = self.reconstructFromPaths(self.paths)
        print(' '.join(self.contigs))
    
    def readData(self):
        data = list(sys.stdin.read().strip().split())
        adj = self.DeBrujin(data) 
        return adj

    def DeBrujin(self, patterns):
        k = len(patterns[0])
        adjdb = dict()
        for p in patterns:
            pl = p[:k-1]
            pr = p[1:]
            if pl in adjdb:
                adjdb[pl].append(pr)
            else:
                adjdb[pl] = []
                adjdb[pl].append(pr)
            if pr not in adjdb:
                adjdb[pr] = []
        return adjdb
    
    def reconstructFromPaths(self, paths):
        contigs = []
        for p in paths:
            contigs.append(p[0] + ''.join([seq[-1] for seq in p[1:]]))
        return contigs

if __name__ == "__main__":
    ContigGeneration()