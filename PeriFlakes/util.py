import cornflakes as cf
import numpy as np
from collections import defaultdict


def Make_Bond_Adjacency(H):
    edges = H.view()[0]
    newgraph = defaultdict(lambda : [])
    for l in edges:
        newgraph[l[0]].append( (l[1],l[2] ) )
        newgraph[l[1]].append( (l[0],l[2] ) )
    HA = cf.Hypergraph()
    for center,bonds in newgraph.iteritems():
        egg = np.empty( 2*len(bonds)+1 , dtype=np.intc )
        egg[0] = center
        for i,b in enumerate(bonds):
            egg[1+i] = b[0]
            egg[1+len(bonds)+i] = b[1]
        HA.Push_Edge(egg)
    return HA
