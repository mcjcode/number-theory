import itertools
from collections import defaultdict, OrderedDict
import numpy as np
from heaps import heapset

import pynauty

from utilities import argmin


def dijkstra(start, end, costs, neighbors):
    visited = set()
    h = heapset()
    h.heappush(costs(start), start)
    while h:
        cost, node = h.heappop()
        if node==end:
            return cost
        for nbr in neighbors(node):
            if nbr not in visited:
                this_cost = cost + costs(nbr)
                if nbr in h:
                    if this_cost < h[nbr]:
                        h.raise_priority(this_cost, nbr)
                else:
                    h.heappush(this_cost, nbr)
        visited.add(node)

def dijkstra2(G, start, end):
    visited = set()
    h = heapset()
    h.heappush(0, start)
    while h:
        cost, node = h.heappop()
        if node==end:
            return cost
        for nbr, dist in G[node]:
            if nbr not in visited:
                this_cost = cost + dist
                if nbr in h:
                    if this_cost < h[nbr]:
                        h.raise_priority(this_cost, nbr)
                else:
                    h.heappush(this_cost, nbr)
        visited.add(node)

def exit_probabilities(G, end, before):
    """
    Inputs

    G: a graph (a dict of node:neighbors)
    end: a node of graph G
    before: a set or list of nodes of G

    Output:

    A dictionary with keys equal to the nodes of G, and values equal
    to the probability that a (uniform) random walk on G will arrive
    at 'end' before first arriving at one of the nodes in 'before'.
    """
    I = np.identity(len(G))
    M = np.array([[ 1/len(G[v])*(w in G[v])*(v!=end)*(v not in before) for w in G] for v in G])
    b = np.array([v==end for v in G])
    return dict(zip(G,np.linalg.solve(I-M,b)))


def exit_probabilities_weighted(G, end, before):
    """
    Inputs

    G: a graph (a dict of node:(neighbors, probability) pairs)
    end: a node of graph G
    before: a set or list of nodes of G

    Output:

    A dictionary with keys equal to the nodes of G, and values equal
    to the probability that a (weighted) random walk on G will arrive
    at 'end' before first arriving at one of the nodes in 'before'.
    """
    n = len(G)
    I = np.identity(n)
    indices = {v:i for (i,v) in enumerate(G)}
    M = np.zeros((n, n),dtype=float)
    for v, i in indices.items():
        if v!=end and v not in before:
            for w, p in G[v]:
                j = indices[w]
                M[i,j] = p
    b = np.array([v==end for v in G])
    return dict(zip(G,np.linalg.solve(I-M,b)))


def expected_exit_time(G, end):
    """
    Inputs

    G: a graph (a dict of node:neighbors)
    end: a node of graph G

    Output:

    A dictionary with keys equal to the nodes of G, and values equal
    to the expected number of steps that a (uniform) random walk on G
    will take to arrive at 'end'.
    """

    I = np.identity(len(G))
    M = np.array([[1/len(G[v])*(w in G[v])*(v!=end) for w in G] for v in G])
    b = np.array([v!=end for v in G])
    return dict(zip(G,np.linalg.solve(I-M, b)))


def conditional_expected_exit_time(G, end, other):
    """
    Return the expected time that a random walk on G
    will take to arrive at node 'end' without first
    landing at one of the nodes of 'other'

    We implement the approach outlined here:

    https://math.stackexchange.com/questions/2852105/expected-time-till-absorption-in-specific-state-of-a-markov-chain

    which essentially entails defining a conditioned transition
    matrix that is a Bayesean update of the original one.
    """
    # first compute the probabilities that a random
    # walk
    n = len(G)
    probs = exit_probabilities(G, end, other)
    I = np.identity(n)
    M = np.zeros((n,n),dtype=float)
    indices = {v:i for (i,v) in enumerate(G)}
    for i, v in enumerate(G):
        if v==end or v in other:
            M[i,i] = 0.0
        else:
            for w in G:
                if w in G[v]:
                    j = indices[w]
                    M[i,j] = probs[w]/probs[v]*(1/len(G[v]))
    b = np.array([v!=end for v in G])
    return dict(zip(G, np.linalg.solve(I-M, b)))


def isconnected(G):
    """
    Input: G - a graph
    Output: is the graph connected?
    """
    n = len(G)

    V = defaultdict(lambda:True)
    components = []
    for v in G:
        if V[v]:
            component = []
            stack = [v]
            while stack:
                v = stack.pop()
                component.append(v)
                V[v] = False
                for w in G[v]:
                    if V[w]:
                        stack.append(w)
            C = {v:G[v] for v in component}
            return len(C)==n
    return True # the empty graph is connected



#
# Code related to finding the earliest (lexicographic)
# labelling of a simple graph.
#

def delete_node(G, node):
    return {n:[nbr for nbr in nbrs if nbr!=node] for (n, nbrs) in G.items() if n!=node}


def relabel(Gt, perm):
    return tuple(
        sorted(
            tuple(
                (perm[nd], tuple(sorted(tuple(perm[nbr] for nbr in nbrs))))
                for nd, nbrs in Gt))
    )


def normalizer(Gt, node):
    """
    Return the first (lexicographically) graph labeling
    of the graph 'Gt' with the distinguished base node 'node'
    """

    n = len(Gt)
    nds = tuple(nd for nd, _ in Gt)
    perm = {nd:i for i, nd in enumerate(nds)}
    Gt = relabel(Gt, perm)
    node = perm[node]
    f = lambda p: (relabel(Gt, {i:j for i,j in enumerate(p)}), p[node])
    p = argmin(f, itertools.permutations(range(n)))
    return {nd:p[perm[nd]] for nd in nds}


def normalizer2(Gt, node):
    n = len(Gt)
    #
    # make sure that our graph starts out
    # as a graph on {0, . . ., n-1}
    #
    nds = tuple(nd for nd, _ in Gt)
    perm = {nd:i for i, nd in enumerate(nds)}
    Gt = relabel(Gt, perm); node = perm[node]
    g = pynauty.Graph(n)
    for nd, nbrs in Gt:
        g.connect_vertex(nd, list(nbrs))
    g.set_vertex_coloring([{node}])

    labels = pynauty.canon_label(g)
    return {nd:labels[perm[nd]] for nd in nds}


def nhamiltonians0(G, start):
    def tuplify(G):
        return tuple((nd,tuple(nbrs)) for nd, nbrs in G.items())

    def dictify(tG):
        return {nd:list(nbrs) for (nd, nbrs) in tG}

    h = {(tuplify(G), start):1}
    for j in range(len(G)-1):
        h2 = defaultdict(int)
        for (g, last), count in h.items():
            g = dictify(g)
            for nd in g[last]:
                gm = delete_node(g, last)
                if nconnected_components(gm)==1:
                    gt = tuplify(gm)
                    p = normalizer2(gt, nd)
                    gt2 = relabel(gt, p)
                    nd2 = p[nd]
                    h2[gt2, nd2] += count
        h = h2
        if j>=len(G)-4:
            for k, v in h.items():
                print(k, v)
        print(j+1, len(h))
    return sum(h.values())


def nhamiltonians(g: pynauty.Graph,
                  start: int):
    """
    Count the number of hamiltonians of the graph 'g'
    starting at the node 'start'.
    """
    g.set_vertex_coloring([{start}])
    h = {pynauty.certificate(g):(g, start, 1)}
    for j in range(len(g.adjacency_dict)-1):
        h2 = {}
        for cert, (g, last, count) in h.items():
            n = len(g.adjacency_dict)
            for nd in g.adjacency_dict[last]:

                new_adj_dict = {nd:(set(g.adjacency_dict[nd]) - {last}) for nd in g.adjacency_dict if nd!=last}
                perm = {nd:i for i, nd in enumerate(new_adj_dict.keys())}
                new_adj_dict = {perm[nd]:list(map(perm.__getitem__, new_adj_dict[nd])) for nd in new_adj_dict}
                nd = perm[nd]

                new_g = pynauty.Graph(number_of_vertices=len(new_adj_dict),
                                      adjacency_dict=new_adj_dict,
                                      vertex_coloring=[{nd}])
                cert = pynauty.certificate(new_g)
                if cert in h2:
                    h2[cert][2] += count
                else:
                    h2[cert] = [new_g, nd, count]
        h = h2
        # uncomment this line if you want to see how many
        # states there are at each stage.
        # print(j+1, len(h))
    return sum(count for _, _, count in h.values())
