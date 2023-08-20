import numpy as np
from heaps import heapset

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
    M = np.zeros((n, n),dtype=np.float)
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
    M = np.zeros((n,n),dtype=np.float)
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
