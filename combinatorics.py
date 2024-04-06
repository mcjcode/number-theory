#!/usr/bin/env python

import numpy as np
import itertools
from collections import defaultdict
from polynomial import intpoly1d # for chromatic_polynomial
from utilities import powerset

def comb(k, a, b):
    """
    Yield the sequence of all k element subsets of the
    consecutive list of integers [a,..,b]
    """
    if k <= b-a+1:  # there are enough left to choose k
        if k == 0:
            yield []
        else:
            for tl in comb(k-1, a+1, b):
                yield [a] + tl
            for tl in comb(k, a+1, b):
                yield tl


def nchoosek(n, ks):
    """
    n: a non-negative integer
    ks: either an integer in the range 0..n,
        or a list of non-negative integers
        k1,...,km whose sum is n
    returns: The multinomial coefficient
        (n choose k1 . . . km)
    """
    if type(ks) == int:
        if ks < 0:
            return 0
        ks = [ks, n-ks]

    retval = 1
    i = n
    while i > 0:
        for k in ks:
            j = 1
            while j <= k:
                retval *= i
                i -= 1
                retval //= j
                j += 1
    return retval


def stirling2(n, k):
    """
    n: a non-negative integer
    k: a non-negative integer
    returns: the number of equivalence relations of
        the numbers [1,...,n] with k classes.

    We use count surjections from [1,..,n] to [1,..,k] by inclusion
    exclusion, then divide by k!
    """
    retval = 0
    nchk = 1
    for j in range(k+1):
        retval += (-1)**(j%2)*nchk*(k-j)**n
        nchk = nchk * (k-j) // (j+1)
    for j in range(1, k+1):
        retval //= j
    return retval


def eqrels(N):
    """
    Return all possible equivalence relations on [0,n)
    """
    h = [[]] # the equivalence relations on [0,1)
    for n in range(N):
        h2 = []
        for r in h:
            for j1 in range(len(r)):
                h2.append([ss+([n]*(j2==j1)) for j2, ss in enumerate(r)])
            h2.append(r + [[n]])
        h = h2
    return h


def monomials(d: int, n: int):
    """
    Yield the exponent vectors of all degree
    d monomials in n variables
    """
    for vv in comb(n-1, 1, d+n-1):
        yield [vv[0]-1] + [vv[i+1]-vv[i]-1 for i in range(n-2)] + [d+n-1-vv[-1]]


def conjugate_partition(p):
    """
    Return the conjugate partition to p
    """
    return [len([pp for pp in p if pp >= k]) for k in range(1, max(p)+1)]


def partitions(d, max_part=None):
    """
    Return a generator yielding the partitions of d
    """
    if max_part is None:
        max_part = d

    if d == 0:
        yield []
    else:
        for pp in range(max_part, 0, -1):
            for i in range(1, d//pp+1):  # number of parts of size pp
                seg = [pp]*i
                for part in partitions(d-i*pp,pp-1):
                    yield seg + part


def peterson_graph():
    """
    Return the adjacency matrix of the
    Peterson graph
    """
    A = np.identity(n=10,dtype=int)
    for i in range(5):
        A[i,(i+1)%5] = 1
        A[5+i,5+((i+2)%5)] = 1
        
        A[(i+1)%5,i] = 1
        A[5+((i+2)%5),5+i] = 1
        
        A[i,5+i] = 1
        A[5+i,i] = 1
    return A


def connected_components(G):
    """
    Input: G - a graph
    Output: the connected_components of G
    """
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
            components.append(C)
    return components


def nconnected_components(G):
    """
    Input: G - a graph
    Output: the # of connected_components of G
    """
    return len(connected_components(G))


def gray(n):
    """
    The nth term in the reflective
    gray code.
    """
    n -= 1
    return n^(n>>1)


def graybits():
    """
    a generator yielding the sequence
    of bit locations needed to flip
    to product the reflective gray code.
    """
    i = 1
    while True:
        a = gray(i)^gray(i+1)
        retval = 0
        while a:
            a //= 2
            retval += 1
        yield retval-1
        i += 1


def chromatic_polynomial2(G):
    vertices = list(G)
    n = len(G)
    edges = set()
    for i in G:
        for j in G[i]:
            if (i, j) not in edges and (j, i) not in edges:
                edges.add((i,j))
    edges = list(edges)
    ne = len(edges)
    retval = intpoly1d([0]*n + [(-1)**(ne)])
    it = powerset(edges); next(it)
    for ss in it:
        na = len(ss)
        G = {v:set() for v in vertices}
        for (i,j) in ss:
            G[i].add(j)
            G[j].add(i)
        cA = nconnected_components(G)
        retval = retval + (-1)**((ne-na)%2) * intpoly1d([0]*cA+[1])
    return (retval.coefs[-1])*retval        


def chromatic_polynomial(G):
    """
    Return the chromatic polynomial p of
    the graph G.

    Algorithm is O(2**(#edges(G))), and spends
    almost all of its time computing the number
    of connected components of the subgraphs
    of G.
    """
    vertices = list(G)
    edges = set()
    for i in G:
        for j in G[i]:
            if (i, j) not in edges and (j, i) not in edges:
                edges.add((i,j))
    edges = list(edges)
    sgn = 1
    retval = [0]*len(G) + [sgn]

    # walk through every subgraph of G, using
    # a graycode to avoid having to construct
    # every subgraph from scratch
    g = graybits()
    G2 = {v:set() for v in vertices}
    for _ in range(2**len(edges)-1):
        i, j = edges[next(g)]
        if j in G2[i]:
            G2[i].remove(j)
            G2[j].remove(i)
        else:
            G2[i].add(j)
            G2[j].add(i)
        sgn = -sgn
        cA = nconnected_components(G2)
        retval[cA] += sgn
        
    return intpoly1d(retval)

def newtons_identity(ps):
    """
    :param ps: the values of the power symmetric polynomials
    :return: the corresponding values of the elementary symmetric polynomials

    Example:

    x1, x2, x3 = 1, 2, 3
    ps = [ x1**n + x2**n + x3**n for n in range(3)]
    es = newtons_identities(ps)
    assert es == [1, 6, 11, 6]
    
    """
    es = [1]
    for j in range(1, len(ps)):
        z = 0
        sign = 1
        for i in range(1, j+1):
            z += sign*ps[i]*es[j-i]
            sign *= -1
        es.append(z//j)
    return es
