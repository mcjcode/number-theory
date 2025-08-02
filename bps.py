#!/usr/bin/env python
#
"""
Routines for looping over numbers formed from a
given sequence of primes below a fixed bound.
"""

from utilities import sqrtInt, prod


def bps(n, xs):
    """
    For an increasing sequence xs of relatively
    prime numbers and an upper bound n>=1, compute
    all of the products <=n of subsets of xs.
    
    This function proceeds depth first in the tree and
    is non-recursive (it does not call it self and is not
    limited by python's recursion depth bound)
    """

    s = [(1, 0)]
    lenxs = len(xs)
    while s:
        prd, idx = s.pop()
        if idx < lenxs:
            x = xs[idx]
            if x*prd <= n:                
                s.append((prd,   idx+1))
                s.append((prd*x, idx+1))
                continue
        yield prd


def bps_w_rep(n, ps, i=0, only_powerful=False):
    """
    Return all numbers in [1..n] whose prime
    factors are all in the list ps[i:] (repetitions
    allowed)
    """
    yield 1
    for pi in range(i, len(ps)): # pi indexes the first prime
        p = ps[pi]
        pk = p*p if only_powerful else p
        while pk<=n:
            for x in bps_w_rep(n//pk, ps, pi+1, only_powerful=only_powerful):
                yield pk*x
            pk *= p


def bps_facts_w_rep(n, ps, i=0, only_powerful=False):
    yield []
    for pi in range(i, len(ps)): # pi indexes the first prime
        p = ps[pi]
        k = 1 if only_powerful else 0
        pk = p**k
        if pk*p>n:
            break
        while (pk:=pk*p)<=n:
            k += 1
            for x in bps_facts_w_rep(n//pk, ps, pi+1, only_powerful=only_powerful):
                yield [(p, k)] + x
    

def bps_w_sign(ps, n):
    """
    For an increasing sequence ps of relatively
    prime numbers ps and an upper bound n>=1, compute
    all of the products <=n of subsets of xs, along
    with a parity: +1 if an even number of elements
    are used, -1 if an odd number are used.
    
    This function proceeds depth first in the tree and
    is non-recursive (it does not call it self and is not
    limited by python's recursion depth bound)
    
    Also, have found this 20x faster than bps_w_sign
    """
    s = [(+1, 1, 0)]
    lenps = len(ps)
    while s:
        sgn, prd, idx = s.pop()
        if idx < lenps:
            p = ps[idx]
            if p*prd <= n:                
                s.append(( sgn, prd,   idx+1))
                s.append((-sgn, prd*p, idx+1))
                continue
        yield sgn, prd


def bpsk(ps, N, k, i=0):
    """
    Yield k-tuples of strictly increasing elements of the
    list increasing list ps[i:] whose product is less than
    or equal to N.
    """
    if k>len(ps)-i or N<1: # even the empty tuple has a product of '1'.
        return
    if k==0:
        yield (), 1
    else:
        pi = i
        p = ps[pi]
        while p**k <= N:
            for t, x in bpsk(ps, N//p, k-1, pi+1):
                yield (p,)+t, x*p
            pi += 1
            if pi==len(ps):
                break
            p = ps[pi]  


def bpsk_w_rep(ps, N, k, i=0):
    """
    Yield k-tuples of strictly increasing elements of the
    list increasing list ps[i:] whose product is less than
    or equal to N.
    """
    if k>len(ps)-i or N<1: # even the empty factorization has a product of '1'.
        return
    elif k==0:
        yield (), 1
    else:
        pi = i
        p = ps[pi]
        while p**k <= N:
            j = 1
            pj = p**j
            while pj <= N:
                for t, x in bpsk_w_rep(ps, N//pj, k-1, pi+1):
                    yield ((p, j),)+t, x*pj
                j += 1
                pj = pj*p
            pi += 1
            if pi==len(ps):
                break
            p = ps[pi]
