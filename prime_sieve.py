#!/usr/bin/env python
#
# prime_sieve.py
#

from collections import Counter
import math
import numpy as np
from utilities import sqrtInt



def zeroit(arr, starti, stepi):
    arr[starti::stepi] = 0


def segmented_sieve(n, trace=False):
    """
    Return the list of all prime numbers
    less than or equal to n.  Runtime
    complexity of O(n) (or a bit better)
    with space complexity of O(sqrt(n)).
    """
    if n<=2:
        if n==2:
            yield 2
        return

    seglen = int((n+1)**(2/3))
    #seglen = sqrtInt(n+1)

    # always make seglen a multiple of 6
    # add 6 so that the first segment contains sqrt(n)
    seglen = max(6, 6*(seglen//6))+6

    wheel = [4, 2]
    wheel_len = len(wheel)
    wheelis = [1, 0]
    
    nsegs = (n+1)//seglen
    if seglen*nsegs < n+1:
        nsegs += 1

    ps = [2, 3]
    for p in ps:
        yield p

    seg_start = 0
    seg = np.ones((seglen, ), dtype=np.int32)
    seg[:4] = [0, 0, 1, 1]

    seg_end = seg_start+seglen-1

    ub = math.sqrt(seg_end)+1
    for p in ps:
        if p > ub:
            break
        seg[(-seg_start)%p::p] = 0

    starti = 1
    endi = min(seglen, ub-seg_start+1)
    stepi = 2

    i = starti
    wheeli = 0
    while i < endi: 
        if seg[i]:
            yield (p:=seg_start+i)
            ps.append(p)
            seg[p*p-seg_start::p] = 0
        i += wheel[wheeli]  # stepi
        wheeli = wheelis[wheeli]

    maxi = min(seglen, n-seg_start+1)
    while i < maxi:
        if seg[i]:  # seq_start+i is prime
            yield (p:=seg_start+i)
            ps.append(p) 
        i += wheel[wheeli]
        wheeli = wheelis[wheeli]

    seg_start += seglen
    seg[:] = 1

    for segi in range(1, nsegs):
        seg_end = seg_start+seglen-1
        ub = math.sqrt(seg_end)
        for p in ps:
            if p > ub:
                break
            seg[(-seg_start)%p::p] = 0
        starti = 1
        endi = min(seglen, ub-seg_start+1)
        stepi = 2
        i = starti
        wheeli = 0
        while i < endi: 
            if seg[i]:
                yield (p:=seg_start+i)
                seg[p*p-seg_start::p] = 0
            i += wheel[wheeli]  # stepi
            wheeli = wheelis[wheeli]
        maxi = min(seglen, n-seg_start+1)
        while i < maxi:
            if seg[i]:  # seq_start+i is prime
                yield (p:=seg_start+i)
            i += wheel[wheeli]
            wheeli = wheelis[wheeli]
        seg_start += seglen
        seg[:] = 1


def indicator(n):
    nm1 = n-1
    while True:
        yield True
        i = 0
        while i < nm1:
            yield False
            i += 1


def sqfree_parts():
    """
    For each positive integer n,
    yield the list of primes dividing
    n, w/o multiplicity.  Horribly
    inefficient.  Just a proof of
    concept.
    """
    yield []  # 1 has no prime divisors
    indicators = []
    i = 2
    while True:
        ps = []
        for (p, pi) in indicators:
            if next(pi):
                ps.append(p)
        if not ps:
            ps = [i]
            g = indicator(i) 
            _ = next(g)
            indicators.append((i, g))
        yield ps
        i += 1


def primeFactors(n):
    """
    Return a length n+1 list pf such that pf[k]
    is the list of primes dividing k
    """
    retval = [[] for _ in range(n+1)]
    for p in range(2, n+1):
        if not retval[p]:  # p is prime
            for j in range(p, n+1, p):
                retval[j].append(p)
    return retval

    
def spf(n: int):
    """
    Return the smallest prime factors of all integers
    up to and including n
    """
    retval = [0]*(n+1)
    for p in range(2, n+1):
        if not retval[p]: # p is prime
            for j in range(p, n+1, p):
                if not retval[j]:
                    retval[j]=p
    return retval


def lpf(n: int):
    """
    Return the largest prime factor of all integers
    up to and including n
    """
    retval = np.array([0]*(n+1), dtype=object)
    for p in range(2, n+1):
        if not retval[p]: # p is prime
            retval[p:n+1:p] = p
    return retval


def factorizations(N: int):
    """
    Returns a length N+1 list of Counters, such that the
    element in position n is the factorization of n.
    """
    lgpf = lpf(N)
    retval = [None]*(N+1)
    retval[1] = Counter()
    for n in range(2, N+1):
        p = lgpf[n]
        fact = retval[n//p].copy()
        fact[p] += 1
        retval[n] = fact
    return retval


def prime_factors(N: int):
    lgpf = lpf(N)
    retval = [None]*(N+1)
    retval[1] = set()
    for n in range(2, N+1):
        p = lgpf[n]
        ps = retval[n//p].copy()
        ps.add(p)
        retval[n] = ps
    return retval


def phiarray(n: int):
    """
    Return [_, phi(1), phi(2),..., phi(n)]
    """
    pass
