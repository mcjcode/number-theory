#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from utilities import isprime, factorize2, sqrtInt
from itertools import chain


def primesum_naive(n):
    """
    A slow, reference implementation.  Add up
    the all of the primes in the range [0, n]
    """
    retval = 0
    k = 2
    while k <= n:
        if isprime(k):
            retval += k
        k += 1
    return retval


def sieve(n):
    """
    Return the list of primes in [2..n]
    """
    S = np.ones((n+1, ), dtype=bool)
    S[0] = False
    S[1] = False
    p = 2
    sqrtn = sqrtInt(n)
    while p <= sqrtn:
        if S[p]:
            S[2*p::p] = False
        p += 1
    return filter(lambda pp: S[pp], range(n + 1))


def primesum_naive2(n):
    return sum(sieve(n))


sievesum_hash = {}


def sievesum1(n, p):
    """
    Return the sum of all of the numbers not
    excluded by a sieve of [1..n] by all primes
    up to and including p.
    
    Uses the lucy_hedgehog algorithm for summing
    the primes in [0, n] (faster than a sieve)

    In theory this works.  But for large n we
    hit maximum stack depth.  And we are using
    prime() to check primality when we should
    be able to determine it from our own
    computation.
    """

    if n == 0:
        retval = 0
    elif p == 0:
        retval = n*(n+1)//2 - 1
    elif (n, p) in sievesum_hash:
        return sievesum_hash[(n, p)]
    else:
        while not isprime(p):
            p -= 1
        retval = sievesum(n, p-1)
        if isprime(p) and p*p <= n:
            num_excluded = sievesum(n//p, p-1) - sievesum(p-1, p-1)
            retval -= p*num_excluded

    sievesum_hash[(n, p)] = retval
    return retval


def sievesum(n, p=None):
    """
    Return the sum of all numbers not excluded
    by a sieve of [1..n] by all primes up to
    and including p.
    """
    if p is None:
        p = sqrtInt(n)
    V = [n//i for i in range(1, p+1)]
    V += list(range(V[-1]-1, 0, -1))
    S = {i: i*(i+1)//2-1 for i in V}
    for p in range(2, p+1):
        if S[p] > S[p-1]:  # p is prime
            sp = S[p-1]  # sum of primes smaller than p
            p2 = p*p
            for v in V:
                if v < p2:
                    break
                S[v] -= p*(S[v//p] - sp)
    return S[n]


def sievecnt(n, p=None):
    """
    Return the count of numbers not excluded
    by a sieve of [1..n] by all primes up to
    and including p.  If p is not given then
    sieve up to sqrt(n). 
    """
    if p is None:
        p = sqrtInt(n)

    V = [n//i for i in range(1, p+1)]
    V += list(range(V[-1]-1, 0, -1))
    S = {i: i-1 for i in V}
    for p in range(2, p+1):
        if S[p] > S[p-1]:  # p is prime
            sp = S[p-1]  # number of primes smaller than p
            p2 = p*p
            for v in V:
                if v < p2:
                    break
                S[v] -= (S[v//p] - sp)
    return S[n]


def sievecntsum(n, p=None):
    """
    Return both the count and the sum of the
    numbers not excluded by a sieve of [1..n]
    by all primes up to and including p.
    """
    if p is None:
        p = sqrtInt(n)
    V = [n//i for i in range(1, p+1)]
    V += list(range(V[-1]-1, 0, -1))
    S0 = {i: i-1 for i in V}
    S1 = {i: i*(i+1)//2-1 for i in V}
    SP = {i: i*(i+1)//2-1 for i in V}

    for p in chain([2], range(3, p+1, 2)):
        if S0[p] > S0[p-1]:  # p is prime
            p2 = p*p
            for v in V:
                if v < p2:
                    break
                vmodp = v//p
                D0 = S0[vmodp] - S0[p-1]
                D1 = S1[vmodp] - S1[p-1]
                S0[v] -= D0
                S1[v] -= p*D1
                SP[v] -= p*(D1-D0)
    return S0[n], S1[n], SP[n]


def smpfsum(n):
    """
    Reference implementation for summing
    the smallest prime factor of all integers
    in range [2..n]
    """
    retval = 0
    for k in range(2, n+1):
        smpf = next(factorize2(k))[0]
        retval += smpf
    return retval
