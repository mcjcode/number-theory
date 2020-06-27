#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

import numpy as np

from utilities import (
    prod,
    factorize2,
    gcd,
    multiplicities,
    sqrtInt
)


def phi(nn):
    """
    :param nn: a positive integer
    :return: the number of natural numbers 1 <= kk < nn
             that are relatively prime to nn.
    """
    factors = factorize2(nn)
    return prod([pp**(kk-1)*(pp-1) for (pp, kk) in factors])

def divisor_function(kk, nn):
    """
    return the sum of d^k over all divisors of nn
    """
    factors = factorize2(nn)
    return prod([(ee+1) if kk==0 else (pp**(kk*(ee+1))-1)//(pp**kk-1) for (pp, ee) in factors])

def sumrange(a, b):
    return b*(b+1)//2 - a*(a+1)//2
    
def sum_sigma0(n):
    """
    return the sum of sigma0(k) for
    all k in [1..n]
    """
    sqrtk=sqrtInt(n)
    part1=sum(n//d for d in range(1, n//sqrtk+1))
    part2=sum(((n//d)-(n//(d+1)))*d for d in range(1, sqrtk))
    return part1+part2
     
def sum_sigma1(n):
    """
    return the sum of sigma1(k) for
    all k in [1..n]
    """
    sqrtk=sqrtInt(n)
    part1=sum(d*(n//d) for d in range(1, n//sqrtk+1))
    part2=sum(sumrange(n//(d+1), n//d)*d for d in range(1, sqrtk))
    return part1+part2

_maxp = 16
_sgns = [1]*(2**_maxp)
c=1
for i in range(_maxp):
    for j in range(c):
        _sgns[j+c] = -_sgns[j]
    c <<= 1
def partial_totient(n, k):
    """
    :param n: a positive integer 
    :param k: a positive integer
    :return: the number of integers in [1..n]
        that are relatively prime to k.  Note
        that the number of distinct prime factors
        of k must be <= _maxp (currently 16) 
    """

    if n==0: return 0
    ps = [p for (p, _) in factorize2(k)]
    nps = len(ps)

    assert nps <= _maxp

    pps = [1]*(2**nps)
    c=1
    for p in ps:
        j = 0
        while j<c:
        #for j in xrange(c): 
            pps[j+c] = p*pps[j]
            j += 1
        c <<= 1

    retval = 0
    ub = 2**nps
    i = 0
    while i < ub:
        retval += _sgns[i]*(n//pps[i])
        i += 1
    return retval

def _partial_totient_alternate(n, k):
    """
    Return the number of integers in [1..n]
    that are relatively prime to k
    """
    if n==0: return 0
    kps = [p for (p, _) in factorize2(k)]
    V = [n//i for i in range(1, sqrtInt(n)+1)]
    V += list(range(V[-1]-1, -1, -1))
    S1 = {i:i for i in V}
    for p in kps:
        for v in V:
            if v<p: break
            S1[v] -= S1[v//p]
    return S1[n]


def coprime(lb, ub, pfacts, i=0, prd=1):
    """
    Return the # of ints in (lb, ub] not divisible by the primes in pfacts.
    """
    if i == len(pfacts):
        return ub//prd - lb//prd
    return coprime(lb, ub, pfacts, i+1, prd) - coprime(lb, ub, pfacts, i+1, prd*pfacts[i])

class PhiTest(unittest.TestCase):
    def test_phi(self):
        for nn in range(1, 100):
            nresid = sum(1 for kk in range(1, nn+1) if gcd(nn, kk) == 1)
            self.assertEqual(nresid, phi(nn))
    def partial_totient_test(self):
        partial_sum = 0
        for nn in range(1, 101):
            partial_sum += phi(nn)
            self.assertEqual(partial_sum, partial_totient(nn))
            self.assertEqual(partial_sum._partial_totient_alternate(nn))
    def test_coprime(self):
        pfacts = [2, 3, 5, 7]
        xx = 2*3*5*7
        for lb in range(1, 151):
            for ub in range(lb+1, 151):
                self.assertEqual(coprime(lb, ub, pfacts), partial_totient(ub, xx)-partial_totient(lb, xx), '(%d, %d]'%(lb, ub))

class SumSigmaTest(unittest.TestCase):
    def test_sum_sigma1(self):
        for nn in range(1, 100):
            ss1 = sum([divisor_function(1, kk) for kk in range(1, nn+1)])
            self.assertEqual(ss1, sum_sigma1(nn))
    def test_sum_sigma0(self):
        for nn in range(1, 100):
            ss0 = sum([divisor_function(0, kk) for kk in range(1, nn+1)])
            self.assertEqual(ss0, sum_sigma0(nn))

