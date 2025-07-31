#!/usr/bin/env python -i
# -*- coding: utf-8 -*-
"""
Various multiplicative arithmetic functions
and their summations.
"""

import numpy as np

from utilities import (
    prod,
    factorize2,
    gcd,
    sqrtInt,
)


def phi(nn):
    r"""
    :param nn: a positive integer
    :return: the number of natural numbers 1 <= kk < nn
             that are relatively prime to nn.
    """
    factors = factorize2(nn)
    return prod([pp**(kk-1)*(pp-1) for (pp, kk) in factors])


def mu(n: int) -> int:
    n = sum(1 for _ in factorize2(n))
    return (-1)**(n%2)


def divisor_function(kk, nn):
    r"""
    :param kk: the exponent
    :param nn: a positive integer
    :return: :math:`\sum_{d|n} d^k`
    """
    factors = factorize2(nn)
    return prod([(ee+1) if kk == 0 else (pp**(kk*(ee+1))-1)//(pp**kk-1)
                 for (pp, ee) in factors])


def sumrange(a, b):
    return b*(b+1)//2 - a*(a+1)//2


def sum_sigma0(n):
    r"""
    :param n: a positive integer
    :return: :math:`\displaystyle\sum_{k\in[1\dots n]} d(k)`
    """
    u = sqrtInt(n)
    return 2*sum(n//k for k in range(1, u+1)) - u**2

def sum_sigma1(n):
    r"""
    :param n: a positive integer
    :return: :math:`\displaystyle\sum_{k\in[1\dots n]} \sigma(k)`
    """
    sqrtk = sqrtInt(n)
    part1 = sum(d*(n//d) for d in range(1, n//sqrtk+1))
    part2 = sum(sumrange(n//(d+1), n//d)*d for d in range(1, sqrtk))
    return part1+part2



def partial_totient(n: int, k: int, ps: list[int] = []) -> int:
    r"""
    :return: how many k in [1,n] are relatively prime to k.
    n and k should be positive
    """
    ps = ps or [p for p,_ in factorize2(k)]
    nps = len(ps)
    xs = [0]*(2**nps)
    xs[0] = (n, +1)
    retval = n
    i = 1
    for p in ps:
        for j in range(i):
            x, s = xs[j]
            if x>=p:
                retval -= s*(xp:=x//p)
                xs[i] = (xp, -s)
                i += 1
    return retval
    #    xs += [(x//p, -s) for x,s in xs if x>=p]
    #return sum(x*s for x,s in xs)


def _partial_totient_alternate(n: int, k: int) -> int:
    r"""
    :return: how many k in [1,n] are relatively prime to k.
    n and k should be positive
    """
    V = [n//i for i in range(1, sqrtInt(n)+1)]
    V += list(range(V[-1]-1, -1, -1))
    S1 = {i:i for i in V}
    for p, _ in factorize2(k):
        for v in V:
            if v < p:
                break
            S1[v] -= S1[v//p]
    return S1[n]


def coprime(lb: int, ub: int, ps: list[int]) -> int:
    r"""
    :param lb: a positive integer
    :param ub: a positive integer greater than lb
    :param ps: a list of prime numbers
    :return: how many k in (lb, ub] are not divisible by the p's in ps
    """    
    def f(i=0, prd=1):
        if i==len(ps):
            return ub//prd - lb//prd
        return f(i+1, prd) - f(i+1, prd*ps[i])
    return f()

def coprime0(ub: int, ps: list[int]) -> int:
    r"""
    :param ub: a positive integer greater than lb
    :param ps: a list of prime numbers
    :return: how many k <= ub are not divisible by the p's in ps
    """
    nps = len(ps)
    def f(i=0, prd=1):
        if i==nps:
            return ub//prd
        return f(i+1, prd) - f(i+1, prd*ps[i])
    return f()
