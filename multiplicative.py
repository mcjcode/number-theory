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
    return prod([pp**(kk-1)*(pp-1) for (pp,kk) in factors])

def divisor_function(kk,nn):
    """
    return the sum of d^k over all divisors of nn
    """
    factors = factorize2(nn)
    return prod([(ee+1) if kk==0 else (pp**(kk*(ee+1))-1)//(pp**kk-1) for (pp,ee) in factors])

def sumrange(a,b):
    return b*(b+1)//2 - a*(a+1)//2
    
def sum_sigma0(n):
    """
    return the sum of sigma0(k) for
    all k in [1..n]
    """
    sqrtk=sqrtInt(n)
    part1=sum(n//d for d in range(1,n//sqrtk+1))
    part2=sum(((n//d)-(n//(d+1)))*d for d in range(1,sqrtk))
    return part1+part2
     
def sum_sigma1(n):
    """
    return the sum of sigma1(k) for
    all k in [1..n]
    """
    sqrtk=sqrtInt(n)
    part1=sum(d*(n//d) for d in range(1,n//sqrtk+1))
    part2=sum(sumrange(n//(d+1),n//d)*d for d in range(1,sqrtk))
    return part1+part2
    
class PhiTest(unittest.TestCase):
    def test_phi(self):
        for nn in range(1,100):
            nresid = sum(1 for kk in range(1,nn+1) if gcd(nn,kk) == 1)
            self.assertEqual(nresid, phi(nn))

class SumSigmaTest(unittest.TestCase):
    def test_sum_sigma1(self):
        for nn in range(1,100):
            ss1 = sum([divisor_function(1,kk) for kk in range(1,nn+1)])
            self.assertEqual(ss1, sum_sigma1(nn))
    def test_sum_sigma0(self):
        for nn in range(1,100):
            ss0 = sum([divisor_function(0,kk) for kk in range(1,nn+1)])
            self.assertEqual(ss0, sum_sigma0(nn))
            
