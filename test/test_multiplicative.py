#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from utilities import gcd
from multiplicative import (
    phi,
    mu,
    little_omega,
    big_omega,
    partial_totient,
    _partial_totient_alternate,
    coprime,
    divisor_function,
    sum_sigma0,
    sum_sigma1,
)
from prime_sieve import segmented_sieve


class PhiTest(unittest.TestCase):

    def test_phi(self):
        for nn in range(1, 100):
            nresid = sum(1 for kk in range(1, nn+1) if gcd(nn, kk) == 1)
            self.assertEqual(nresid, phi(nn))

    def test_partial_totient(self):
        for nn in range(1, 101):
            self.assertEqual(phi(nn), partial_totient(nn, nn))
            self.assertEqual(phi(nn), _partial_totient_alternate(nn, nn))
            for kk in range(1, nn):
                self.assertEqual(partial_totient(nn, kk),
                                 _partial_totient_alternate(nn, kk))
            
    def test_coprime(self):
        pfacts = [2, 3, 5, 7]
        xx = 2*3*5*7
        for lb in range(1, 151):
            for ub in range(lb+1, 151):
                pt = partial_totient(ub, xx)-partial_totient(lb, xx)
                msg = '(%d, %d]' % (lb, ub)
                self.assertEqual(coprime(lb, ub, pfacts), pt, msg)


class SumSigmaTest(unittest.TestCase):

    def test_sum_sigma1(self):
        for nn in range(1, 100):
            ss1 = sum([divisor_function(1, kk) for kk in range(1, nn+1)])
            self.assertEqual(ss1, sum_sigma1(nn))

    def test_sum_sigma0(self):
        for nn in range(1, 100):
            ss0 = sum([divisor_function(0, kk) for kk in range(1, nn+1)])
            self.assertEqual(ss0, sum_sigma0(nn))

def test_mu():
    assert mu(1)==1
    for i in range(1, 100):
        assert abs(mu(i))<=1

    for p in segmented_sieve(100):
        assert mu(p)==-1
        assert mu(p**2)==0

        
def test_big_omega():
    assert big_omega(1)==0
    for p in segmented_sieve(100):
        assert big_omega(p)==1
        assert big_omega(p**2)==2

        
def test_little_omega():
    for p in segmented_sieve(100):
        assert little_omega(p)==1
        assert little_omega(p**2)==1
        