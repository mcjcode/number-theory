#!/usr/bin/env python -i
# -*- coding: utf-8 -*-
"""
General purpose, factorization and modular
arithmetic routines.
"""

import unittest
import math

from utilities import (
    argmax,
    powerset,
    isprime,
    isprime_miller_rabin,
    squarefree,
    gcd,
    bezout,
    modinv,
    sqrtInt,
    cbrtInt,
    issq,
    step_modp_pascal,
    order,
    factorize2,
    primitive_root,
    cardano,
    divisors,
)

from prime_sieve import segmented_sieve


class UtilitiesTest(unittest.TestCase):

    def runTest(self):
        pass

    def test_argmax(self):
        f = math.sin
        xs = [2*math.pi*x/100 for x in range(100)]
        self.assertEqual(f(argmax(f, xs)), max(map(f, xs)))
        
    def test_powerset(self):
        nn = len(list(powerset([1, 2, 3, 4, 5])))
        self.assertEqual(nn, 2**5, 'number of subsets should = 32. was %d' % (nn, ))

    def test_isprime(self):
        self.assertEqual(False, isprime(0), 'zero is not a (maximal) prime')
        self.assertEqual(False, isprime(1), 'units are not prime')
        self.assertEqual(True, isprime(7), 'yes, 7 is a prime')
        self.assertEqual(False, isprime(49), 'a non-unit square is not prime')
        self.assertEqual(False, isprime(91), '91=7*13 is not prime')
        self.assertEqual(True, isprime(-7), '(some) negative numbers are prime')
        self.assertRaises(TypeError, isprime, 7.0)

    def test_miller_rabin(self):
        for n in range(1000):
            self.assertEqual(isprime_miller_rabin(n), isprime(n))
            
    def test_squarefree(self):
        self.assertEqual(True, squarefree(1), '1 is square free')
        self.assertEqual(True, squarefree(-1), '-1 is square free')
        self.assertEqual(False, squarefree(4), '4 is not square free')
        self.assertEqual(False, squarefree(18), '18 is not square free')

    def test_gcd(self):
        self.assertEqual(gcd(2*3*5, 3*5*7), 3*5)

    def test_bezout(self):
        for a in range(-100, 100):
            for b in range(-100, 100):
                x, y = bezout(a, b)
                g = gcd(a, b)
                self.assertEqual(g, x*a + y*b, f'gcd ({g}) != x*a+y*b ({x}*{a}+{y}*{b})')

    def test_modinv(self):
        for m in range(2, 100):
            for a in range(1,m):
                if gcd(m, a)==1:
                    b = modinv(m, a)
                    self.assertEqual((a*b)%m, 1, f'{a}*{b}%{m} != 1')

    def test_sqrtInt(self):
        for n in range(1, 1000):
            n2 = n*n
            self.assertEqual(sqrtInt(n2-1), n-1)
            self.assertEqual(sqrtInt(n2),   n)
            self.assertEqual(sqrtInt(n2+1), n)
            
    def test_cbrtInt(self):
        for n in range(1000):
            cbrtn = cbrtInt(n)
            self.assertTrue(cbrtn**3 <= n < (cbrtn+1)**3)

    def test_issq(self):
        for n in range(1, 10000):
            self.assertTrue(issq(n*n))
            self.assertTrue(not issq(n*n+1))
            
    def test_step_modp_pascal(self):
        for p in [2, 3, 5]:
            for k in range(1,3):
                row = [(0,1)]
                N = p**k
                for i in range(N-1):
                    row = step_modp_pascal(row, p)
                self.assertEqual(len(row), N, f'row should be all non-zero p={p}, k={k}, {row}')
            
    def test_factorize2(self):
        for n0 in range(2,10_000):
            n = n0*n0-1
            n2 = 1
            for p, e in factorize2(n):
                self.assertTrue(isprime(p), f'non-prime {p} found in factorization')
                n2 *= p**e
            self.assertEqual(n, n2)

    def test_primitive_root(self):
        for p in segmented_sieve(100):
            with self.subTest(p=p):
                self.assertEqual(order(primitive_root(p), p), p-1)
                
    def test_trial_division_exception(self):
        with self.assertRaises(Exception):
            f = factorize2(41*43, 10)
            next(f)

    def test_cardano(self):
        for x in range(1, 10):
            self.assertAlmostEqual(cardano(0.0, -x), math.cbrt(x), places=10)

    def test_divisors(self):
        assert sum(divisors(10))==1+2+5+10
        assert sum(divisors({2:3, 3:2}))==(1 + 2 + 2**2 + 2**3)*(1 + 3 + 3**2)
