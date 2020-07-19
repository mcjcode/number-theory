#!/usr/bin/env python -i
# -*- coding: utf-8 -*-
"""
General purpose, factorization and modular
arithmetic routines.
"""

import unittest

from utilities import (
    powerset,
    isprime,
    squarefree,
    gcd,
    euclidean_algorithm,
)


class UtilitiesTest(unittest.TestCase):

    def runTest(self):
        pass

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

    def test_squarefree(self):
        self.assertEqual(True, squarefree(1), '1 is square free')
        self.assertEqual(True, squarefree(-1), '-1 is square free')
        self.assertEqual(False, squarefree(4), '4 is not square free')
        self.assertEqual(False, squarefree(18), '18 is not square free')

    def test_gcd(self):
        self.assertEqual(gcd(2*3*5, 3*5*7), 3*5)

    def test_euclidean_algorithm(self):
        a = 89
        b = 55
        x, y = euclidean_algorithm(a, b)
        self.assertEqual(abs(gcd(a, b)), abs(x*a + y*b))

    def test_euclidean_algorithm2(self):
        for a in range(-100, +100):
            for b in range(-100, +100):
                x, y = euclidean_algorithm(a, b)
                self.assertEqual(gcd(a, b), x*a + y*b, 'gcd != x*a+y*b')
