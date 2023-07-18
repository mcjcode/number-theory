#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

import sympy

from cyclotomic_fields import (
    CyclotomicInteger,
    cyclotomic_polynomial,
)

from multiplicative import (
    phi,
)


class CyclotomicIntegerTest(unittest.TestCase):

    def test_axioms(self):
        nn = 13
        zero = CyclotomicInteger.zero(nn)
        one = CyclotomicInteger.one(nn)

        xx = CyclotomicInteger.random(nn)
        yy = CyclotomicInteger.random(nn)
        zz = CyclotomicInteger.random(nn)

        self.assertEqual(one * xx, xx)
        self.assertEqual(zero + xx, xx)
        self.assertEqual(xx + yy, yy + xx)
        self.assertEqual(xx * yy, yy * xx)
        self.assertEqual(xx * (yy + zz), xx * yy + xx * zz)

    def test_norm(self):
        nn = 13

        one = CyclotomicInteger.one(nn)
        zeta = CyclotomicInteger.zeta(nn)

        xx = one - zeta
        yy = one + zeta

        self.assertEqual(zeta ** nn, one)
        self.assertEqual(xx.norm(), nn)
        self.assertEqual(yy.norm(), 1)

class CyclotomicPolynomialTest(unittest.TestCase):

    def test_degree(self):
        """
        The degree of Phi_n(x) should equal phi(n)
        """
        for n in range(1,100):
            p = cyclotomic_polynomial(n)
            d = phi(n)
            self.assertEqual(p.degree(), d)

    def test_products(self):
        """
        Product_{d|n} Phi_d(n) == x**n - 1.

        Check that this holds.
        """
        x = sympy.abc.x
        for n in range(1, 100):
            lhs = sympy.poly(x**n - 1, x, domain='ZZ')
            rhs = sympy.poly(1, x, domain='ZZ')
            for d in range(1, n+1):
                if n%d==0:
                    rhs *= cyclotomic_polynomial(d)
            self.assertEqual(lhs, rhs)
            
