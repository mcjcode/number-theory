#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

import sympy

from utilities import (
    factorize2,
)

from cyclotomic_fields import (
    CyclotomicInteger,
    cyclotomic_polynomial,
    lucas_formula,
    gauss_formula,
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
            
    def test_lucas_formula(self):
        """
        If n=1(4) is squarefree, then the lucas polynomials
        C and D satisfy
        
        C_n^2 - n*x*D_n^2 = Phi_n

        Check that this holds for a range of n values

        If n=3(4) is squarefree, then the lucas polynomials
        C_n and D_n satisfy

        C_n^2 - n*x*D_n^2 = Phi_{2n}
        """

        x = sympy.abc.x
        for n in range(5, 100, 4):            
            if all(e==1 for _, e in factorize2(n)):
                C, D = lucas_formula(n)
                Phi = cyclotomic_polynomial(n)
                self.assertEqual(C**2 - n*x*D**2, Phi)

        for n in range(3, 100, 4):
            if all(e==1 for _, e in factorize2(n)):
                C, D = lucas_formula(n)
                Phi = cyclotomic_polynomial(2*n)
                self.assertEqual(C**2 - n*x*D**2, Phi)

    def test_gauss_formula(self):
        x = sympy.abc.x
        for n in range(5, 100, 4):
            if all(e==1 for _, e in factorize2(n)):
                A, B = gauss_formula(n)
                Phi = cyclotomic_polynomial(n)
                self.assertEqual(A**2 - n*B**2, 4*Phi)
