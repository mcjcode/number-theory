#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from cyclotomic_fields import CyclotomicInteger


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
