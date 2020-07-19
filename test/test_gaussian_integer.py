#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from gaussian_integer import GaussianInteger


class GaussianIntegerTest(unittest.TestCase):
    def test_axioms(self):
        one_elem = GaussianInteger.one()
        zero_elem = GaussianInteger.zero()
        xx = GaussianInteger.random()
        yy = GaussianInteger.random()
        zz = GaussianInteger.random()

        self.assertEqual(one_elem * xx, xx)
        self.assertEqual(zero_elem + xx, xx)
        self.assertEqual(xx + yy, yy + xx)
        self.assertEqual(xx * yy, yy * xx)
        self.assertEqual(xx * (yy + zz), xx * yy + xx * zz)

    def test_two(self):
        zero_elem = GaussianInteger.zero()
        xx = GaussianInteger.random()
        self.assertEqual(xx.conj().conj(), xx)
        self.assertEqual(xx - zero_elem, xx)
        self.assertEqual(-(-xx), xx)
        self.assertEqual(xx**2, xx*xx)
        self.assertEqual((xx * xx.conj()).imag(), 0)
