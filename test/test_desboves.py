#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from desboves import DesbovesCurvePoint


class DesbovesUnitTest(unittest.TestCase):
    def test_one(self):
        zero = DesbovesCurvePoint(13, (1, -1, 0))
        p = DesbovesCurvePoint(13, (2, 7, 3))
        self.assertEqual(zero, zero)
        self.assertEqual(zero+p, p)
        self.assertEqual(p+zero, p)
        self.assertEqual(p+(-p), zero)

    def test_integerMult(self):
        zero = DesbovesCurvePoint(13, (1, -1, 0))
        p = DesbovesCurvePoint(13, (2, 7, 3))
        self.assertEqual(0*p, zero)
        self.assertEqual(1*p, p)
        self.assertEqual(2*p, p+p)
        self.assertEqual(3*p, p+p+p)
