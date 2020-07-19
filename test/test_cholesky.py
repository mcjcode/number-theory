#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

import numpy as np
from fractions import Fraction
import cholesky


class CholeskyTest(unittest.TestCase):
    def runTest(self):
        aa = list(map(lambda xs: list(map(Fraction, xs)), [[5, 3, 3],
                                                           [3, 3, 3],
                                                           [3, 3, 8]]))
        m = np.array(aa)
        d, a = cholesky.cholesky(m)
        atma = a.transpose().dot(m).dot(a)
        self.assertEqual(sum((atma-d).ravel()), Fraction(0, 1))
