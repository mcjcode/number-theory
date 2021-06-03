#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from gauss_last_entry import gauss_polynomial, normalize
from finite_field import (
    FiniteField,
    count_curve_points_affine
    )


class GaussLastEntryTest(unittest.TestCase):
    def runTest(self):
        for p in [5, 13, 17, 29]:
            ff = FiniteField(p, 1)
            q = gauss_polynomial(ff)
            np = count_curve_points_affine(q, ff)
            js = ff.jacobi_sum(4)
            a, _ = normalize(js[0], js[1])
            self.assertEqual(np+2, p-1-2*a)
