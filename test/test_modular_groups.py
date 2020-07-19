#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

import numpy as np

from modular_groups import (
    psl2q_order,
    coset_reps,
    coset_reps_alt,
    random_complex_gaussian,
    poincare_disk_to_halfplane,
    halfplane_to_poincare_disk,
)


class CosetRepsUnitTest(unittest.TestCase):
    def runTest(self):
        for q in range(2, 8):
            expected = psl2q_order(q)
            actual1 = len(list(coset_reps(q))) / (1 if q == 2 else 2)
            actual2 = len(coset_reps_alt(q))
            self.assertEqual(expected, actual1)
            self.assertEqual(expected, actual2)


class HalfplaneToDiskUnitTest(unittest.TestCase):
    def runTest(self):
        np.random.seed(17)
        for i in range(100):
            z1 = random_complex_gaussian()
            z1 = z1.real + 1.0j*abs(z1.imag)
            z2 = poincare_disk_to_halfplane(halfplane_to_poincare_disk(z1))
            self.assertAlmostEqual(z1, z2)
        for i in range(100):
            th = np.random.uniform()*2*np.pi
            rr = np.sqrt(np.random.uniform())
            z1 = rr * (np.cos(th) + 1.0j*np.sin(th))
            z2 = halfplane_to_poincare_disk(poincare_disk_to_halfplane(z1))
            self.assertAlmostEqual(z1, z2)
