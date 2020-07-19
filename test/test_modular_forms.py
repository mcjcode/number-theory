#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
import numpy as np
from modular_forms import (
    _rectangle_n_points,
    eisenstein,
    unrestricted_eisenstein,
    poincare,
)


class RectangleNPointsTest(unittest.TestCase):

    def runTest(self):
        for nn in range(1, 10):
            pts = list(_rectangle_n_points(nn))
            npts = len(pts)
            nuniqpts = len(list(set(pts)))

            self.assertEqual(npts, 8 * nn)
            self.assertEqual(nuniqpts, 8 * nn)


class EisensteinTest(unittest.TestCase):

    def runTest(self):
        pass

    @staticmethod
    def _random_pt():
        while True:
            x = np.random.rand() - 0.5
            y = np.random.rand()
            if x ** 2 + y ** 2 <= 1 and (x - 1) ** 2 + y ** 2 >= 1 and (x + 1) ** 2 + y ** 2 > 1:
                break
        return x + 1.0j * y

    def _test_modularity(self, f):
        z = self._random_pt()
        for kk in [4, 5]:
            self.assertAlmostEqual(f(kk, z + 1), f(kk, z))
            self.assertAlmostEqual(f(kk, -1 / z), z ** (2 * kk) * f(kk, z))

    def test_eisenstein(self):
        """
        Verify that the Eisenstein series are modular forms.
        """
        self._test_modularity(lambda kk, zz: eisenstein(kk, zz))
        self._test_modularity(lambda kk, zz: eisenstein(kk, zz))

    def test_unrestricted_eisenstein(self):
        """
        Verify that the unrestricted Eisenstein series are modular forms.
        """
        self._test_modularity(lambda _k, _z: unrestricted_eisenstein(_k, _z))
        self._test_modularity(lambda _k, _z: unrestricted_eisenstein(_k, _z))

        zz = self._random_pt()
        for kk in range(3, 7):
            zeta2k = sum(1. / n ** (2 * kk) for n in range(100000, 0, -1))
            self.assertAlmostEqual(unrestricted_eisenstein(kk, zz), zeta2k * eisenstein(kk, zz), places=6)

    def test_poincare(self):
        """
        Verify that the Poincare series are modular forms.
        """
        self._test_modularity(lambda kk, zz: poincare(kk, 0, zz))
        self._test_modularity(lambda kk, zz: poincare(kk, 1, zz))
