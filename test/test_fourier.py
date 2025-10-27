#!/usr/bin/env python -i
#
# test_fourier.py

import unittest
import numpy as np
from fourier import (
    generalized_hadamard_transform,
    hadamard_transform,
    or_transform,
    inverse_or_transform,
    fourierZ6,
    fourierZ10,
    dft_reference,
    _dft,
    dft,
)


class FourierTests(unittest.TestCase):

    @staticmethod
    def mul(v, w):
        n = len(v)
        retval = [0]*n
        for i in range(n):
            for j in range(n):
                ij = i*j % n
                retval[ij] = (retval[ij] + v[i]*w[j])
        return retval

    @staticmethod
    def check_orthogonality(m):
        m = m.transpose()
        # check that the system of vs is correct
        for v in m:
            for w in m:
                x = FourierTests.mul(v, w)
                if np.linalg.norm(x):
                    if np.any(v!=w) or np.linalg.norm(np.array(v)-np.array(x)):
                        return False
        return True

    def test_fourierZ6(self):
        self.assertTrue(FourierTests.check_orthogonality(fourierZ6()))

    def test_fourierZ10(self):
        self.assertTrue(FourierTests.check_orthogonality(fourierZ10()))


class GeneralizedHadamardTest(unittest.TestCase):
    """
    Make sure that the generalize_hadamard_transform
    gives the same results as the individual
    hard-coded transforms
    """
    def test_transforms(self):
        for transform, m in ((hadamard_transform,   [[1,  1],
                                                     [1, -1]]),
                             (or_transform,         [[1,  0],
                                                     [1,  1]]),
                             (inverse_or_transform, [[1,  0],
                                                    [-1,  1]])):
            m = np.array(m, dtype=np.int8)
            n = 2
            for _ in range(6):
                vals = []
                for method in transform, (lambda a: generalized_hadamard_transform(a, m)):
                    vs = np.identity(n, dtype=np.int8)
                    val = []
                    for v in vs:
                        method(v)
                        val.append(v.copy())
                    m = np.array(val)
                    vals.append(np.array(val))
                self.assertTrue(np.all(vals[0]==vals[1].transpose()))
                n *= 2


def test_dft_reference_convolve():
    def convolve(u, v, p):
        n = u.shape[0]
        w = np.zeros((n,), dtype=object)
        for i in range(n):
            for j in range(n):
                w[i] = (w[i] + u[j]*v[(n+i-j)%n]) % p
        return w

    for p in [7, 13, 19, 31]:

        p = 7
        n = p-1
        u = np.random.randint(0, p, n)
        v = np.random.randint(0, p, n)

        uv1 = convolve(u, v, p)

        ut = dft_reference(u, p)
        vt = dft_reference(v, p)

        uv2 = (dft_reference(ut*vt, p, inverse=True)*pow(n, -1, p))%p

        assert np.all(uv1==uv2)


def test_inner_dft():
    """
    Test the _dft routine against the reference dft implementation
    """

    n = 30
    p = 31
    qs = [2, 3, 5]
    a = 3  # a primitive root mod 30
    for inverse in False, True:
        for v1 in np.eye(n):
            v2 = _dft(v1, qs=qs, p=p, a=a, inverse=inverse)
            v3 = dft_reference(v1, p, inverse)
            assert np.all(v2==v3)


def test_dft():
    """
    Test dft against the reference implementation
    """
    n = 30
    p = 31
    for inverse in False, True:
        for v in np.eye(n):
            w1 = dft_reference(v, p, inverse=inverse)
            w2 = dft(v, p, inverse=inverse)
            assert np.all(w1==w2)
