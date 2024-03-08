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
        for transform, m in ((  hadamard_transform, [[ 1, 1],
                                                     [ 1,-1]]),
                             (        or_transform, [[ 1, 0],
                                                     [ 1, 1]]),
                             (inverse_or_transform, [[ 1, 0],
                                                     [-1, 1]])):
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

