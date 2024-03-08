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
)

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

