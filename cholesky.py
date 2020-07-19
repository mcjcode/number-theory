#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import numpy as np
from fractions import Fraction


def cholesky(mm):
    """
    Take a symmetric matrix m and calculate an upper
    triangular matrix a s.t. tr(a)*m*a is diagonal.
    """    
    mat = mm.copy()
    n = mat.shape[0]
    a = np.array([[Fraction(0, 1)]*n]*n, dtype=Fraction)
    for i in range(n):
        a[i, i] = Fraction(1, 1)
    for i in range(n-1):
        for j in range(i+1, n):
            r1 = Fraction(mat[j, i], mat[i, i])
            r2 = Fraction(mat[i, j], mat[i, i])
            mat[j, :] = mat[j, :] - mat[i, :] * r1
            mat[:, j] = mat[:, j] - mat[:, i] * r2
            a[:, j] = a[:, j] - a[:, i] * r1
    return mat, a
