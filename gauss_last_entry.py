#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

from __future__ import print_function

from utilities import isprime
from finite_field import (
    FiniteField,
    count_curve_points_affine
    )

print(u"""
For the curve given by the affine equation

x\u00B2 + y\u00B2 + x\u00B2y\u00B2 = 1

calculate the number of points on curve over
the finite field Z/pZ.

If p\u22611(4), and we write p=a\u00B2+b\u00B2 with
a odd and b even and

a\u22613(4) if b\u2261 2(4)
a\u22611(4) if b\u2261 0(4)

then Gauss conjectured that the number of points
on this curve, including the 2 add'l points at
infinity, should equal

p-1-2a.

Do a brute force calculation and verify this for
such primes up to 199.
""")

for p in range(5, 200, 4):
    if isprime(p):
        ff = FiniteField(p, 1)
        Np = count_curve_points_affine(lambda x, y: x**2+y**2+x**2*y**2-ff.one(), ff)
        js = ff.jacobi_sum(4)
        print('%4d %5d+2 = p-1%+3d, J=%+2.0f%+2.0fi' % (p, Np, (Np+2-p+1), js[0], js[1]))
