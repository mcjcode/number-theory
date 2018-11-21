#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import numpy as np
import utilities; from utilities import isprime
import finite_field; from finite_field import (
	FiniteField,
	count_diagonal_cubic_points) 

print u"""
For a fixed characteristic p (where 3 divides p-1), there
is a unique (up to conjugation) non-trivial cubic
character \u03C7 on F_{p^k}, and so the real part of
the Jacobi sum J(\u03C7,\u03C7) only depends on the
field F_{p^k}.

Calculate J(\u03C7,\u03C7) for range of values of k
to get an idea of how these values are related to
each other.
"""

theta = 2*np.pi/3
zeta = np.cos(theta)+1j*np.sin(theta)
pi=-1-3*zeta
for p in [7,13] :
	if isprime(p) :
		for ii in xrange(1,5) :
			ff = FiniteField(p,ii)
			#print
			#print ff.rpoly
			#print
			c0,c1=ff.J(3)
			#print c0,c1
			#Nk = count_projective_curve_points(lambda x,y,z:x**3+y**3+z**3, ff)
			Nk = count_diagonal_cubic_points(ff)
			#print Nk
			print u'%6d    %+.0f%+.0f\u03c9   %4.0f   %4d' % (p**ii, c0,c1, (p**ii+1)-(-1)**ii*(pi**ii + pi.conj()**ii).real, Nk)
			
