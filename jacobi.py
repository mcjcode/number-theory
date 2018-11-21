#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import numpy as np
import utilities; from utilities import isprime, primitive_root

k = 3

# Calculate the Jacobi sums of a cubic character
# for some finite fields of order p^2
#
		
zeta3 = np.cos(2*np.pi/3) + 1j*np.sin(2*np.pi/3)

# These are all going to
for p in (pp for pp in xrange(7,200,12) if isprime(pp)) :

	Np = p if (p % 4 == 1) else p*p

	#print 'Finite Field with %d^2=%d elements' % (p,p*p)
	def mul(a,b) :
		# Here's how we multiply Gaussian integers
		return ((a[0]*b[0]-a[1]*b[1]) % p, (a[0]*b[1]+a[1]*b[0]) % p)
	
	def add(a,b) :
		# Here's how we add Gaussian integers
		return ((a[0]+b[0]) % p), ((a[1]+b[1]) % p)
	
	def ord(a) :
		retval = 1
		b = a
		while b != (1,0) :
			b = mul(b,a)
			retval += 1
		return retval
	
	def prim_elt() :
		for x in xrange(p) :
			for y in xrange(p) :
				elt = (x,y)
				if elt == (0,0) :
					continue
				if ord(elt) == p*p-1 :
					return elt

	xx =  prim_elt()
	powers = [(0,0),(1,0)]
	cnt = 0
	while cnt < Np-2 :
		powers.append(mul(powers[-1],xx))
		cnt += 1
	
	chars = [0.0] + map(lambda n : zeta3**n, xrange(p*p))
	avals = powers
	bvals = map( lambda xx : add((1,0),mul((-1,0),xx)), avals)
	bindices = [avals.index(bb) for bb in bvals]
	bchars = [chars[bi] for bi in bindices]
	jsum = sum( aa*bb for (aa,bb) in zip(chars,bchars) )
	
	#express jsum as an integral combination of 1 and zeta3=(-1+sqrt(3))/2
	#
	c2 = jsum.imag/(np.sqrt(3.0)/2.0)
	tmp = jsum - c2*zeta3
	c1 = tmp.real
	
	print '%6d %10d %35s %7.1f %7.1f %35s %s' % (p, p*p, jsum, c1, c2, c1+c2*zeta3, '%d+%di' % xx)
		
def f(x,y,p) :
	return (y*y - (x-1)*x*(x+1)) % p

for p in []: #xrange(2,1000) :
	if not isprime(p) :
		continue
		#[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61] :
	cnt = 0
	for x in range(p) :
		for y in range(p) :
			cnt += f(x,y,p)==0
	bound = 2*int(np.sqrt(p))
	if abs(cnt-p) >= bound :
		print '%5d %5d %5d %5d' % (p, cnt, cnt - p, bound)
	
			
