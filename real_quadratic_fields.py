#!/usr/bin/env python -i
# -*- coding: utf-8 -*-
import math; from math import sqrt, floor
import utilities; from utilities import factorize, prod, isprime, gcd
import fractions; from fractions import Fraction
import itertools; from itertools import islice

pi = 3.14159265358979323

"""
"""

def cont_frac(m) :
	"""
	For a positive floating point number 'm',
	return a generator for the sequence of
	integers that give the continued fraction
	representation of
	the number.
	
	          1    1    1
	m = a0 + —-   —-   —-
	         a1 + a2 + a3 + ...
	         
	You probably shouldn’t use this
	if you really care about all of the terms
	since rounding eventually corrupts the
	process.
	
	For numbers in quadratic number fields,
	use ‘cont_frac_quad’ instead.
	"""
	while True :
		f = int(floor(m))
		yield f
		if m==f :
			return
		m = 1./(m-f)

def divint(num,denom,prec=15) :
	"""
	Return a floating point number representing num/denom
	accurate to prec decimal digits of accuracy.
	
	Use this when the integers num and denom are not
	representable as floating point numbers (too big) but
	the quotient is.
	
	Don’t bother using prec>15 as the arithmetic is
	single precision floats.
	"""
	ans = 0.0
	for i in xrange(prec) :
		q, r = divmod(num,denom)
		ans += q/(10.0**i)
		num = (10*r)
	return ans
	
def floor_rem(x) :
	fl = floor(x)
	return int(fl), x - fl
	
def floor_custom(a,b,c,d) :
	j, epsj = floor_rem(sqrt(d))
	k, epsk = floor_rem(floor((a+b*j)/c))
	
	
def cont_frac_quad(a,b,c,d) :
	"""
	Return the continued fraction of (a+b*sqrt(d))/c
	"""
	while True :
		## a,b and c may be large, but a/c and b/c
		## should be reasonable, so do int division
		## first before trying to convert to floats
		
		m = int(divint(a,c,15) + sqrt(d)*divint(b,c,15))
		## m = int(floor( (a+b*sqrt(d))/c ))
		yield m
		a1 = c*(a-c*m)
		b1 = -c*b
		c1 = ((a-c*m)**2 - d*b**2)
		a, b, c = a1, b1, c1

def approximants(d) :
	h0, k0 = 0, 1
	h1, k1 = 1, 0
	
	m = sqrt(d)
	
	for ai in cont_frac_quad(0,1,1,d) :
		h = ai*h1 + h0	
		k = ai*k1 + k0
		c = h**2 - d*k**2
		yield QuadInt(d,h-k,2*k) if ((d%4) == 1) else QuadInt(d,h,k)

		if c == 1 :
			break
		h0, k0 = h1, k1
		h1, k1 = h,  k
		
def discriminant( d ) :
	"""
	Return the discriminant of the quadratic field Q[\sqrt{d}]
	"""
	factors = factorize(abs(d))
	if d % 4 == 1 :
		return abs(d)
	else :
		return 4*abs(d)
		
def minkowski_bound( d ) :
	"""
	For an integer d let K=Q[\sqrt{d}].
	
	All ideal classes in O_K have a representative J
	with ||J|| <= the minkowski bound.
	"""
	disc = discriminant( d )
	if d > 0 :
		return int( 0.5 * sqrt( disc ) )
	else :
		return int( 0.5 * (4.0/pi) * sqrt(disc) )

def factorize_in(p,d) :
	"""
	If d is a squarefree integer, factorize the
	prime integer n in the ring of integers in the
	quadratic field extension Q[\sqrt{d}].
	"""
	if d % 4 != 1 :
		## the minimal polynomial is x**2 - d = 0
		x = 0
		while x < p :
			if (x**2-d) % p == 0 :
				return (p,QuadInt(d,-x,1))
			x += 1
		return (p,)
	else :
		## the minimal polynomial is x**2 - x - (d-1)/4
		zz = (d-1)/4
		x = 0
		while x < p :
			if (x*(x-1) - zz) % p == 0 :
				return (p,QuadInt(d,-x,1))
			x += 1
		return (p,)

def squares_mod_d(d) :
	return sorted({(n**2)%d for n in xrange(d)})

def issquare(nn) :
	return (nn >= 0) and (int(sqrt(nn))**2 == nn)

class QuadInt(object) :
	def __init__(self, d, a, b) :
		self.d = d
		self.a = a
		self.b = b
	def __str__(self) :
		if self.d % 4 != 1 :
			part1 = u"%d" % (self.a)
			sgn   = "+" if self.b > 0 else "-"
			if abs(self.b) == 1 :
				part2 = u"\u221A%d" % (self.d,)
			else :
				part2 = u"%d\u221A%d" % (abs(self.b),self.d)
			return u" ".join([part1,sgn,part2])
		else :
			if self.b == 1 :
				return u"%d + (1+\u221A%d)/2" % (self.a,self.d)
			else :
				return u"%d + %d(1+\u221A%d)/2" % (self.a,self.b,self.d)
	def __repr__(self) :
		return self.__str__()
	def __add__(self,other) :
		if self.d != other.d :
			raise ValueError("integers are not from same field")
		return QuadInt(self.d,self.a+other.a,self.b+other.b)
	def __mul__(self,other) :
		if self.d != other.d :
			raise ValueError("integers are not from same field")
		else :
			return QuadInt(self.d,self.a*other.a-self.d*self.b*other.b, self.a*other.b-self.b*other.b)


	def norm(self) :
		if self.d % 4 != 1 :
			return self.a**2 - self.d * self.b**2
		else :
			return self.a**2 + self.a * self.b + self.b**2 * ((1-self.d)/4)

def norm_search(p,d) :
	md = None
	if (d % 4) != 1 :
		for nn in range(1,20000) :
			for mm in [nn**2*d+p,nn**2*d-p] :
				if issquare(mm) :
					m2 = int(sqrt(mm))
					if (gcd(nn,p)==1) or (gcd(m2,p)==1) :
						md = QuadInt(d,m2,nn)
						break
	return md
					
def class_group_info(d) :
	mb = minkowski_bound(d)
	disc = discriminant(d)
	print("Discriminant = %d"%disc)
	print("Minkowski Bound = %d" % (mb,))
	split_primes = []
	for p in [pp for pp in range(2,mb+1) if isprime(pp)] :
		fact = factorize_in(p,d)
		if fact == (p,) :
			print(fact)
		else :
			md = norm_search(p,d)
			print(fact,md)
			if md == None :
				split_primes.append(p)
	
	if (d%4) != 1 and d>0 :
		print("approximant based elts and norms")
		print([(xx,xx.norm()) for xx in approximants(d)])
	
	if split_primes == [] :
		print "No non principal primes under the Minkowski bound!"
		return
	print("Split primes")
	print(split_primes)
	print("Split primes that are not squares mod %d" % (d,))
	sq = squares_mod_d(d)
	print( [ pp for pp in split_primes if (pp not in sq) and ((d-pp) not in sq)])

	for ii in range(len(split_primes)) :
		for jj in range(ii,len(split_primes)) :
			p = split_primes[ii]*split_primes[jj]
			md = norm_search(p,d)
			print(p,md)

	for ii in range(len(split_primes)) :
		for jj in range(ii,len(split_primes)) :
			for kk in range(jj,len(split_primes)) :		
				p = split_primes[ii]*split_primes[jj]*split_primes[kk]
				md = norm_search(p,d)
				print(p,md)
	
