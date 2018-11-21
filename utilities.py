#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import math; from math import sqrt
import gaussian_integer
from gaussian_integer import GaussianInteger
import numpy as np

def assert_equal(expected, actual, msg) :
	if expected==actual :
		print '.',
		return True
	else :
		raise Exception('testing "%s": %s expected, %s actual' % (msg, expected, actual))

def assert_exception(f,msg) :
	"""
	Test whether the code in the function f
	throws an exception
	"""
	try :
		f()
	except Exception :
		print '.',
		return
	raise Exception(msg)

def isprime(p) :
	if type(p) != int :
		raise ValueError('%s is %s, but should be int' % (p,type(p)))
	
	if p < 0  :
		p = -p

	if p == 0 or p == 1 : 
		return False

	trdiv = 2
	while trdiv*trdiv <= p :
		if p % trdiv == 0 :
			return False
		trdiv += 1
	return True

def issq(nn) :
	return nn >= 0 and nn == int(sqrt(nn))**2
	
def factorize(n) :
	"""
	Return a non-decreasing list of primes whose product is n.
	"""
	if n==1 :
		return []
	q = 2
	while q*q <= n :
		if n % q == 0 :
			return [q] + factorize(n/q)
		q += 1
	return [n]

def squarefree(mm) :
	"""
	Return True if mm is square free, False otherwise.
	"""
	factors = factorize(abs(mm))
	for i in xrange(len(factors)-1) :
		if factors[i] == factors[i+1] :
			return False
	return True

def order(a,p) :
	retval = 1
	one = type(a)(1)
	cnt = 1
	b = a
	while not b == one :
		b = (b * a) % p
		cnt += 1
	return cnt

def primitive_root(p) :
	for a in xrange(2,p) :
		if order(a,p) == p-1 :
			return a

def prod(lst) :
	return reduce(lambda x, y: x * y, lst, 1)

def modpow(a,k,p) :
	"""
	Return a**k (mod p)
	
	a can be of any type that has a multiplicative
	identity and supports multiplication and modding.
	"""
	retval = type(a)(1)
	cnt = 0
	while cnt < k :
		retval = (retval * a) % p
		cnt += 1
	return retval

def modpow2(a,k,p) :
	if k == 0 :
		return type(a)(1)
	tail = (modpow2(a,k/2,p) ** 2) % p
	if k % 2 == 1 :
		return (a * tail) % p
	else :
		return tail % p
		
	
def legendre_ch(p) :
	"""
	Return the mod p legendre character, a function 'ch' such that
	
	ch(a) =  0  if p divides a
	        +1  if a is a quadratic residue
	        -1  if a is a quadratic non-residue
 	"""
	if not isprime(p) or p==2 :
		raise ValueError("%d is not an odd prime." % (p,))
	def ch(a) :
		if a % p == 0 :
			return 0
		rr = modpow2(a,(p-1)/2,p)
		return (-1) if (rr==p-1) else +1
	return ch


	
def chi(pi) :
	N = pi.norm()
	xp = (N-1)/4
	retval = GaussianInteger(1)
		
def JacobiSum(chi1,chi2,p) :
	"""
	Return the Jacobi sum of two mod p characters
	with values in the GaussianIntegers (so the
	characters are either the trivial character,
	the non-trivial quadratic character or one of
	the two non-trivial quartic characters.)
	"""
	retval = GaussianInteger(0)
	for a in xrange(p) :
		retval = retval + chi1(a)*chi1(1-a)
	return retval
	
def round(x) :
	
	s=np.sign(x)
	return s*int(s*x+0.5)
	
def J(p) :
	"""
	Returns the Jacobi sum J(chi,chi) associated with the
	quartic character which maps the smallest primitive
	root mod p to the pure imaginary unit i.
	"""
	a = primitive_root(p)
	zeta = GaussianInteger(0,1)
	avals = [0] + [modpow(a,k,p) for k in xrange(1,p)]
	achars = [GaussianInteger(0)] + [zeta**(k%4) for k in xrange(1,p)]
	bvals = [(1-aa)%p for aa in avals]
	binds = [avals.index(bb) for bb in bvals]
	bchars = [achars[bi] for bi in binds]	
	retval = GaussianInteger(0)
	for val in (ac*bc for (ac,bc) in zip(achars,bchars)) :
		retval = retval + val
	return retval

def gcd(a,b) :
	if a==0 or b==0 :
		return a+b
	while a%b != 0 :
		a, b = b, a%b
	return b
	
if __name__=='__main__' :
	assert_equal(False,isprime(0),'zero is not a (maximal) prime')
	assert_equal(False,isprime(1),'units are not prime') 
	assert_equal(True,isprime(7),'yes, 7 is a prime')
	assert_equal(False,isprime(49),'a non-unit square is not prime')
	assert_equal(False,isprime(91),'91=7*13 is not prime')
	assert_equal(True,isprime(-7),'(some) negative numbers are prime')
	assert_equal(True,squarefree(1),'1 is square free')
	assert_equal(True,squarefree(-1),'-1 is square free')
	assert_equal(False,squarefree(4),'4 is not square free')
	assert_equal(False,squarefree(18),'18 is not square free')
	assert_exception(lambda : isprime(7.0), 'isprime should throw an exception if passed a non-int')
	print

