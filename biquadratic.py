#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import gaussian_integer
from gaussian_integer import GaussianInteger

import utilities
from utilities import J, isprime, modpow

def maptochar(x,pi) :
	zero = GaussianInteger(0)
	roots = map(lambda nn : GaussianInteger(0,1)**nn, range(4))
	for root in roots :
		if (x-root) % pi == zero :
			return root
			
maxp = 75


flip = lambda xx : GaussianInteger(xx.real(),abs(xx.imag()))
#flip = lambda xx : xx if xx.imag()>0 else GaussianInteger(xx.real(),-xx.imag())

qprimes = [(flip(J(qq)) if (qq%4==1) else GaussianInteger(qq))
					 for qq in xrange(3,maxp,2) if isprime(qq)]
					 
pprimes = qprimes #[GaussianInteger(0,1), GaussianInteger(1,1)] + qprimes

#qprimes = [qq for qq in qprimes if qq.imag()==0]

def biquad(pi1,pi2) :
	if type(pi2) == int :
		pi2 = GaussianInteger(pi2)
		
	if pi2 == GaussianInteger(0) :
		return GaussianInteger(0)
		
	z = pi2 % pi1
	ch = modpow(z,(pi1.norm()-1)/4,pi1)
	return maptochar(ch,pi1)
	
class Character(object) :
		
	def __init__(self,xx) :
		if type(xx) == GaussianInteger :
			pi = xx
			self.f = lambda a : biquad(pi,a)
		else :
			self.f = xx

	def __call__(self,a) :
		return self.f(a)
		
	def __mul__(self,other) :
		def f(a) :
			return self.f(a)*other.f(a)
		return Character(f)		



def JacobiSum(chi1,chi2,p) :
	"""
	Return the Jacobi sum of two mod p characters.
	"""
	retval = GaussianInteger(0)
	for a in xrange(p) :
		a = GaussianInteger(a)
		retval = retval + chi1(a)*chi2(GaussianInteger(1)-a)
	return retval

pi = GaussianInteger(-3,2)

q2primes = [flip(J(qq)) for qq in xrange(5,maxp,4) if isprime(qq)]

for pi1 in q2primes :
	for pi in [pi1,GaussianInteger(pi1.real(),-pi1.imag()),
							GaussianInteger(-pi1.imag(), pi1.real()),
							GaussianInteger(-pi1.real(),-pi1.imag()),
							GaussianInteger( pi1.imag(),-pi1.real()),
							GaussianInteger(pi1.real(),-pi1.imag()),
							GaussianInteger(pi1.imag(),pi1.real()),
							GaussianInteger(-pi1.real(),pi1.imag()),
							GaussianInteger(-pi1.imag(),-pi1.real())] :
		
		chi = Character(pi)
		chi2 = chi*chi
		#chi3 = chi*chi*chi
		#chi4 = chi*chi*chi*chi
		
		for a in xrange(pi.norm()) :
			a = GaussianInteger(a)
			#print a, chi(a)
			#print '%10s %10s %10s %10s %10s' % (a,chi(a),chi2(a),chi3(a),chi4(a))
		print '%10s %10s %10s %10s' %(pi,chi(GaussianInteger(-1)),JacobiSum(chi,chi,pi.norm()),JacobiSum(chi,chi2,pi.norm()))

#q2primes = [flip(J(qq)) for qq in xrange(5,maxp,4) if isprime(qq)]

#for pi in q2primes :
#	print pi
#	chi = Character(GaussianInteger(-1,2))
#	chi2 = chi*chi
#	print '%10s %10s' %(pi,JacobiSum(chi,chi,pi.norm()))
# print '%10s %10s %10s' %(pi,JacobiSum(chi,chi,pi.norm()),JacobiSum(chi,chi2,pi.norm()))

if False :
	print """A table of biquadratic residues.  Each column shows
	the values of the character associated with the prime at the
	top of the column
	"""
	print '%5s' % ' ',
	for pi2 in qprimes :
		print '%4s' % (GaussianInteger(pi2.real(),0),),
	print
	
	print '%5s' % ' ',
	for pi2 in qprimes :
		if pi2.imag()==0 :
			print '%4s' % (' '),
		else :
			print '%+4s' % (GaussianInteger(0,pi2.imag()),),
	print
	
	
	for pi1 in pprimes :
		print '%5s' % (pi1,),
		for pi2 in qprimes :
			if pi1 == pi2 :
				print '%4s' % ' ',
			else :
				print '%4s' % (biquad(pi2,pi1),),
		print


