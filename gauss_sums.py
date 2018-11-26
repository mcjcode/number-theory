#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import numpy as np

from itertools import islice

import utilities
from utilities import (
	primitive_root,
	isprime,
	modpow,
	gcd,
	squarefree,
	)
import quadratic_extensions
from quadratic_extensions import legendre_ch

import real_quadratic_fields
from real_quadratic_fields import (
	approximants, 
	fundamental_unit,
	fundamental_unit_old,
	discriminant,
	)

def gs_numerical(ch,m,k) :
	"""
	Compute the Gauss sum 'sum_a(ch(a)w**(ak))'
	of the mod m characher ch.
	"""
	ww = np.exp(2.0j*np.pi/m)
	return sum (ch(a)*ww**(k*a) for a in xrange(1,m))

def L_one_chi(ch,m) :
	"""
	Compute the value of the Dirichlet L-series L(s,ch) at s=1.
	"""
	zz = np.exp(2.0j*np.pi/m)
	return -(1.0/m) * sum( gs_numerical(ch,m,k)*np.log(1-zz**(-k)) for k in xrange(1,m) if gcd(k,m)==1 )

def all_modm_chars(m) :
	"""
	Returns a list of all mod m characters,
	starting with the trivial character.  Only
	works with prime m for now.
	"""
	if not isprime(m) :
		raise ValueError("%d is not a prime." % (m,))
	a = primitive_root(m)
	zz = np.exp(2.0j*np.pi/(m-1))
	
	vals = [None]*(m-1)
	pows = [0] + [ modpow(a,k,m) for k in xrange(m-1) ]
	for k in xrange(m-1) :
		pows2 = [0] + [ (zz**k)**i for i in xrange(m-1) ]
		vals[k] = dict(zip(pows,pows2))
	##return vals
	return map(lambda val : (lambda i:val[i%m]), vals)

## Calculate the value of kappa for the QQ[z_5], the cyclotomic
## field of 5th roots of unity.  This has 4 complex embeddings,
## no real embeddings, (s=2,r=0,n=r+s=2), (1+z_5) is a fundamental
## unit, the regulator = (1/2)log(||(1+z_5)||), the discriminant
## of the QQ[z_5] is 5^(5-2) (i.e. disc(QQ[z_p])=p^(p-2)), and there
## are 5 roots of unity in the field.
##
##>>> (2.0*np.pi)**2.0*(0.5)*np.log((1+zz)*(1+zz).conj())/(5.0*np.sqrt(abs(5**3)))
##(0.33983727824052351+0j)
##
## now compare to the product of the L(1,ch) as ch runs over
## all non-trivial mod 5 characters:
##
##>>> np.product([L_one_chi(ch,5) for ch in all_modm_chars(5)[1:]])
##(0.33983727824052373+0j)

def kappa(m) :
	""" 
	Return kappa = lim_{t->infty} #{J | ||J|| < t}/t for Q[√m].
	
	|L(chi,1)| = h*kappa, where h is the class number of Q[√m],
	chi is the 'legendre character' and L is the associated
	Dirichlet L function.
	""" 
	if not squarefree(m) :
		raise ValueError('%d is not square free.' % (m,))
	
	if m == 1 :
		return 1
		
	D = discriminant(m)
	if m > 0 :
		u = fundamental_unit(m).real()
		return 2*np.log(u)/np.sqrt(abs(D))
	else :
		if m == -1 :
			w = 2
		elif m == -3 :
			w = 3
		else :
			w = 1
		return np.pi / (w*np.sqrt(abs(D)))
		
def ideal_class_number(d) :
	"""
	Return the ideal class number of Q[√d], for squarefree integer d.
	"""
	if not squarefree(d) :
		raise ValueError("%d is not squarefree." % (d,))

	if d==1 :
		return 1
		
	D = abs(discriminant(d))
		
	ch = legendre_ch(d)
	
	z=np.exp(2.0j*np.pi/D)
	rho = np.abs(1.0/np.sqrt(abs(D)) * sum(
		ch(k) * np.log(1-z**(-k)) for k in xrange(1,D) if gcd(k,D)==1)
		)
	
	k = kappa(d)
	
	# magically rho should be an integral multiple of kappa.
	# make sure this is true before we return the rounded answer.
	assert( abs(round(rho/k) - rho/k) < 0.001 )
	
	return int(round(rho/k))

def test_ideal_class_number() :
	for d in range(3000) :
		if squarefree(d) and d != 0 :
			print '%4d %3d ' % (d, ideal_class_number(d))

