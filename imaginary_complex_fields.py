#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import math
from math import sqrt
import numpy as np
import utilities
from utilities import (
	gcd,
	issq,
	isprime,
	squarefree,
	)
from itertools import islice
import quadratic_extensions
from quadratic_extensions import (
	QuadInt, 
	factorize_in,
	minkowski_bound,
	)

ordinar_chars_lower = "`1234567890-=qwertyuiop[]\asdfghjkl;'zxcvbnm,./"
special_chars_lower = "`¡™£¢∞§¶•ªº–≠œ∑´®†¥¨ˆøπ“‘«åß∂ƒ©˙∆˚¬…æΩ≈ç√∫˜µ≤≥÷"

def form_disc(a,b,c) :
	return b**2 - 4*a*c

def genus(a,b,c) :
	f = lambda x,y : a*x**2 + b*x*y + c*y**2
	D = form_disc(a,b,c)
	val = sorted(list(set([f(x,y)%(-D) for x in xrange(-D) for y in xrange(-D) ])))
	return filter(lambda x : gcd(x,-D)==1, val)

def proper_reduced_form(a,b,c) :
	"""
	Return a reduced form properly equivalent to
	
	ax^2 + bxy + cy^2
	"""
	
	## first check that the discriminant is negative
	disc = form_disc(a,b,c)
	if disc >= 0 :
		raise ValueError("discriminant is non-negative")
	if a < 0 or c < 0 :
		raise ValueError("form is not positive_definite")
		
	while abs(b) > a or abs(b) > c :
		if abs(b) > a :
			r = -np.sign(b)*int(abs(b) / (2.0*a) + 0.5)
			b, c = b + 2*a*r, a*r*r + b*r + c
		elif abs(b) > c :
			r = -np.sign(b)*int(abs(b) / (2.0*c) + 0.5)
			a, b = c*r*r + b*r + a, b + 2*c*r
	if c < a :
		a, b, c = c, -b, a
	if a == -b :
		r = 1
		b, c = b + 2*a, c + a + b
	if a==c :
		b = abs(b)
	return a, b, c

def all_reduced_forms(D) :
	"""
	Return a list of all reduced forms with
	negative discriminant D.
	"""
	if D >= 0 :
		raise ValueError("discriminant must be negative")
	retval = []
	for a in xrange(1,int(sqrt(-D/3.0)+1.0)) :
		for b in xrange(-a+1,a+1) :
			# D = b**2 - 4*a*c
			c, r = divmod(b**2-D, 4*a)
			if r == 0 :
				if gcd(gcd(a,abs(b)),c) == 1 :
					##retval.append(proper_reduced_form(a,b,c))
					if (a <= c) and ((a != c) or (b>=0)) :
						retval.append((a,b,c))
	return sorted(list(set(retval)))

def pids() :
	D = -1
	while True :
		nforms = len(all_reduced_forms(D))
		if nforms == 1 :
			if D % 4 == 0 :
				yield -D // 4
			else :
				yield (1-D) // 4
		D -= 1

def idoneal(N) :
	for D in xrange(-4,-N,-4) :
		forms = all_reduced_forms(D)
		nforms = len(forms)
		ngenera = len(list(set([genus(*form)[0] for form in forms])))
		if nforms == ngenera :
			print D,

def intsqrt(k) :
	return int(sqrt(k))
	
def principal_representations(n,m) :
	"""
	yield a sequence of positive (a,b) such that
	a**2 + n * b**2 = m
	"""
	ub = intsqrt( m / n )
	ii = 0
	while ii <= ub :
		xx = m - n*ii**2
		if issq(xx) :
			yield (intsqrt(xx),ii)
		ii += 1

def repmod11(p) :
	princ_reps = list(principal_representations(11,p))
	if len(princ_reps) > 0 :
		return u"%d\u00B2 + 11\u00B7%d\u00B2" % princ_reps[0]
	else :
		princ_reps = list(principal_representations(11,3*p))
		(a,b) = princ_reps[0]
		if (a-b) % 3 == 0 :
			c = (a-b)/3
			if c > 0 : 
				return u"3\u00B7%d\u00B2 + 2\u00B7%d\u00B7%d + 4\u00B7%d\u00B2" % (c,c,b,b)
			else :
				c = -c
				return u"3\u00B7%d\u00B2 - 2\u00B7%d\u00B7%d + 4\u00B7%d\u00B2" % (c,c,b,b)				
		else :
			c = (a+b)/3
			if c > 0 :
				return u"3\u00B7%d\u00B2 - 2\u00B7%d\u00B7%d + 4\u00B7%d\u00B2" % (c,c,b,b)
			else :
				c = -c
				return u"3\u00B7%d\u00B2 + 2\u00B7%d\u00B7%d + 4\u00B7%d\u00B2" % (c,c,b,b)
										
def repsmod11() :
	p = 1
	gen = genus(1,0,11)
	while True :
		if (p % 44) in gen and isprime(p) :
			yield p, repmod11(p)
		p += 1

def makelist() :
	for xx in list(islice(repsmod11(),100)) :
		print "%5d =" % (xx[0],),
		print xx[1]
		
def foo() :
	for D in range(-4,-8000,-4) :
		forms = all_reduced_forms(D)
		if all( [(a==b or a==c or b==0) for (a,b,c) in forms]) :
		##if all([form[1] == 0 for form in forms]) :
		##if (len(forms) in [1,2,4,8,16,32,64,128,256]) :
			#for form in all_reduced_forms(D) :
			#	print D, #form#, genus(*form)
			yield -D / 4
			
def class_groups_of_size(k) :
	for D in range(-4,-5000,-4) :
		forms = all_reduced_forms(D)
		##if all( [(a==b or a==c or b==0) for (a,b,c) in forms]) :
		##if all([form[1] == 0 for form in forms]) :
		##if (len(forms) in [1,2,4,8,16,32,64,128,256]) :
		if len(forms) == k :
			for form in forms :
				print D, form, genus(*form)
			print
		
def ideal_class_group_info(d) :
	"""
	Calculate info relevant to the class group of Q[√-d].
	"""
	mb = minkowski_bound(d)
	
	print("Minkowski bound.  All ideal classes contain an ideal of norm ≤ %d." % (mb,))
	
	split_primes = []
	for p in filter(isprime,range(2,mb+1)) :
		fact = factorize_in(p,-d)
		print(fact)

