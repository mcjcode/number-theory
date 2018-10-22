import math
from math import sqrt
import utilities
from utilities import gcd

def form_disc(a,b,c) :
	return b**2 - 4*a*c


def sgn(b) :
	if b > 0 :
		return 1
	elif b == 0 :
		return 0
	else :
		return -1
		
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
			r = -sgn(b)*int(abs(b) / (2.0*a) + 0.5)
			b, c = b + 2*a*r, a*r*r + b*r + c
		elif abs(b) > c :
			r = -sgn(b)*int(abs(b) / (2.0*c) + 0.5)
			a, b = c*r*r + b*r + a, b + 2*c*r
	if c < a :
		a, b, c = c, -b, a
	if a == -b :
		r = 1
		b, c = b + 2*a, c + a + b
		
	return a, b, c
	
def all_reduced_forms(D) :
	retval = []
	for a in xrange(1,int(sqrt(-D/3.0)+1.0)) :
		for b in xrange(-a,a+1) :
			# D = b**2 - 4*a*c
			if (b**2 - D) % (4*a) == 0 :
				c = (b**2 - D) / (4*a)
				if gcd(gcd(a,abs(b)),c) == 1 :
					retval.append(proper_reduced_form(a,b,c))
	return sorted(list(set(retval)))
