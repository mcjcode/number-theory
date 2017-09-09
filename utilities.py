import gaussian_integer
reload(gaussian_integer)
from gaussian_integer import GaussianInteger

def isprime(p) :
	trdiv = 2
	while trdiv*trdiv <= p :
		if p % trdiv == 0 :
			return False
		trdiv += 1
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
	
def modpow(a,k,p) :
	retval = type(a)(1)
	cnt = 0
	while cnt < k :
		retval = (retval * a) % p
		cnt += 1
	return retval

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
