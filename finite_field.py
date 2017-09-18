import numpy as np
import itertools
from utilities import factorize
from numpy import poly1d

def total_poly(p,n) :
	"""
	Return the polynomial p(x) that is satisfied
	by all elements of the finite field of degree
	p^n.  (i.e. x^(p^n)-x)
	"""
	return poly1d([1,0])**(p**n) - poly1d([1,0])

def product_of_irreducibles(d,p) :
	"""
	Return the product of all of the degree d 
	irreducible polynomials.
	"""
	prime_factors = list(set(factorize(d)))
	
	numer = total_poly(p,d)
	denom = poly1d([1])
	
	for ii in xrange(1,len(prime_factors)+1) :
		for subset in itertools.combinations(prime_factors,ii) :
			fact = total_poly(p,d/np.prod(subset))
			if ii % 2 == 0 :
				numer *= fact
			else :
				denom *= fact
	return (numer/denom)[0]
	
class FiniteFieldElement(object) :
	
	def __init__(self,value,ff) :
		self.value = value
		self.ff = ff
		
	def __repr__(self) :
		return 'FiniteFieldElement(%s,%s)' % (self.value,self.ff)
	
	def __repr__(self) :
		return '%s' % self.value

	def __add__(self,other) :
		return FiniteFieldElement((self.value+other.value) % self.ff.p, self.ff)
	
	def __iadd__(self,other) :
		self.value = (self.value + other.value) % self.ff.p
		return self
		
	def __sub__(self,other) :
		return FiniteFieldElement((self.value-other.value) % self.ff.p, self.ff)
	
	def __neg__(self) :
		return FiniteFieldElement(-self.value % self.ff.p, self.ff)
		
	def __pow__(self,xp) :
		retval = self.ff.one()
		runpow = self
		while xp :
			if xp & 1 :
				retval *= runpow
			runpow = runpow * runpow
			xp >>= 1
		return retval
		
	def __mul__(self,other) :
		return FiniteFieldElement(self.ff.mul(self.value,other.value), self.ff)
	
	def __imul__(self,other) :
		self.value = self.ff.mul(self.value,other.value)
		return self
		
	def __eq__(self,other) :
		return (self.ff == other.ff) and (self.value == other.value).all()

	def __ne__(self,other) :
		return not (self==other)

class FiniteField(object) :
	
	def __init__(self,p,n) :
		self.p = p
		self.n = n
		self.size = self.p**self.n
		self.rpoly = self.defining_polynomial()
	
	def polynomial_remainder_modp(self,dvd,dvs) :
		pc = poly1d(dvd.c)
		x = poly1d([1,0])
		while pc.order >= dvs.order :
			pc = pc - dvs*(x**(pc.order-dvs.order))*pc.c[0]
			pc = poly1d([coef % self.p for coef in pc.c])
		return pc
		
	def defining_polynomial(self) :
		indicators = np.zeros((self.p,)*(self.n),dtype=int)
		for ii in xrange(1,self.n//2+1) :
			for acoefs in itertools.product(*[xrange(self.p)]*(ii)) :
				for bcoefs in itertools.product(*[xrange(self.p)]*(self.n-ii)) :
					ccoefs = (poly1d((1,)+acoefs)*poly1d((1,)+bcoefs)).c
					ccoefs = tuple(xx % self.p for xx in ccoefs)
					#print ccoefs
					indicators[tuple(ccoefs[1:])] = 1
		#print indicators
		for kk in itertools.product(*[xrange(self.p)]*(self.n)) :
			if indicators[kk] == 0 :
				return poly1d((1,)+kk)
					
	def defining_polynomial_old(self) :
		##
		## find the product of all the monic irreducible
		## polynomials of degree d over Z/pZ.
		##
		pp = product_of_irreducibles(d=self.n,p=self.p)
		##
		## now loop over all degree d monic polynomials
		## and find a factor.  This will be guaranteed
		## to be irreducible.
		##
		for coefs in itertools.product(*[xrange(self.p)]*(self.n-1)) :
			coefs = [1] + list(coefs) + [1] # monic polynomial of degree self.n
			trial_poly = poly1d(coefs)
			r = self.polynomial_remainder_modp(pp,trial_poly)
			if all(r.c==0) :
				return trial_poly
							
	def __eq__(self,other) :
		
		return self.p == other.p and self.n == other.n and self.rpoly == other.rpoly
			
	def order(self,a) :
		cnt = 1
		b = a
		while not b==self.one() :
			if b==-self.one() :
				return 2*cnt
			b = b * a
			cnt += 1
		return cnt
	
	def primitive_element(self) :
		for elt in self :
			if elt == self.zero():
				continue
			if self.order(elt)==self.size-1:
				return elt
				
	def __repr__(self) :
		return 'FiniteField(%d,%d)' % (self.p,self.n)

	def __str__(self) :
		return self.__repr__()
	
	def mul(self,x,y) :
		z = poly1d(list(reversed(x)))*poly1d(list(reversed(y)))
		for ii in xrange(z.order+1) :
			z[ii] %= self.p
		z = self.polynomial_remainder_modp(z,self.rpoly)
		return np.array(list(reversed(z.c))+[0]*(self.n-len(z.c)))

	def zero(self) :
		return FiniteFieldElement(np.array([0]*self.n,dtype=int),self)
	
	def one(self) :
		return FiniteFieldElement(np.array([1]+[0]*(self.n-1),dtype=int),self)
		
	def __iter__(self) :
			return (FiniteFieldElement(np.array(val,dtype=int),self) for val in itertools.product(*[xrange(self.p)]*self.n))
			
	def J(self,ord) :
		"""
		Returns the Jacobi sum J(chi,chi) associated with the
		character which maps the smallest primitive
		root to the ord-root of unity
		"""
		N = self.p**self.n
		a = self.primitive_element()
		theta = 2*np.pi/ord
		zeta = np.cos(theta)+1j*np.sin(theta)
		
		chars = np.zeros((self.p,)*self.n,dtype=type(0j))
		chars[(0,)*self.n] = 0j
		pa = self.one()
		char_val = 1+0j
		for k in xrange(1,N) :
			chars[tuple(pa.value)] = char_val
			pa = pa * a
			char_val *= zeta
		
		charsum = 0j
		for aindex in itertools.product(*[xrange(self.p)]*self.n) :
			bindex = list((self.p-ai)%self.p for ai in aindex)
			bindex[0] = (bindex[0]+1)%self.p
			bindex = tuple(bindex)
			charsum += chars[aindex]*chars[bindex]

		c2 = charsum.imag/zeta.imag
		tmp = charsum - c2*zeta
		c1 = tmp.real
		return round(c1),round(c2)		

	
def count_projective_curve_points(f,ff) :
	
	N = 0
	zero = ff.zero()
	for xx in ff :
		for yy in ff :
			for zz in ff :
				if xx != zero or yy != zero or zz != zero :
					N += (f(xx,yy,zz) == zero)
	N //= ff.p**ff.n-1
	return N

def count_curve_points_affine(f,ff) :
	
	N = 0
	zero = ff.zero()
	for xx in ff :
		for yy in ff :
			if f(xx,yy) == zero :
				N += 1
	return N
	
def count_diagonal_cubic_points(ff) :
	N = 0
	zero = ff.zero()
	cubes = []
	for xx in ff :
		if xx**3 not in cubes :
			cubes.append(xx**3)
	#cubes = list(set(xx**3 for xx in ff))
	if not len(cubes) == (len(list(ff))-1)/3 + 1 :
		print cubes
	for xx in cubes :
		cx = (1 if xx==zero else 3)
		for yy in cubes :
			cy = (1 if yy==zero else 3)
			for zz in cubes :
				cz = (1 if zz==zero else 3)
				if xx+yy+zz == zero :
					N+=cx*cy*cz
	N -= 1 #(0,0,0) is not a projective solution
	assert N % (ff.p**ff.n-1) == 0
	N //= ff.p**ff.n-1
	return N
	
def assert_equal(expected, actual, msg) :
	if expected==actual :
		print '.',
		return True
	else :
		raise Exception('testing %s: %s expected, %s actual' % (msg, expected, actual))
		
if __name__ == '__main__' :
	import utilities; from utilities import isprime
	
	ff = FiniteField(7,2)
	L = list(ff)
	
	for xx in L :
		assert_equal(xx,ff.one()*xx,msg='1*X==X')
		assert_equal(xx,xx*ff.one(),msg='X*1==X')
		assert_equal(ff.zero()*xx,ff.zero(),msg='0*X=0')
		assert_equal(xx*ff.zero(),ff.zero(),msg='X*0=0')
		
	assert_equal(L[8]**2,L[8]*L[8],msg='X**2==X*X')
	
	print count_diagonal_cubic_points(ff)
	
