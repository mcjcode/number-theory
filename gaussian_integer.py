class GaussianInteger(object) :
	
	def __init__(self,a,b=0) :
		self.a = a
		self.b = b

	def __str__(self) :
		if self.b == 0 :
			return '%d' % self.a
		elif self.a == 0 :
			if self.b == 1 :
				return 'i'
			elif self.b == -1 :
				return '-i'
			return '%di' % self.b
		else :
			if self.b == 1 :
				return '%d+i' % (self.a,)
			else :
				return '%d%+di' % (self.a,self.b)
	
	def __repr__(self) :
		return 'GaussianInteger(%d,%d)' % (self.a,self.b)
	
	def real(self) :
		return self.a
		
	def imag(self) :
		return self.b
	
	def norm(self) :
		return (self * self.conj()).real()
		
	def __eq__(self,other) :
		return self.a == other.a and self.b == other.b
		
	def __ne__(self,other) :
		return not (self==other)
		
	def __add__(self,other) :
		return GaussianInteger(self.a+other.a,self.b+other.b)
	
	def __sub__(self,other) :
		return GaussianInteger(self.a-other.a,self.b-other.b)
		
	def __neg__(self) :
		return GaussianInteger(-self.a,-self.b)
		
	def __mul__(self,other) :
		return GaussianInteger(self.a*other.a - self.b*other.b, self.a*other.b+self.b*other.a)
	
	def __pow__(self,n) :
		retval = GaussianInteger(1,0)
		cnt = 0
		while cnt < n :
			retval = retval * self
			cnt += 1
		return retval
		
	def __div__(self,other) :
		N = other.norm()
		numer = self * other.conj()
		return GaussianInteger(numer.real()/N, numer.imag()/N)
	
	def conj(self) :
		return GaussianInteger(self.a,-self.b)
		
	def __mod__(self,other) :
		"""
		Return the unique gaussian integer congruent
		to self mod other that is in the fundamental
		square of C containing the number 1
		"""
		
		if other == GaussianInteger(0) :
			return GaussianInteger(self.a,self.b)
		
		if other.a < 0 :
			other = other * GaussianInteger(-1,0)
		if other.b < 0 :
			other = other * GaussianInteger(0,1)
		
		perp = other * GaussianInteger(0,-1)
		
		## Now other, perp form an orthogonal basis
		## in the positive (counterclockwise)
		## orientation.
		
		N = other.norm()
		u = (self * other.conj()).real() / N
		v = (self * perp.conj()).real() / N
		
		div = self / other
		u1 = div.real()
		u2 = div.imag()
		
		retval = self - GaussianInteger(u,0) * other - GaussianInteger(v,0) * perp
		retval2 = self - div * other
		
		return retval2

