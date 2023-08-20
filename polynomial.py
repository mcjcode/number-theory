class intpoly1d():
    def __init__(self,c):
        self.coefs = c
        
    def __add__(self, other):
        c = [0]*max(len(self.coefs),len(other.coefs))
        for i in range(len(self.coefs)):
            c[i] += self.coefs[i]
        for i in range(len(other.coefs)):
            c[i] += other.coefs[i]
        return intpoly1d(c)
    
    def __sub__(self, other):
        c = [0]*max(len(self.coefs),len(other.coefs))
        for i in range(len(self.coefs)):
            c[i] += self.coefs[i]
        for i in range(len(other.coefs)):
            c[i] -= other.coefs[i]
        return intpoly1d(c)
    
    def __rmul__(self, scalar):
        deg = len(self.coefs)
        c = [0]*deg
        for i in range(deg):
            c[i] = self.coefs[i] * scalar
        return intpoly1d(c)
    
    def __call__(self, n):
        deg = len(self.coefs)
        a = 1
        i = 0
        retval = 0
        while i < deg:
            retval += self.coefs[i]*a
            a *= n
            i += 1
        return retval
    
    def __repr__(self):
        return f'intpoly1d({self.coefs})'