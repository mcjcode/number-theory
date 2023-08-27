#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import numpy as np
import math
from fractions import Fraction

import sympy
import sympy.abc

from sympy.polys.polytools import (
    compose,
)

from utilities import (
    prod,
    gcd,
    primitive_root,
    modpow,
    symmetric_function,
    factorize2,
)

from multiplicative import (
    phi,
    mu,
)

import jacobi

class CyclotomicInteger(object):

    @staticmethod
    def one(nn):
        return CyclotomicInteger([1] + [0]*(nn-1))

    @staticmethod
    def zero(nn):
        return CyclotomicInteger([0]*nn)

    @staticmethod
    def random(nn):
        coefs = np.random.randint(0, 20, (nn, )).tolist()
        return CyclotomicInteger(coefs)

    @staticmethod
    def zeta(nn):
        return CyclotomicInteger([0, 1] + [0]*(nn-2))

    @staticmethod
    def periods(nn, kk):
        """
        :param nn:
        :param kk: the order of the periods to return
        :return: a list of the order n gaussian periods
        """

        zero = CyclotomicInteger.zero(nn)
        one = CyclotomicInteger.one(nn)

        order = phi(nn)
        r = primitive_root(nn)
        powers = [modpow(r, k, nn) for k in range(1, order+1)]
        rts = [CyclotomicInteger.zeta(nn)**rr for rr in powers]

        if kk == 1:
            return rts

        arr = np.array(rts)
        arr.resize((kk, nn//kk))
        arr = arr.transpose().tolist()

        return list(map(sum, arr))

    def __init__(self, coefs):
        self.coefs = coefs

    def __repr__(self):
        return 'CyclotomicInteger(%s)' % (self.coefs, )

    def conjugates(self):
        """
        Yield the conjugates of the cyclotomic integer.
        """
        nc = len(self.coefs)
        for k in range(nc):
            if gcd(k, nc) == 1:
                yield CyclotomicInteger([self.coefs[ii*k % nc] for ii in range(nc)])

    def norm(self):
        """
        Return the product of the conjugates, as an integer.
        """
        xx = np.prod(list(self.conjugates()))
        assert len({yy for yy in xx.coefs[1:]}) == 1
        return xx.coefs[0] - xx.coefs[1]

    def __eq__(self, other):
        return len({xx - yy for (xx, yy) in zip(self.coefs, other.coefs)}) == 1

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        if type(other) == int:
            other = CyclotomicInteger([other]+[0]*(len(self.coefs)-1))
        return CyclotomicInteger([xx+yy for (xx, yy) in zip(self.coefs, other.coefs)])

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return CyclotomicInteger([xx-yy for (xx, yy) in zip(self.coefs, other.coefs)])

    def __neg__(self):
        return CyclotomicInteger([-xx for xx in self.coefs])

    def __mul__(self, other):
        if type(other) == int:
            other = CyclotomicInteger([other]+[0]*(len(self.coefs)-1))
        if len(self.coefs) != len(other.coefs):
            raise ValueError('values are from different cyclotomic fields')

        nc = len(self.coefs)
        ret_coefs = [0]*nc
        for ii in range(nc):
            ret_coefs[ii] = sum([self.coefs[jj]*other.coefs[(ii-jj) % nc] for jj in range(nc)])
        return CyclotomicInteger(ret_coefs)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, n):
        retval = CyclotomicInteger([1]+[0]*(len(self.coefs)-1))
        cnt = 0
        while cnt < n:
            retval = retval * self
            cnt += 1
        return retval


def cyclotomic_polynomial(nn, kk=1):
    """
    Return the defining polynomial for the degree phi(nn)/kk
    subfield of Q[zeta_nn].

    kk must divide phi(nn) and if kk!=1, only works for nn prime.

    For kk=1, we just use the moebius inversion for Phi_n(x),

    For kk>1, we 'program by algebra' (i.e. spelling out the Gaussian
    periods and computing the elementary symmetric functions of them)
    and is very very slow for large nn (and small kk).
    """
    x = sympy.abc.x
    if kk==1:
        ans = sympy.poly(x-1, x, domain='ZZ')
        for p, e in factorize2(nn):
            ans = compose(ans, x**(p**e)) // compose(ans, x**(p**(e-1)))
        return ans
    else:
        pds = CyclotomicInteger.periods(nn, kk)
        coefs_ci = [(-1)**dd * symmetric_function(dd, pds) for dd in range(1, phi(nn)//kk + 1)]
        coefs = [ci.coefs[0]-ci.coefs[1] for ci in coefs_ci]
        return sympy.Poly([1]+coefs, x)
        

#
# We implement the algorithm of R.P.Brent (On Computing Factors of Cyclotomic Polynomials)
# for computing the Gauss and Lucas identities.
#

#
# put my own version here representing elements of Q[sqrt{d}]
# as Q-linear combinations of 1 and sqrt(d) (not (1+sqrt(d))/2)
#
class QuadInt(object):
    def __init__(self, d, a, b):
        self.d = d
        self.a = a
        self.b = b

    def __str__(self):
        part1 = u"%s" % (self.a, )
        sgn = "+" if self.b > 0 else "-"
        if abs(self.b) == 1:
            part2 = u"sqrt(%d)" % (self.d, )
        else:
            part2 = u"%ssqrt(%d)" % (abs(self.b), self.d)
        return u" ".join([part1, sgn, part2])

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        if type(other)==int or type(other)==Fraction:
            return QuadInt(self.d, self.a+other, self.b)
        
        if self.d != other.d:
            raise ValueError("integers are not from same field")
        return QuadInt(self.d, self.a+other.a, self.b+other.b)

    def __radd__(self, other):
        if type(other)==int or type(other)==Fraction:
            return QuadInt(self.d, self.a+other, self.b)
        else:
            return self.__add__(other)
        
    def __sub__(self, other):
        if self.d != other.d:
            raise ValueError("integers are not from same field")
        return QuadInt(self.d, self.a-other.a, self.b-other.b)
    
    def __mul__(self, other):
        if type(other)==int or type(other)==Fraction:
            return QuadInt(self.d, self.a*other, self.b*other)
        
        if self.d != other.d:
            raise ValueError("integers are not from same field")
        else:
            return QuadInt(self.d, self.a*other.a+self.d*self.b*other.b, self.a*other.b+self.b*other.a)

    def __truediv__(self, other):
        if type(other)==int or type(other)==Fraction:
            return QuadInt(self.d, self.a/other, self.b/other)

    def __rmul__(self, other):
        if type(other)==int or type(other)==Fraction:
            return QuadInt(self.d, self.a*other, self.b*other)
        else:
            return self.__mul__(other)
        

#
# put my own version here that uses 'true' division
# instead of integer division (z/j and not z//j).
#
def newtons_identity(ps):
    es = [1]
    for j in range(1, len(ps)):
        z = 0
        sign = 1
        for i in range(1, j+1):
            z = z + sign*ps[i]*es[j-i]
            sign *= -1
        es.append(z/j)
    return es


def gauss_formula(n):
    """
    Return polynomials A(x) and B(x) such that A**2 +/- B**2
    equals the cyclotomic polynomial Phi_n(x)
    """
    def p(n):
        s = 2 - n%4
        yield QuadInt(s*n, phi(n)//2, 0)
        for k in range(1, phi(n)//2+1):
            g = math.gcd(k, n)
            if g==1:
                yield Fraction(mu(n),2) + jacobi.jacobi2(k, n)*QuadInt(s*n, Fraction(0,2), Fraction(1,2))
            else:
                yield QuadInt(s*n, Fraction(mu(n//g)*phi(g), 2), Fraction(0, 2))

    x = sympy.abc.x
    
    xs = list(p(n))
    ys = newtons_identity(xs)
    ys = [(-1)**i*y for i, y in enumerate(ys)]
    
    Acoeffs = [2] + [(2*z.a).numerator for z in ys[1:]]
    Bcoeffs = [0] + [(2*z.b).numerator for z in ys[1:]]
    
    A = sympy.Poly(Acoeffs, x, domain='ZZ')
    B = sympy.Poly(Bcoeffs, x, domain='ZZ')
    
    return A, B

def lucas_formula(n):
    """
    Return polynomials C(x) and D(x) such that C**2 +/- n*x*D**2
    equals the cyclotomic polynomial Phi_n(x)
    """
    if not(n%2 and all(e==1 for _, e in factorize2(n))): # has to be odd and squarefree
        raise ValueError(f'n has to be odd and squarefree')

    def p(n):
        yield QuadInt(n, phi(n), 0)
        sprime = -1 if n%8==5 else +1
        for k in range(1, phi(n)+1):
            np = n if n%4==1 else 2*n
            gp = math.gcd(k, np)
            if k%2==1:
                yield QuadInt(n, 0, sprime*jacobi.jacobi2(n, k))
            else:
                if ((n-1)//2 * k)%4==0:
                    z = 1
                elif ((n-1)//2 * k)%2==0:
                    z = -1
                else:
                    z = 0
                #z = math.cos((n-1)*k*math.pi/4)
                yield QuadInt(n, mu(np/gp)*phi(gp)*int(z), 0)

    xs = list(p(n))
    ys = newtons_identity(xs)
    ys[0] = QuadInt(n, 1, 0)
    Ccoeffs = [int(z.a) for z in ys[0::2]]
    Dcoeffs = [int(-z.b) for z in ys[1::2]]
        
    x = sympy.abc.x
    C = sympy.Poly(Ccoeffs, x, domain='ZZ')
    D = sympy.Poly(Dcoeffs, x, domain='ZZ')

    return C, D
