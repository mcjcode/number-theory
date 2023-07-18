#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import numpy as np
import sympy
import sympy.abc

from sympy.polys.polytools import (
    compose,
)

from utilities import (
    gcd,
    primitive_root,
    modpow,
    symmetric_function,
    factorize2,
)

from multiplicative import (
    phi,
)


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
        return ans ##.coeffs()
    else:
        pds = CyclotomicInteger.periods(nn, kk)
        print('pds = ', pds)
        coefs_ci = [(-1)**dd * symmetric_function(dd, pds) for dd in range(1, phi(nn)//kk + 1)]
        print('coefs_ci = ', coefs_ci)
        coefs = [ci.coefs[0]-ci.coefs[1] for ci in coefs_ci]
        return sympy.Poly([1]+coefs, x)



        
