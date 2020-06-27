#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

import numpy as np

from utilities import (
    gcd,
    primitive_root,
    modpow,
    symmetric_function,
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

        return map(sum, arr)

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


def cyclotomic_polynomial(nn, kk):
    """
    Return the defining polynomial for the degree phi(nn)/kk
    subfield of Q[zeta_nn].

    Only works for nn prime for now.

    Note that this implementation is 'programming by algebra'
    (i.e. spelling out the Gaussian periods and computing the
    elementary symmetric functions of them) and is very very
    slow for large nn (and small kk)
    """
    pds = CyclotomicInteger.periods(nn, kk)
    coefs_ci = [(-1)**dd * symmetric_function(dd, pds) for dd in range(1, phi(nn)//kk + 1)]
    coefs = [ci.coefs[0]-ci.coefs[1] for ci in coefs_ci]
    return [1]+coefs


class CyclotomicIntegerTest(unittest.TestCase):

    def test_axioms(self):
        nn = 13
        zero = CyclotomicInteger.zero(nn)
        one = CyclotomicInteger.one(nn)

        xx = CyclotomicInteger.random(nn)
        yy = CyclotomicInteger.random(nn)
        zz = CyclotomicInteger.random(nn)

        self.assertEqual(one * xx, xx)
        self.assertEqual(zero + xx, xx)
        self.assertEqual(xx + yy, yy + xx)
        self.assertEqual(xx * yy, yy * xx)
        self.assertEqual(xx * (yy + zz), xx * yy + xx * zz)

    def test_norm(self):
        nn = 13

        one = CyclotomicInteger.one(nn)
        zeta = CyclotomicInteger.zeta(nn)

        xx = one - zeta
        yy = one + zeta

        self.assertEqual(zeta ** nn, one)
        self.assertEqual(xx.norm(), nn)
        self.assertEqual(yy.norm(), 1)
