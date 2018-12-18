#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

import numpy as np

import itertools

from utilities import (
    gcd,
)


class CyclotomicInteger(object):

    @staticmethod
    def one(nn):
        return CyclotomicInteger([1] + [0]*(nn-1))

    @staticmethod
    def zeta(nn):
        return CyclotomicInteger([0, 1] + [0]*(nn-2))

    def __init__(self, coefs):
        self.coefs = coefs

    def __repr__(self):
        return 'CyclotomicInteger(%s)' % (self.coefs,)

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
        return CyclotomicInteger([xx+yy for (xx, yy) in zip(self.coefs, other.coefs)])

    def __sub__(self, other):
        return CyclotomicInteger([xx-yy for (xx, yy) in zip(self.coefs, other.coefs)])

    def __neg__(self):
        return CyclotomicInteger([-xx for xx in self.coefs])

    def __mul__(self, other):
        if len(self.coefs) != len(other.coefs) :
            raise ValueError('values are from different cyclotomic fields')

        nc = len(self.coefs)
        ret_coefs = [0]*nc
        for ii in range(nc):
            ret_coefs[ii] = sum([self.coefs[jj]*other.coefs[(ii-jj) % nc] for jj in range(nc)])
        return CyclotomicInteger(ret_coefs)

    def __pow__(self, n):
        retval = CyclotomicInteger([1]+[0]*(len(self.coefs)-1))
        cnt = 0
        while cnt < n:
            retval = retval * self
            cnt += 1
        return retval


class CyclotomicIntegerTest(unittest.TestCase):

    def test_1(self):
        nn = 13
        one = CyclotomicInteger.one(nn)
        zeta = CyclotomicInteger.zeta(nn)
        self.assertEqual(one * zeta, zeta)
        self.assertEqual(zeta ** nn, one)

    def test_norm(self):
        nn = 13
        xx = CyclotomicInteger.one(nn) - CyclotomicInteger.zeta(nn)
        self.assertEquals(xx.norm(), nn)
