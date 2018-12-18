#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

import numpy as np

from math import sqrt

from utilities import (
    squarefree,
    gcd,
    isprime,
    factorize,
    )


def discriminant(d):
    """
    Return the discriminant of the quadratic field Q[√d]
    """
    return d*(1 if d % 4 == 1 else 4)


def legendre(a, p):
    """
    Return (a/p), the Legendre symbol (p is prime)

    (a/p) = 0 if p|a
          = 1 if a is a quad residue mod p
          = -1 if a is a quad non-residue mod p
    """
    a = a % p
    if a == 0:
        return 0
    if a == 1:
        return 1
    if a == -1:
        if p % 4 == 3:
            return -1
        elif p % 4 == 1:
            return 1
        else:
            return 0
    if a == 2:
        p_mod_8 = p % 8
        if p_mod_8 == 1 or p_mod_8 == 7:
            return 1
        elif p_mod_8 == 3 or p_mod_8 == 5:
            return -1
        else:
            return 0
    if not isprime(a):
        return np.prod([legendre(q, p) for q in factorize(a)])

    if p % 4 == 1 or a % 4 == 1:
        return legendre(p % a, a)
    else:
        return -legendre(p % a, a)


def legendre_ch(d):
    """
    Return the mod |disc Q[√d]| Legendre character.

    This is the character ch with modulus D=|disc Q[√d]|
    such that for any odd prime p, ch(p)=(d/p) (i.e.

    ch(p)=0 if p|d
    ch(p)=1 if d is a square mod p
    ch(p)=-1 if d is not a square mod p
    """
    if not squarefree(d):
        raise ValueError('%d is not square free' % (d,))

    abs_disc = abs(discriminant(d))

    vals = [0]*abs_disc

    for k in range(abs_disc):
        if gcd(abs_disc, k) == 1:
            n = k
            while True:
                if isprime(n) and n % 2 == 1:
                    vals[k] = legendre(d, n)
                    break
                n += abs_disc
    return lambda a: vals[a % abs_disc]


def _ch12(a):
    """
    # the mod-12 character that identifies for which
    # primes p not dividing 12 is 3 a quadratic
    # residue.

    :param a:
    :return: 0/1/-1
    """
    # ----- 0.  1.  2.  3.  4.  5.  6.  7.  8.  9. 10  11
    return [0, +1,  0,  0,  0, -1,  0, -1,  0,  0,  0, +1][a % 12]


def minkowski_bound(d):
    """
    For an integer d let K=Q[√d].

    All ideal classes in O_K have a representative J
    with ||J|| <= the minkowski bound.
    """
    disc = abs(discriminant(d))
    if d > 0:
        return int(0.5 * sqrt(disc))
    else:
        return int(0.5 * (4.0/np.pi) * sqrt(disc))


class QuadInt(object):
    def __init__(self, d, a, b):
        self.d = d
        self.a = a
        self.b = b

    def __str__(self):
        if self.d % 4 != 1:
            part1 = u"%d" % (self.a,)
            sgn = "+" if self.b > 0 else "-"
            if abs(self.b) == 1:
                part2 = u"√%d" % (self.d,)
            else:
                part2 = u"%d√%d" % (abs(self.b), self.d)
            return u" ".join([part1, sgn, part2])
        else:
            if self.b == 1:
                if self.a == 0:
                    return u"(1+√%d)/2" % (self.d,)
                else:
                    return u"%d + (1+√%d)/2" % (self.a, self.d)
            else:
                return u"%d + %d(1+√%d)/2" % (self.a, self.b, self.d)

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        if self.d != other.d:
            raise ValueError("integers are not from same field")
        return QuadInt(self.d, self.a+other.a, self.b+other.b)

    def __mul__(self, other):
        if self.d != other.d:
            raise ValueError("integers are not from same field")
        else:
            return QuadInt(self.d, self.a*other.a-self.d*self.b*other.b, self.a*other.b-self.b*other.b)

    def real(self):
        if self.d % 4 == 1:
            return self.a + self.b * (1+sqrt(self.d))/2.0
        else:
            return self.a + self.b * sqrt(self.d)

    def norm(self):
        if self.d % 4 == 1:
            return self.a**2 + self.a * self.b + (1-self.d)/4 * self.b**2
        else:
            return self.a**2 - self.d * self.b**2


def factorize_in(p, d):
    """
    Return the factorization of (p) in Q[√d].

    If p remains prime in Q[√d], then return (p)
    otherwise return a prime (p,a) lying over (p).
    """
    if d % 4 != 1:
        # the minimal polynomial is x**2 - d = 0
        x = 0
        while x < p:
            if (x**2-d) % p == 0:
                return p, QuadInt(d, -x, 1)
            x += 1
        return p,
    else:
        # the minimal polynomial is x**2 - x - (d-1)/4
        zz = (d-1)/4
        x = 0
        while x < p:
            if (x*(x-1) - zz) % p == 0:
                return p, QuadInt(d, -x, 1)
            x += 1
        return p,


class LegendreCharacterTest(unittest.TestCase):
    def test_12(self):
        ch = legendre_ch(3)
        for xx in range(12):
            self.assertEqual(ch(xx), _ch12(xx), 'mod 3*4 Legendre character incorrect')
