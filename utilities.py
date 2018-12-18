#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

from math import sqrt
from gaussian_integer import GaussianInteger


def isprime(p):
    if type(p) != int:
        raise TypeError('%s is %s, but should be int' % (p, type(p)))

    if p < 0:
        p = -p

    if p == 0 or p == 1:
        return False

    trdiv = 2
    while trdiv * trdiv <= p:
        if p % trdiv == 0:
            return False
        trdiv += 1
    return True


def issq(nn):
    return nn >= 0 and nn == int(sqrt(nn)) ** 2


def factorize(n):
    """
    Return a non-decreasing list of primes whose product is n.
    """
    if n == 1:
        return []
    q = 2
    while q * q <= n:
        if n % q == 0:
            return [q] + factorize(n / q)
        q += 1
    return [n]


def squarefree(mm):
    """
    Return True if mm is square free, False otherwise.
    """
    factors = factorize(abs(mm))
    for i in range(len(factors) - 1):
        if factors[i] == factors[i + 1]:
            return False
    return True


def order(a, p):
    """
    Return the order of a in the multiplicative group (Z/pZ)^*.
    """
    one = type(a)(1)
    cnt = 1
    b = a
    while not b == one:
        b = (b * a) % p
        cnt += 1
    return cnt


def primitive_root(p):
    """
    Return a generator of the (cyclic) multiplicative group (Z/pZ)^*
    """
    for a in range(2, p):
        if order(a, p) == p - 1:
            return a


def modpow(a, k, p):
    """
    Return a**k (mod p).

    a can be of any type that has a multiplicative
    identity and supports multiplication and modding.
    """
    retval = type(a)(1)
    cnt = 0
    while cnt < k:
        retval = (retval * a) % p
        cnt += 1
    return retval


def modpow2(a, k, p):
    """
    Return a**k (mod p).
    """
    if k == 0:
        return type(a)(1)
    tail = (modpow2(a, k / 2, p) ** 2) % p
    if k % 2 == 1:
        return (a * tail) % p
    else:
        return tail % p


def legendre_ch(p):
    """
    Return the mod p legendre character, a function 'ch' such that

    ch(a) =  0  if p divides a
            +1  if a is a quadratic residue
            -1  if a is a quadratic non-residue
    """
    if not isprime(p) or p == 2:
        raise ValueError("%d is not an odd prime." % (p,))

    def ch(a):
        if a % p == 0:
            return 0
        rr = modpow2(a, (p - 1) / 2, p)
        return (-1) if (rr == p - 1) else +1

    return ch


def jacobi_sum(chi1, chi2, p):
    """
    Return the Jacobi sum of two mod p characters
    with values in the GaussianIntegers (so the
    characters are either the trivial character,
    the non-trivial quadratic character or one of
    the two non-trivial quartic characters.)
    """
    retval = GaussianInteger(0)
    for a in range(p):
        retval = retval + chi1(a) * chi2(1 - a)
    return retval


def jacobi_sum_quartic(p):
    """
    Returns the Jacobi sum J(chi,chi) associated with the
    quartic character which maps the smallest primitive
    root mod p to the pure imaginary unit i.
    """
    a = primitive_root(p)
    zeta = GaussianInteger(0, 1)
    avals = [0] + [modpow(a, k, p) for k in range(1, p)]
    achars = [GaussianInteger(0)] + [zeta ** (k % 4) for k in range(1, p)]
    bvals = [(1 - aa) % p for aa in avals]
    binds = [avals.index(bb) for bb in bvals]
    bchars = [achars[bi] for bi in binds]
    retval = GaussianInteger(0)
    for val in (ac * bc for (ac, bc) in zip(achars, bchars)):
        retval = retval + val
    return retval


def gcd(a, b):
    """
    Return the greatest common divisor of a and b.
    """
    if a == 0 or b == 0:
        return a + b
    while a % b != 0:
        a, b = b, a % b
    return b


def euclidean_algorithm(a, b):
    """
    Return x,y such that x*a + y*b = gcd(a,b).
    """
    q, r = divmod(a, b)  # a = q * b + r

    if r == 0:
        return 0, 1
    # d = gcd(a, b) = gcd(b, r)
    #
    # Suppose xb+yr = d
    #
    #  xb + y(a-qb) = d
    #
    #  ya + (x-yq)b = d

    x, y = euclidean_algorithm(b, r)
    return y, (x-y*q)


class UtilitiesTest(unittest.TestCase):

    def test_isprime(self):
        self.assertEqual(False, isprime(0), 'zero is not a (maximal) prime')
        self.assertEqual(False, isprime(1), 'units are not prime')
        self.assertEqual(True, isprime(7), 'yes, 7 is a prime')
        self.assertEqual(False, isprime(49), 'a non-unit square is not prime')
        self.assertEqual(False, isprime(91), '91=7*13 is not prime')
        self.assertEqual(True, isprime(-7), '(some) negative numbers are prime')
        self.assertRaises(TypeError, isprime, 7.0)

    def test_squarefree(self):
        self.assertEqual(True, squarefree(1), '1 is square free')
        self.assertEqual(True, squarefree(-1), '-1 is square free')
        self.assertEqual(False, squarefree(4), '4 is not square free')
        self.assertEqual(False, squarefree(18), '18 is not square free')

    def test_gcd(self):
        self.assertEqual(gcd(2*3*5, 3*5*7), 3*5)

    def test_euclidean_algorithm(self):
        a = 89
        b = 55
        x, y = euclidean_algorithm(a, b)
        self.assertEqual(abs(gcd(a, b)), abs(x*a + y*b))
