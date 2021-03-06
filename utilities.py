#!/usr/bin/env python -i
# -*- coding: utf-8 -*-
"""
General purpose, factorization and modular
arithmetic routines.
"""

import time
import unittest
import itertools

from math import sqrt


def prod(xs, start=1):
    """
    :param xs: a sequence of elements to multiply
    :param start: the value to return if xs is empty
    :return: the product of the xs
    """
    retval = start
    for x in xs:
        retval *= x
    return retval


def timeit(f):
    """
    :param f: a function
    :return: a function with the side effect of printing the time it takes to run f

    A decorator to wrap around functions when we
    want to measure the run time.
    """
    def g(*args):
        t0 = time.time()
        retval = f(*args)
        t1 = time.time()
        print('%.3f seconds' % (t1-t0))
        return retval
    return g


def memoize(f):
    """
    :param f: a function
    :return: a memoized version of f.
    """
    fhash = {}

    def g(*args):
        targs = tuple(args)
        if targs in fhash:
            retval = fhash[targs]
        else:
            retval = f(*targs)
            fhash[targs] = retval
        return retval
    return g


def mymod(n, m):
    """
    :param n: an integer
    :param m: an integer, the modulus
    :return: :math:`n(mod m)`, if m != 0, otherwise n.

    (Kind of what I'd want '%' to do in the first place.)
    """
    return n % m if m else n


def sqrtInt(n):
    """
    :param n: a non-negative integer
    :return: the largest integer smaller than the square root of n

    Note that this depends on `math.sqrt`, and its double precision
    accuracy means that this function should not be trusted for n on
    the order of :math:`10^{52}` and up.
    """
    sqrtn = int(sqrt(n))
    if (sqrtn+1)**2 <= n:
        sqrtn += 1
    return sqrtn


def cbrtInt(n):
    """
    :param n: a non-negative integer
    :return: the largest integer smaller than the cube root of n

    Note that this depends on `math.sqrt`, and its double precision
    accuracy means that this function should not be trusted for n on
    the order of :math:`10^{52}` and up.
    """
    cbrtn = int(n**(1./3.))
    if (cbrtn+1)**2 <= n:
        cbrtn += 1
    return cbrtn


def multiplicities(xs):
    """
    :param xs: a list of elements
    :return: the list of unique items in xs and how often they appear in xs
    """
    items = sorted(list(set(xs)))
    counts = []
    for item in items:
        counts.append(len([xx for xx in xs if xx == item]))
    return items, counts


def symmetric_function(k, xs, zero=0, one=1):
    """
    :param k: the degree of the symmetric function
    :param xs: the list of elements
    :param zero: the zero of the elements
    :param one: the 'one' of the elements
    :return: the value of the kth symmetric function evaluated on xs
    """
    retval = zero
    for comb in itertools.combinations(xs, k):
        term = one
        for elt in comb:
            term = term * elt
        retval = retval + term
    return retval


def isprime(p):
    """
    :param p: an integer
    :return: a boolean saying whether p is a prime number
    """
    if type(p) != int:
        raise TypeError('%s is %s, but should be int' % (p, type(p)))

    if p < 0:
        p = -p

    if p == 0 or p == 1:
        return False

    sqrtp = sqrtInt(p)
    trdiv = 2
    while trdiv <= sqrtp:
        if p % trdiv == 0:
            return False
        trdiv += 1
    return True


def issq(nn):
    """
    :param nn: an integer
    :return: a boolean telling whether nn is a square
    """
    return nn >= 0 and nn == int(sqrt(nn)) ** 2


def factorize(n):
    """
    :param n: a positive integer
    :return: a non-decreasing list of primes whose product is n.
    """
    if n == 1:
        return []
    q = 2
    while q * q <= n:
        if n % q == 0:
            return [q] + factorize(n//q)
        q += 1
    return [n]


def factorize2(n):
    """
    :param n: a positive integer
    :return: yields a sequence of (prime, exponent) pairs where the
             primes are distinct and in increasing order giving
             the prime factorization of n.
    """
    if n < _maxn:  # see below
        g = factorize2_bounded(n)
        for fact in g:
            yield fact
        return

    q = 2
    while q*q <= n:
        e = 0
        while n % q == 0:
            e += 1
            n = n//q
        if e > 0:
            yield q, e
        q += 1
    if n > 1:
        yield n, 1


_maxn = 10**10
_ps = [p for p in range(2, sqrtInt(_maxn)) if isprime(p)]


def factorize2_bounded(n):
    """
    Yield a sequence of (prime, exponent) pairs where the
    primes are distinct and in increasing order giving
    the prime factorization of n.

    prime factors are precomputed.  Only n < _maxn are
    allowed. 
    """
    assert n < _maxn
    for p in _ps:
        if p*p > n:
            break
        e = 0
        while n % p == 0:
            e += 1
            n //= p
        if e:
            yield p, e
    if n > 1:
        yield n, 1


def squarefree(mm):
    """
    :param mm: a positive integer
    :return: True if mm is square free, False otherwise.
    """
    factors = factorize(abs(mm))
    for i in range(len(factors) - 1):
        if factors[i] == factors[i + 1]:
            return False
    return True


def order(a, p):
    """

    :param a: an integer relatively prime to p
    :param p: a positive integer
    :return: the order of a in the multiplicative group :math:`(Z/pZ)^*`.
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
    :param p: a prime number
    :return: a generator of the (cyclic) multiplicative group :math:`(Z/pZ)^*`.
    """
    facts = list(factorize2(p-1))
    a = 2
    while a < p:
        ok = True
        for (q, e) in facts:
            if modpow2(a, (p-1)//q, p) == 1:
                ok = False
                break
        if ok:
            return a
        a += 1


def modpow(a, k, p):
    """
    :param a: the base
    :param k: the exponent (a positive integer)
    :param p: the modulus
    :return: :math:`a^k(p)`.

    a can be of any type that has a multiplicative
    identity and supports multiplication and modding.
    """
    retval = type(a)(1)
    cnt = 0
    while cnt < k:
        retval = (retval * a) % p
        cnt += 1
    return retval


def powerset(xs):
    """
    :param xs: a list of elements
    :return: a generator iterating over all of subsets of
             xs, starting with the smallest and ending with
             the largest subsets
    """
    lengths = range(len(xs)+1)
    return itertools.chain(*[itertools.combinations(xs, nn) for nn in lengths])


def modpow2(a, k, p):
    """
    :param a: the base
    :param k: the exponent (a positive integer)
    :param p: the modulus
    :return: :math:`a^k(p)`.

    a can be of any type that has a multiplicative
    identity and supports multiplication and modding.
    O(log(k))-time algorithm
    """
    retval = 1 % p
    while k:  # k != 0
        if k % 2:  # k odd
            retval = retval*a % p
        a = a*a % p
        k = k >> 1
    return retval


def legendre_ch(p):
    """
    :param p: an odd prime
    :return: the mod p Legendre character
    """
    if not isprime(p) or p == 2:
        raise ValueError("%d is not an odd prime." % (p, ))

    def ch(a):
        if a % p == 0:
            return 0
        rr = modpow2(a, (p - 1) / 2, p)
        return (-1) if (rr == p - 1) else +1

    return ch


def gcd(a, b):
    """
    :param a: an integer
    :param b: an integer
    :return: the greatest common divisor of a and b.

    Uses the Euclidean algorithm.
    """
    if a == 0 or b == 0:
        return abs(a + b)
    while a % b != 0:
        a, b = b, a % b
    return abs(b)


def sgn(a):
    """
    :param a: an integer
    :return: the sign of a (+1,-1, or 0)
    """
    if a > 0:
        return +1
    elif a < 0:
        return -1
    else:
        return 0


def euclidean_algorithm(a, b):
    """
    :param a: an integer
    :param b: an integer
    :return: x, y such that :math:`xa + yb = gcd(a, b)`.
    """

    if b == 0:
        return sgn(a), 0

    q, r = divmod(a, b)  # a = q * b + r

    if r == 0:
        return 0, sgn(b)

    # d = gcd(a, b) = gcd(b, r)
    #
    # Suppose xb+yr = d
    #
    #  xb + y(a-qb) = d
    #
    #  ya + (x-yq)b = d

    x, y = euclidean_algorithm(b, r)
    return y, (x-y*q)


def crt(r1, m1, r2, m2):
    """
    :param r1: the first residue
    :param m1: the first modulus
    :param r2: the second residue
    :param m2: the second modulus
    :return: the smallest positive simultaneous solution 'x' to the congruences
             x = r1 (mod m1)
             x = r2 (mod m2)

    Raises a ValueError if no solution exists.
    Note that we do *not* require m1 and m2 to be
    relatively prime.
    """
    c1, c2 = euclidean_algorithm(m1, m2)
    g = c1*m1+c2*m2
    q, r = divmod(r1-r2, g)
    if r != 0:  # no solution
        raise ValueError()
    else:
        x = r1 - q*(c1*m1)
        # = r2 + q*(c2*m2)
        return x % (m1*m2//g)


def ea3(a, b, c):
    """
    :param a: an integer
    :param b: an integer
    :param c: an integer
    :return: x, y, z such that :math:`xa + yb + zc = gcd(a, b, c)`.

    Uses the Euclidean algorithm twice.
    """
    x, y = euclidean_algorithm(a, b)
    # now x*a+y*b=gcd(a, b)
    s, t = euclidean_algorithm(gcd(a, b), c)
    return s*x, s*y, t


class UtilitiesTest(unittest.TestCase):

    def runTest(self):
        pass

    def test_powerset(self):
        nn = len(list(powerset([1, 2, 3, 4, 5])))
        self.assertEqual(nn, 2**5, 'number of subsets should = 32. was %d' % (nn, ))

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

    def test_euclidean_algorithm2(self):
        for a in range(-100, +100):
            for b in range(-100, +100):
                x, y = euclidean_algorithm(a, b)
                self.assertEqual(gcd(a, b), x*a + y*b, 'gcd != x*a+y*b')
