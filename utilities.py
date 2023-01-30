#!/usr/bin/env python -i
# -*- coding: utf-8 -*-
"""
General purpose, factorization and modular
arithmetic routines.
"""

import time
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


def is_miller_rabin_witness(p, a):
    """
    :param p: an odd number >= 3 whose primality we are testing
    :param a: the potential 'witness' a
    """
    s = 0
    d = p-1
    while d%2==0:
        s += 1
        d //= 2
    # now p-1 = 2**s * d with d odd.
    assert p == 2**s * d + 1
    assert d%2 == 1

    ad = modpow2(a, d, p)    ##p, d, a)
    if ad==1 or ad==p-1:
        return False
    for _ in range(s-1):
        ad *= ad
        ad %= p
        if ad==p-1:
            return False
    return True


def isprime_miller_rabin(p):
    if p<=1:
        return False
    if p==2:
        return True
    if p%2==0:
        return False
    # p is now an odd number >= 3.
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= p:
            break
        if is_miller_rabin_witness(p, a):
            return False
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
    b = a % p
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


def bezout(a:int, b:int) -> (int, int):
    """
    :param a: an integer
    :param b: an integer
    :return: x, y such that :math:`xa + yb = gcd(a, b)`.
    """

    sa = sgn(a)
    sb = sgn(b)
    
    a = abs(a)
    b = abs(b)
    
    if b == 0:
        return sa, 0
    
    x0, y0 = 1, 0 # x0*a + y0*b = a
    x1, y1 = 0, 1 # x1*a + y1*b = b
    q, r = divmod(a, b) # now r = a - q*b

    while r:
        a, b = b, r
        (x0, y0), (x1, y1) = (x1, y1), (x0 - q*x1, y0 - q*y1)
        q, r = divmod(a, b)
        
    return x1*sa, y1*sb
        
def modinv(m:int, a:int) -> int:
    """
    :param m: a positive integer
    :param a: an integer, with (m,a)==1
    :return: the multiplicative inverse of a(mod m)
    """
    x, y = bezout(m, a)

    if x*m + y*a != 1:
        raise ValueError(f'{a} is not relatively prime to {m}')
    #
    # now we have xp+ya=1
    #
    return y % m


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
    c1, c2 = bezout(m1, m2)
    g = c1*m1+c2*m2
    q, r = divmod(r1-r2, g)
    if r != 0:  # no solution
        raise ValueError()
    else:
        x = r1 - q*(c1*m1)
        return x % (m1*m2//g)


def ea3(a, b, c):
    """
    :param a: an integer
    :param b: an integer
    :param c: an integer
    :return: x, y, z such that :math:`xa + yb + zc = gcd(a, b, c)`.

    Uses the Euclidean algorithm twice.
    """
    x, y = bezout(a, b)
    assert x*a + y*b == gcd(a, b)
    s, t = bezout(gcd(a, b), c)
    return s*x, s*y, t


def digitsum(n,base=10):
    """
    :param n: a non-negative integer
    :return: the sum of all the digits of n in the given base
    """
    retval = 0
    while n:
        k = n % base
        retval += k
        n = (n-k)//base
    return retval


def digits(n, base=10):
    """
    :param n: an integer
    :return: a list ds so that the base^i digit of n is ds[i]
    """
    while n:
        n, d = divmod(n, base)
        yield d


def num_from_digits(ds, base=10):
    """
    :param ds: a list of integers from [0,base)
    :return: the integer whose base^i digit is ds[i]
    """
    n = 0
    a = 1
    for d in ds:
        n += a*d
        a *= base 
    return n


def sign(x):
    """
    Return the sign of a number
    """
    if x < 0:
        return -1
    if x > 0:
        return +1
    return 0
