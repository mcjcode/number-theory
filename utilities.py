#!/usr/bin/env python -i
# -*- coding: utf-8 -*-
"""
General purpose, factorization and modular
arithmetic routines.
"""

import time
import itertools
import functools
import operator
import random

import math
from math import sqrt


def prod(xs, start=1):
    """
    :param xs: a sequence of elements to multiply
    :param start: the value to return if xs is empty
    :return: the product of the xs
    """
    return functools.reduce(operator.mul, xs, start)


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


def sqrtInt(a:int) -> int:
    """
    :param a: a non-negative integer
    
    :return: the largest integer smaller than or equal to the square
    root of n
    """
    if a==0:
        return 0

    n = 1 << a.bit_length()//2
    q, r = divmod(a, n)
    n0 = (n + (q+(r>0)))//2
    while n0!=n:
        n = n0
        q, r = divmod(a, n)
        n0 = (n + (q+(r>0)))//2

    if n0*n0>a:
        n0 -= 1
        
    return n0

def cbrtInt(n):
    """
    :param n: a non-negative integer
    :return: the largest integer smaller than the cube root of n

    Note that this depends on `math.sqrt`, and its double precision
    accuracy means that this function should not be trusted for n on
    the order of :math:`10^{52}` and up.
    """
    cbrtn = int(n**(1./3.))
    if (cbrtn+1)**3 <= n:
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

    assert p-1 == 2**s * d
    assert d%2 == 1

    ad = pow(a, d, p)
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


def lucas(n, nattempts=5, factoring_steps=0):
    if n==1:
        return False, None, '1 is a unit'
    if n==2:
        return True, 1, 'Lucas: Found a primitive root'
    assert n>2

    try:
        f = list(factorize2(n-1, factoring_steps))
    except Exception as e:
        return 'Maybe', None, f'Lucas: Too hard to factor n-1'

    for _ in range(nattempts):
        a = random.randint(1, n)
        g = gcd(a, n)
        if g>1 and g!=n:
            return False, g, 'Lucas: Found a factor'
        
        if pow(a, n-1, n) != 1:
            return False, a, f'Lucas: order({a})!=n-1. n not prime'
        
        if all(pow(a, (n-1)//p, n)!=1 for p, _ in f):
            return True, a, 'Lucas: found a primitive root'

    return 'Maybe', None, 'Lucas: test inconclusive'

# if n is a square, then _sqlowbits[n%63] is True
_sqlowbits = [False]*64
for x in range(64):
    _sqlowbits[(x*x)%63] = True

def issq(n: int) -> bool:
    # first check the low 5 bits (52/64=87.5% of
    # numbers will fail this check)
    return n>=0 and _sqlowbits[n%63] and sqrtInt(n)**2==n


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

_wheel = [ 1, 6, 5, 4, 3, 2,
           1, 4, 3, 2, 1, 2,
           1, 4, 3, 2, 1, 2,
           1, 4, 3, 2, 1, 6,
           5, 4, 3, 2, 1, 2 ]


def trial_division(n, ntrials=0, ps=[2, 3, 5]):
    """
    :param n: a positive integer
    :return: yields a sequence of (prime, exponent) pairs where the
             primes are distinct and in increasing order giving
             the prime factorization of n.
    """

    itrial = 0
    
    for p in ps:
        if p*p > n:
            if n>1:
                yield n, 1
            return
        e = 0
        while n%p == 0:
            e += 1
            n //= p
        if e:
            yield p, e
        itrial += 1
        if ntrials and itrial==ntrials:
            raise Exception(f'Exceeded {ntrials} trial divisions.  Giving up.')

    q = ps[-1]+2
    while q*q <= n:
        e = 0
        while n%q == 0:
            e += 1
            n //= q
        if e:
            yield q, e
        q += _wheel[q%30]
        itrial += 1
        if ntrials and itrial==ntrials:
            raise Exception(f'Exceeded {ntrials} trial divisions.  Giving up.')

    if n > 1:
        yield n, 1

factorize2 = trial_division

def factors(f):
    """
    Given a factorization `f` of an integer `n`, return
    a list of all of the factors of `n`
    """
    retval = [1]
    for p, e in f:
        xs = [b*p**ee for b in retval for ee in range(1, e+1)]
        retval += xs
    return retval


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
        for q, e in facts:
            if pow(a, (p-1)//q, p)==1:
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


def fastpow(a, k, mul, one):
    """
    :param a: the base
    :param k: the exponent (a positive integer)
    :param mul: the binary multiplication operator
    :param one: the multiplicative identity for mul
    :return: :math:`a^k`.
    """
    retval = one
    while k:  # k != 0
        if k % 2:  # k odd
            retval = mul(retval, a)
        a = mul(a, a)
        k = k >> 1
    return retval


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
        rr = pow(a, (p-1)//2, p)
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


def bezout(a:int, b:int) -> (int, int):
    """
    :param a: an integer
    :param b: an integer
    :return: x, y such that :math:`xa + yb = gcd(a, b)`.
    """

    sa = sign(a)
    sb = sign(b)
    
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
    :param a: an integer
    :return: the sign of a (+1,-1, or 0)
    """
    if x < 0:
        return -1
    if x > 0:
        return +1
    return 0


def step_modp_pascal(row, p):
    """
    :param row: a row of Pascal's triangle mod p, in compressed form
    :param p: the modulus, not necessary prime
    :return: the next row of Pascal's triangle mod p, in compressed form
    
    compressed mod p pascal's triangle generation.
    
    Given the nth row of pascal's triangle, generate the
    (n+1)st row, where both are represented in compressed
    form enumerating just the non-zero entries.
    
    e.g. row=[(0,1),(1,1),(4,1),(5,1)] is the 5th row
    of the mod 2 triangle, then the result would be
    [(0,1),(2,1),(4,1),(6,1)], the 6th row.
    """

    new_row = [(0,1)]
    for k in range(len(row)-1):
        if row[k][0]+1 == row[k+1][0]:
            newval = (row[k][1] + row[k+1][1]) % p
            if newval:
                new_row.append((row[k+1][0], newval)) # if p=2, you never land here
        else:
            new_row.append((row[k][0]+1,row[k][1]))
            new_row.append((row[k+1][0],row[k+1][1]))
    new_row.append((row[-1][0]+1,1))
    return new_row


