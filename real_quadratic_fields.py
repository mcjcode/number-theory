#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

from math import sqrt, floor
import mpmath

from utilities import (
    isprime,
    gcd,
    sqrtInt,
    issq,
    prod,
    )

from quadratic_extensions import (
    QuadInt,
    factorize_in,
    discriminant,
    minkowski_bound,
    )


def cont_frac(m):
    r"""
    :param m: a positive floating point number 'm',
    :return: a generator for the sequence of integers :math:`a_0,a_1,a_2,\dots`
             that give the continued fraction representation of m:

    :math:`m = a_0 + \cfrac{1}{a_1 + \cfrac{1}{a_2 + \cfrac{1}{a_3 + \cdots}}}`
    You probably shouldn't use this if you really care about all of the terms
    since rounding eventually corrupts the process.  For numbers in quadratic
    number fields use ‘cont_frac_quad’ instead.
    """
    while True:
        f = int(floor(m))
        yield f
        if m == f:
            return
        m = 1./(m-f)


_prec = 56
_multint = 2**_prec
_multflt = 2.0**_prec


def divint(num, denom, prec=15):
    r"""
    Return a floating point number representing num/denom
    accurate to prec decimal digits of accuracy.

    Use this when the integers num and denom are not
    representable as floating point numbers (too big) but
    the quotient is.

    Don’t bother using prec>15 as the arithmetic is
    single precision floats.
    """

    return (num*_multint)/denom/_multflt


def floor_rem(x):
    fl = floor(x)
    return int(fl), x - fl


def cont_frac_quad(a, b, c, d):
    r"""
    Yield the continued fraction for (a+b*sqrt(d))/c
    """
    sqd = mpmath.sqrt(d)
    while True:
        m = int((mpmath.mpf(a)+sqd*mpmath.mpf(b))/mpmath.mpf(c))
        yield m
        a1p = (a-c*m)
        a1 = c*a1p
        b1 = -c*b
        c1 = a1p**2 - d*b**2
        g = gcd(c*gcd(a1p, b), c1)
        a, b, c = a1//g, b1//g, c1//g


def cont_frac_quad2(P, Q, n):
    """
    Yield the repeating simple continued fraction
    for the number (P+sqrt(n))/Q.  The integer Q
    must divide n-P*P.
    """
    A0, B0 = 1, 0
    q0 = sqrtInt(n)
    q, A, B = q0, q0, 1
    yield P, Q, q, A, B
    while True:
        assert (n-P*P) % Q == 0
        P = Q*q-P
        Q = (n-P*P)//Q
        q = (P+q0)//Q

        A, A0 = q*A + A0, A
        B, B0 = q*B + B0, B
        yield P, Q, q, A, B
        if q==2*q0:
            break


def approximants(d):
    r"""
    :param d: a positive non-square integer
    :return: Yields the approximants for :math:`\sqrt{d}`
    """
    h0, k0 = 0, 1
    h1, k1 = 1, 0

    for ai in cont_frac_quad(0, 1, 1, d):
        h = ai*h1 + h0
        k = ai*k1 + k0
        c = h**2 - d*k**2
        yield QuadInt(d, h-k, 2*k) if (d % 4) == 1 else QuadInt(d, h, k)

        if c == 1:
            break
        h0, k0 = h1, k1
        h1, k1 = h,  k


def approximants2(d):
    r"""
    :param d: a positive non-square integer
    :return: Yields the approximants for :math:`\sqrt{d}`

    Here you just get (h, k) pairs, not the QuadInts
    that approximants returns.
    """
    h0, k0 = 0, 1
    h1, k1 = 1, 0
    for ai in cont_frac_quad(0, 1, 1, d):
        h = ai*h1 + h0
        k = ai*k1 + k0
        yield h, k
        if h**2 - d*k**2 == 1:
            break
        h0, k0 = h1, k1
        h1, k1 = h,  k


def pell(d):
    r"""
    :param d: a positive non-square integer
    :return: positive integers h, k such that :math:`h^2 - dk^2 = 1`.

    Uses Lagrange's method of continued fractions.
    """
    h0, k0 = 0, 1
    h1, k1 = 1, 0

    h0s, k0s = 0, 1
    h1s, k1s = 1, 0

    modulus = prod([2, 3, 5, 7, 11, 13, 17, 19, 23, 29])

    for ai in cont_frac_quad(0, 1, 1, d):
        h = ai*h1 + h0
        k = ai*k1 + k0

        hs = (ai*h1s + h0s) % modulus
        ks = (ai*k1s + k0s) % modulus

        if (hs**2 - d*ks**2) % modulus == 1:
            if h**2 - d*k**2 == +1:
                return h, k
        h0, k0 = h1, k1
        h1, k1 = h,  k
        h0s, k0s = h1s, k1s
        h1s, k1s = hs,  ks


def squares_mod_d(d):
    return sorted({(n**2) % d for n in range(d)})


def norm_search(p, d):
    md = None
    if (d % 4) != 1:
        for n in range(1, 20000):
            for m in [n**2*d+p, n**2*d-p]:
                if issq(m):
                    m2 = int(sqrt(m))
                    if (gcd(n, p) == 1) or (gcd(m2, p) == 1):
                        md = QuadInt(d, m2, n)
                        break
    return md


def fundamental_unit_old(d):
    if d <= 1:
        raise ValueError('%d is not >= 2' % (d, ))
    b = 1
    while True:
        if issq(b*b*d-4):
            a = sqrtInt(b*b*d-4)
            break
        if issq(b*b*d+4):
            a = sqrtInt(b*b*d+4)
            break
        b += 1
    if d % 4 == 1:
        return QuadInt(d, (a-b)/2, b)
    else:
        return QuadInt(d, a/2, b/2)


def fundamental_unit(d):
    r"""
    :param d: a non-square positive integer
    :return: return a fundamental unit (a QuadInt) for the ring of integers
             in :math:`\mathbb{Q}(\sqrt{d})`.
    """
    if d <= 1:
        raise ValueError('%d is not >= 2' % (d, ))

    if d % 4 == 3 or d % 4 == 2:
        return next(a for a in approximants(d) if abs(a.norm()) == 1)
    else:
        if d == 5:
            return QuadInt(5, 0, 1)
        a = next(a for a in approximants(d) if abs(a.norm()) in [1, 4])
        if abs(a.norm()) == 1:
            return a
        else:
            if a.a % 2 == 0 and a.b % 2 == 0:
                a.a /= 2
                a.b /= 2
            return a


def class_group_info(d):
    r"""
    :param d: a non-square integer
    :return: prints information relevant to the class group of the ring
             of integers in :math:`\mathbb{Q}(\sqrt{d})`
    """
    mb = minkowski_bound(d)
    disc = discriminant(d)
    print('Discriminant = %d' % disc)
    print('Minkowski Bound = %d' % (mb, ))
    split_primes = []
    for p in [p for p in range(2, mb+1) if isprime(p)]:
        fact = factorize_in(p, d)
        if fact == (p, ):
            print(fact)
        else:
            md = norm_search(p, d)
            print(fact, md)
            if md is None:
                split_primes.append(p)

    if (d % 4) != 1 and d > 0:
        print('approximant based elements and norms')
        print([(x, x.norm()) for x in approximants(d)])

    if not split_primes:
        print('No non principal primes under the Minkowski bound!')
        return
    print('Split primes')
    print(split_primes)
    print('Split primes that are not squares mod %d' % (d, ))
    sq = squares_mod_d(d)
    print([p for p in split_primes if (p not in sq) and ((d-p) not in sq)])

    for i in range(len(split_primes)):
        for j in range(i, len(split_primes)):
            p = split_primes[i]*split_primes[j]
            md = norm_search(p, d)
            print(p, md)

    for i in range(len(split_primes)):
        for j in range(i, len(split_primes)):
            for k in range(j, len(split_primes)):
                p = split_primes[i] * split_primes[j] * split_primes[k]
                md = norm_search(p, d)
                print(p, md)
