#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

from __future__ import print_function
import numpy as np
from utilities import primitive_root, modpow, squarefree, isprime
from gaussian_integer import GaussianInteger


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
    Returns the Jacobi sum J(chi, chi) associated with the
    quartic character which maps the smallest primitive
    root mod p to the pure imaginary unit i.
    """
    a = primitive_root(p)
    zeta = GaussianInteger(0, 1)
    avals = [0] + [modpow(a, k, p) for k in range(1, p)]
    achars = [GaussianInteger(0)] + [zeta ** (k % 4) for k in range(1, p)]
    bvals = [(1 - a) % p for a in avals]
    binds = [avals.index(b) for b in bvals]
    bchars = [achars[bi] for bi in binds]
    retval = GaussianInteger(0)
    for val in (ac * bc for (ac, bc) in zip(achars, bchars)):
        retval = retval + val
    return retval


def run1():
    """
    Calculate the Jacobi sums of a cubic character
    for some finite fields of order p^2
    """

    zeta3 = np.cos(2*np.pi/3) + 1j*np.sin(2*np.pi/3)

    # These are all going to
    for p in filter(isprime, range(7, 200, 12)):

        nnp = p if (p % 4 == 1) else p*p

        # print 'Finite Field with %d^2=%d elements' % (p, p*p)
        def mul(a, b):
            # Here's how we multiply Gaussian integers
            return (a[0]*b[0]-a[1]*b[1]) % p, (a[0]*b[1]+a[1]*b[0]) % p

        def add(a, b):
            # Here's how we add Gaussian integers
            return ((a[0]+b[0]) % p), ((a[1]+b[1]) % p)

        def order(a):
            retval = 1
            b = a
            while b != (1, 0):
                b = mul(b, a)
                retval += 1
            return retval

        def prim_elt():
            for x in range(p):
                for y in range(p):
                    elt = (x, y)
                    if elt == (0, 0):
                        continue
                    if order(elt) == p*p-1:
                        return elt

        x = prim_elt()
        powers = [(0, 0), (1, 0)]
        cnt = 0
        while cnt < nnp-2:
            powers.append(mul(powers[-1], x))
            cnt += 1

        chars = [0.0] + list(map(lambda n: zeta3**n, range(p*p)))
        avals = powers
        bvals = map(lambda y: add((1, 0), mul((-1, 0), y)), avals)
        bindices = [avals.index(b) for b in bvals]
        bchars = [chars[bi] for bi in bindices]
        jsum = sum(a*b for (a, b) in zip(chars, bchars))

        # express jsum as an integral combination of 1 and zeta3=(-1+sqrt(3))/2
        #
        c2 = jsum.imag/(np.sqrt(3.0)/2.0)
        tmp = jsum - c2*zeta3
        c1 = tmp.real

        print('%6d %10d %35s %7.1f %7.1f %35s %s' % (p, p*p, jsum, c1, c2, c1+c2*zeta3, '%d+%di' % x))


def f(x, y, p):
    return (y*y - (x-1)*x*(x+1)) % p


def run2():
    for p in []:
        if not isprime(p):
            continue
            # [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61] :
        cnt = 0
        for x in range(p):
            for y in range(p):
                cnt += f(x, y, p) == 0
        bound = 2*int(np.sqrt(p))
        if abs(cnt-p) >= bound:
            print('%5d %5d %5d %5d' % (p, cnt, cnt - p, bound))


def jacobi2(a, n):
    if n==1:
        return 1
    assert n % 2 == 1
    a %= n
    if a == 0:
        return 0
    t = 1
    while a:  # != 0:
        while a % 2 == 0:
            a /= 2
            r = n % 8
            if r == 3 or r == 5:
                t = -t
        a, n = n, a
        if a % 4 == n % 4 == 3:
            t = -t
        a %= n
    return t if n == 1 else 0


def partial_jacobi_sum(modulus, ub):
    """
    add up (j/modulus) for j in [0..ub].  use a sieve
    approach (though it doesn't seem to be buying us
    much over a good jacobi symbol implementation.
    """

    # modulus must be odd and squarefree
    assert squarefree(modulus) and (modulus % 2)

    m = modulus

    J = np.ones((ub+1, ), dtype=np.int8)  # this will only store 0, +1, -1's
    P = np.ones((ub+1, ), dtype=np.int8)
    P[0] = 0
    P[1] = 0
    J[0] = 0
    for p in range(2, ub+1):
        if P[p]:  # p is prime
            P[p*p::p] = 0  # strike higher multiples
            if m % p:  # p does not divide m
                poverm = jacobi2(p, m)
                if poverm == 1:
                    pass
                else:
                    pk = p
                    while pk <= ub:
                        J[pk::pk] *= -1
                        pk *= p
            else:
                J[p::p] = 0
    return J
