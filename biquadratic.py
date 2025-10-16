# -*- coding: utf-8 -*-

from __future__ import print_function

from gaussian_integer import GaussianInteger
from utilities import isprime, modpow2
from jacobi import jacobi_sum_quartic


def maptochar(x, pi):
    zero = GaussianInteger(0)
    roots = map(lambda n: GaussianInteger(0, 1)**n, range(4))
    for root in roots:
        if (x-root) % pi == zero:
            return root


def biquad(pi1, pi2):
    if type(pi2) == int:
        pi2 = GaussianInteger(pi2)

    if pi2 == GaussianInteger(0):
        return GaussianInteger(0)

    z = pi2 % pi1
    ch = modpow2(z, (pi1.norm()-1)/4, pi1)
    return maptochar(ch, pi1)


class Character(object):

    def __init__(self, x):
        if type(x) == GaussianInteger:
            self.f = lambda a: biquad(x, a)
        else:
            self.f = x

    def __call__(self, a):
        return self.f(a)

    def __mul__(self, other):
        def f(a):
            return self.f(a)*other.f(a)
        return Character(f)


def jacobi_sum(chi1, chi2, p):
    """
    Return the Jacobi sum of two mod p characters.
    """
    retval = GaussianInteger(0)
    for a in range(p):
        a = GaussianInteger(a)
        retval = retval + chi1(a)*chi2(GaussianInteger(1)-a)
    return retval


def flip(gint):
    return GaussianInteger(gint.real(), abs(gint.imag()))


def table2():
    maxp = 75
    q2primes = [flip(jacobi_sum_quartic(q)) for q in range(5, maxp, 4) if isprime(q)]

    for pi1 in q2primes:
        for pi in [pi1, GaussianInteger(pi1.real(), -pi1.imag()),
                   GaussianInteger(-pi1.imag(), pi1.real()),
                   GaussianInteger(-pi1.real(), -pi1.imag()),
                   GaussianInteger(pi1.imag(), -pi1.real()),
                   GaussianInteger(pi1.real(), -pi1.imag()),
                   GaussianInteger(pi1.imag(), pi1.real()),
                   GaussianInteger(-pi1.real(), pi1.imag()),
                   GaussianInteger(-pi1.imag(), -pi1.real())]:

            chi = Character(pi)
            chi2 = chi*chi

            print('%10s %10s %10s %10s' % (pi, chi(GaussianInteger(-1)), jacobi_sum(chi, chi, pi.norm()),
                                           jacobi_sum(chi, chi2, pi.norm())))


def biquadratic_residue_table():
    print("""A table of biquadratic residues.  Each column shows
    the values of the character associated with the prime at the
    top of the column
    """)

    maxp = 75
    qprimes = [(flip(jacobi_sum_quartic(q)) if (q % 4 == 1) else GaussianInteger(q))
               for q in range(3, maxp, 2) if isprime(q)]
    pprimes = qprimes

    print('%5s' % ' ', end=' ')
    for pi2 in qprimes:
        print('%4s' % (GaussianInteger(pi2.real(), 0), ), end=' ')
    print('')

    print('%5s' % ' ', end=' ')
    for pi2 in qprimes:
        if pi2.imag() == 0:
            print('%4s' % ' ', end=' ')
        else:
            print('%+4s' % (GaussianInteger(0, pi2.imag()), ), end=' ')
    print('')

    for pi1 in pprimes:
        print('%5s' % (pi1, ), end=' ')
        for pi2 in qprimes:
            if pi1 == pi2:
                print('%4s' % ' ', end=' ')
            else:
                print('%4s' % (biquad(pi2, pi1), ), end=' ')
        print('')
