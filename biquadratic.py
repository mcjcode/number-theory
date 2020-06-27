#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

from __future__ import print_function

from gaussian_integer import GaussianInteger
from utilities import isprime, modpow
from jacobi import jacobi_sum_quartic

def maptochar(x, pi):
    zero = GaussianInteger(0)
    roots = map(lambda nn: GaussianInteger(0, 1)**nn, range(4))
    for root in roots:
        if (x-root) % pi == zero:
            return root


def biquad(pi1, pi2):
    if type(pi2) == int:
        pi2 = GaussianInteger(pi2)

    if pi2 == GaussianInteger(0):
        return GaussianInteger(0)

    z = pi2 % pi1
    ch = modpow(z, (pi1.norm()-1)/4, pi1)
    return maptochar(ch, pi1)


class Character(object):

    def __init__(self, xx):
        if type(xx) == GaussianInteger:
            self.f = lambda aa: biquad(xx, aa)
        else:
            self.f = xx

    def __call__(self, aa):
        return self.f(aa)

    def __mul__(self, other):
        def f(aa):
            return self.f(aa)*other.f(aa)
        return Character(f)


def jacobi_sum(chi1, chi2, p):
    """
    Return the Jacobi sum of two mod p characters.
    """
    retval = GaussianInteger(0)
    for aa in range(p):
        aa = GaussianInteger(aa)
        retval = retval + chi1(a)*chi2(GaussianInteger(1)-aa)
    return retval


def flip(gint):
    return GaussianInteger(gint.real(), abs(gint.imag()))


def table2():
    maxp = 75
    q2primes = [flip(jacobi_sum_quartic(qq)) for qq in range(5, maxp, 4) if isprime(qq)]

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
    qprimes = [(flip(jacobi_sum_quartic(qq)) if (qq % 4 == 1) else GaussianInteger(qq))
               for qq in range(3, maxp, 2) if isprime(qq)]
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
