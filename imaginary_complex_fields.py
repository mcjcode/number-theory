# -*- coding: utf-8 -*-


from math import sqrt
from fractions import Fraction
import numpy as np
from utilities import (
    gcd,
    sqrtInt,
    issq,
    isprime,
    )
from itertools import islice
from quadratic_extensions import (
    factorize_in,
    minkowski_bound,
    )

ordinar_chars_lower = "`1234567890-=qwertyuiop[]\asdfghjkl;'zxcvbnm, ./"
special_chars_lower = "`¡™£¢∞§¶•ªº–≠œ∑´®†¥¨ˆøπ“‘«åß∂ƒ©˙∆˚¬…æΩ≈ç√∫˜µ≤≥÷"


def form_disc(a, b, c):
    return b**2 - 4*a*c


def genus(a, b, c):
    def f(x, y):
        return a*x**2 + b*x*y + c*y**2
    disc = form_disc(a, b, c)
    val = sorted(list(set([f(x, y) % (-disc) for x in range(-disc) for y in range(-disc)])))
    return filter(lambda x: gcd(x, -disc) == 1, val)


def proper_reduced_form(a, b, c):
    """
    Return a reduced form properly equivalent to

    ax^2 + bxy + cy^2
    """

    # first check that the discriminant is negative
    disc = form_disc(a, b, c)
    if disc >= 0:
        raise ValueError("discriminant is non-negative")
    if a < 0 or c < 0:
        raise ValueError("form is not positive_definite")

    while abs(b) > a or abs(b) > c:
        if abs(b) > a:
            r = -np.sign(b)*int(abs(b) / (2.0*a) + 0.5)
            b, c = b + 2*a*r, a*r*r + b*r + c
        elif abs(b) > c:
            r = -np.sign(b)*int(abs(b) / (2.0*c) + 0.5)
            a, b = c*r*r + b*r + a, b + 2*c*r
    if c < a:
        a, b, c = c, -b, a
    if a == -b:
        b, c = b + 2*a, c + a + b
    if a == c:
        b = abs(b)
    return a, b, c


def all_reduced_forms(disc):
    """
    Return a list of all reduced forms with
    negative discriminant D.
    """
    if disc >= 0:
        raise ValueError("discriminant must be negative")
    retval = []
    for a in range(1, int(sqrt(-disc / 3.0) + 1.0)):
        for b in range(-a+1, a+1):
            c, r = divmod(b ** 2 - disc, 4 * a)
            if r == 0:
                if gcd(gcd(a, abs(b)), c) == 1:
                    if a <= c and ((a != c) or (b >= 0)):
                        retval.append((a, b, c))
    return sorted(list(set(retval)))


def pids():
    disc = -1
    while True:
        nforms = len(all_reduced_forms(disc))
        if nforms == 1:
            if disc % 4 == 0:
                yield -disc // 4
            elif disc % 4 == 1:
                yield (1-disc) // 4
            else:
                pass
        disc -= 1


def idoneal(max_abs_disc):
    for D in range(-4, -max_abs_disc, -4):
        forms = all_reduced_forms(D)
        nforms = len(forms)
        ngenera = len(list(set([genus(*form)[0] for form in forms])))
        if nforms == ngenera:
            print(D, end=' ')


def principal_representations(n, m):
    """
    yield a sequence of positive (a, b) such that
    a**2 + n * b**2 = m
    """
    ub = sqrtInt(m/n)
    i = 0
    while i <= ub:
        x = m - n*i**2
        if issq(x):
            yield sqrtInt(x), i
        i += 1


def repmod11(p):
    princ_reps = list(principal_representations(11, p))
    if len(princ_reps) > 0:
        return u"%d\u00B2 + 11\u00B7%d\u00B2" % princ_reps[0]
    else:
        princ_reps = list(principal_representations(11, 3*p))
        (a, b) = princ_reps[0]
        if (a-b) % 3 == 0:
            c = (a-b)/3
            if c > 0:
                return u"3\u00B7%d\u00B2 + 2\u00B7%d\u00B7%d + 4\u00B7%d\u00B2" % (c, c, b, b)
            else:
                c = -c
                return u"3\u00B7%d\u00B2 - 2\u00B7%d\u00B7%d + 4\u00B7%d\u00B2" % (c, c, b, b)
        else:
            c = (a+b)/3
            if c > 0:
                return u"3\u00B7%d\u00B2 - 2\u00B7%d\u00B7%d + 4\u00B7%d\u00B2" % (c, c, b, b)
            else:
                c = -c
                return u"3\u00B7%d\u00B2 + 2\u00B7%d\u00B7%d + 4\u00B7%d\u00B2" % (c, c, b, b)


def repsmod11():
    p = 1
    gen = genus(1, 0, 11)
    while True:
        if (p % 44) in gen and isprime(p):
            yield p, repmod11(p)
        p += 1


def makelist():
    for x in list(islice(repsmod11(), 100)):
        print("%5d =" % (x[0], ), end='')
        print(x[1])


def foo():
    for D in range(-4, -8000, -4):
        forms = all_reduced_forms(D)
        # if all([form[1] == 0 for form in forms]):
        # if (len(forms) in [1, 2, 4, 8, 16, 32, 64, 128, 256]):
        if all([(a == b or a == c or b == 0) for (a, b, c) in forms]):
            # for form in all_reduced_forms(D):
            #     print D, #form#, genus(*form)
            yield -D / 4


def class_groups_of_size(k):
    for D in range(-4, -5000, -4):
        forms = all_reduced_forms(D)
        # if all( [(a==b or a==c or b==0) for (a, b, c) in forms]):
        # if all([form[1] == 0 for form in forms]):
        # if (len(forms) in [1, 2, 4, 8, 16, 32, 64, 128, 256]):
        if len(forms) == k:
            for form in forms:
                print(D, form, genus(*form))
            print('')


def ideal_class_group_info(d):
    """
    Calculate info relevant to the class group of Q[√-d].
    """
    mb = minkowski_bound(d)

    print("Minkowski bound.  All ideal classes contain an ideal of norm ≤ %d." % (mb, ))

    # split_primes = []
    for p in filter(isprime, range(2, mb+1)):
        fact = factorize_in(p, -d)
        print(fact)


def cf_from_rational(n, d):
    """
    Yield the terms of the continued fraction
    representation of n/d
    """
    while n % d:
        q, r = divmod(n, d)
        yield n//d
        n, d = d, r
    yield n


def rational_from_cf(cf):
    """
    Return the rational number (a Fraction)
    corresponding to a continued fraction.
    """
    if len(cf) == 0:
        return Fraction(1, 1)
    elif len(cf) == 1:
        return Fraction(cf[0], 1)
    else:
        return Fraction(cf[0], 1) + 1 / rational_from_cf(cf[1:])


def sum_sq_rep(p):
    """
    A prime p=2 or congruent to 1 (mod 4) has a 
    representation as a sum of 2 squares.
    Return such a representation
    """
    if p == 2:
        return 1, 1

    for m in range(int(np.sqrt(p)), p//2+1):
        cf = list(cf_from_rational(p, m))
        ln = len(cf)
        if cf == list(reversed(cf)):
            if ln % 2 == 0:
                a = rational_from_cf(cf[:ln//2]).numerator
                b = rational_from_cf(cf[:ln//2-1]).numerator
                return a, b
