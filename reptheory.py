#!/usr/bin/env python

import sympy
from combinatorics import partitions


def hook_length(partition, i, j):
    return (partition[i]-j) + (len([1 for part in partition if part > j])-i) - 1


def print_hook_length_diagram(p):
    s = ''
    for i,pi in enumerate(p):
        for j in range(pi):
            s += '%3d ' % (hook_length(p,i,j))
        s += '\n'
    print(s)


def hook_length_formula(p):
    n = sum(p)
    return sympy.factorial(n) // sympy.prod([hook_length(p, i, j) for i, pi in enumerate(p) for j in range(pi)])


def nvariables(lbl, n):
    s = ' '.join([('%s%d' % (lbl, i+1)) for i in range(n)])
    retval = sympy.abc.symbols(s)
    if type(retval) != tuple:
        retval = (retval,)
    return retval


def discriminant(xs):
    n = len(xs)
    retval = sympy.sympify('1')
    for i in range(n-1):
        for j in range(i+1,n):
            retval = retval * (xs[i] - xs[j])
    return retval


def power_function(xs, k):
    """
    Return the kth symmetric power_function
    in the variables xs
    """
    retval = sympy.sympify('0')
    for x in xs:
        retval = retval + (x**k)
    return retval


def multiplicities(mu):
    n = max(mu)
    retval = [0]*(n+1)
    for mui in mu:
        retval[mui] += 1
    return retval


def frobenius_formula(lmda, mu):
    """
    Calculate the value of the character for the
    representation associated with the partition
    lmda on the conjugacy class given by mu.
    """
    k = len(lmda)
    ells = [lmda[j-1] + k - j for j in range(1, k+1)]
    xs = nvariables('x', k)
    disc = discriminant(xs)
    retval = disc
    mults = multiplicities(mu)
    for j in range(len(mults)):
        retval = retval * (power_function(xs, j)**mults[j])
    term = sympy.sympify('1')
    for x, ell in zip(xs,ells):
        term = term * (x**ell)
    return retval.expand().coeff(term)


def frobenius_polynomial(mu):
    """
    """
    dim = sum(mu)
    xs = nvariables('x', dim)
    disc = discriminant(xs)
    retval = disc
    for mm in mu:
        retval = retval * power_function(xs,mm)
    return xs, retval.expand()


def char_table(n):
    mus = list(reversed(list(partitions(n))))
    lmbdas = mus
    wdth = max(map(lambda x: len(str(x)), mus))
    fmt_str = '%%%ds' % wdth
    outstr = ''
    for mu in mus:
        outstr += fmt_str % mu
        outstr += ' '
        for lmbda in lmbdas:
            v = frobenius_formula(lmbda, mu)
            outstr += '%+5d' % (v,)
        outstr += '\n'
    return outstr


def char_table2(d):
    mus = list(reversed(list(partitions(d))))
    lmbdas = [list(mu) for mu in mus]
    wdth = max(map(lambda x: len(str(x)), mus))
    fmt_str = '%%%ds' % (wdth,)
    outstr = ''
    for mu in mus:
        xs, frob_poly = frobenius_polynomial(mu)
        outstr += fmt_str % (mu,)
        outstr += ' '
        for lmbda in lmbdas:
            lmbda += [0]*(d-len(lmbda))
            ells = [lmbda[j-1] + d - j for j in range(1, d+1)]
            term = sympy.sympify('1')
            for x, ell in zip(xs, ells):
                term = term * (x**ell)
            v = frob_poly.coeff(term)
            outstr += '%+5d' % (v,)
        outstr += '\n'
    return outstr
