
import numpy as np
import itertools
from collections import defaultdict
from utilities import powerset

def comb(k, a, b):
    """
    Yield the sequence of all k element subsets of the
    consecutive list of integers [a,..,b]
    """
    if k <= b-a+1:  # there are enough left to choose k
        if k == 0:
            yield []
        else:
            for tl in comb(k-1, a+1, b):
                yield [a] + tl
            for tl in comb(k, a+1, b):
                yield tl


def nchoosek(n, ks, p=0):
    """
    n: a non-negative integer
    ks: either an integer in the range 0..n,
        or a list of non-negative integers
        k1,...,km whose sum is n
    returns: The multinomial coefficient
        (n choose k1 . . . km)
    """
    if type(ks) == int:
        if ks<0 or ks>n:
            return 0
        ks = [ks, n-ks]

    if p==0:
        retval = 1
        i = n
        while i > 0:
            for k in ks:
                j = 1
                while j <= k:
                    retval *= i
                    i -= 1
                    retval //= j
                    j += 1
        return retval
    else:
        retval = 1
        i = n
        while i > 0:
            for k in ks:
                j = 1
                while j<=k:
                    retval = (retval*i)  % p
                    i -= 1
                    retval = (retval*pow(j, -1, p)) % p
                    j += 1
        return retval


def stirling2(n, k):
    """
    n: a non-negative integer
    k: a non-negative integer
    returns: the number of equivalence relations of
        the numbers [1,...,n] with k classes.

    We use count surjections from [1,..,n] to [1,..,k] by inclusion
    exclusion, then divide by k!
    """
    retval = 0
    nchk = 1
    for j in range(k+1):
        retval += (-1)**(j%2)*nchk*(k-j)**n
        nchk = nchk * (k-j) // (j+1)
    for j in range(1, k+1):
        retval //= j
    return retval


def eqrels(N):
    """
    Return all possible equivalence relations on [0,n)
    """
    h = [[]] # the equivalence relations on [0,1)
    for n in range(N):
        h2 = []
        for r in h:
            for j1 in range(len(r)):
                h2.append([s+([n]*(j2==j1)) for j2, s in enumerate(r)])
            h2.append(r + [[n]])
        h = h2
    return h


def monomials(d: int, n: int):
    """
    Yield the exponent vectors of all degree
    d monomials in n variables
    """
    for v in comb(n-1, 1, d+n-1):
        yield [v[0]-1] + [v[i+1]-v[i]-1 for i in range(n-2)] + [d+n-1-v[-1]]


def conjugate_partition(partition):
    """
    Return the conjugate partition
    """
    return [len([1 for part in partition if part >= k])
            for k in range(1, max(partition)+1)]


def partitions(d, max_part=None):
    """
    Return a generator yielding the partitions of d
    """
    if max_part is None:
        max_part = d

    if d == 0:
        yield []
    else:
        for p in range(max_part, 0, -1):
            for i in range(1, d//p+1):  # number of parts of size p
                seg = [p]*i
                for part in partitions(d-i*p,p-1):
                    yield seg + part


def gray(n):
    """
    The nth term in the reflective
    gray code.
    """
    n -= 1
    return n^(n>>1)


def graybits():
    """
    a generator yielding the sequence
    of bit locations needed to flip
    to product the reflective gray code.
    """
    i = 1
    while True:
        a = gray(i)^gray(i+1)
        retval = 0
        while a:
            a //= 2
            retval += 1
        yield retval-1
        i += 1


def newtons_identity(ps):
    """
    :param ps: the values of the power symmetric polynomials
    :return: the corresponding values of the elementary symmetric polynomials

    Example:

    x1, x2, x3 = 1, 2, 3
    ps = [ x1**n + x2**n + x3**n for n in range(3)]
    es = newtons_identities(ps)
    assert es == [1, 6, 11, 6]
    
    """
    es = [1]
    for j in range(1, len(ps)):
        z = 0
        sign = 1
        for i in range(1, j+1):
            z += sign*ps[i]*es[j-i]
            sign *= -1
        es.append(z//j)
    return es
