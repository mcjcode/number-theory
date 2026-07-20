from itertools import pairwise
import numpy as np
import sympy
from sympy import (
    resultant, symbols, Matrix, degree, sqf_list, expand
)

from utilities import listmap, listtake, fastpow, prod


def _matrix_modpow(a, k, p):
    one = np.identity(a.shape[0], dtype=object)
    return fastpow(a, k, mul=lambda x, y: x.dot(y) % p, one=one)


def fibonacci(a=0, b=1, modulus=0):
    """
    Yield the sequence of fibonacci numbers, reduced
    by the given modulus.
    """
    if modulus:
        while True:
            yield a
            a, b = b, (a+b) % modulus
    else:
        while True:
            yield a
            a, b = b, a+b


def tribonacci(a=0, b=0, c=1, modulus=2):
    """
    Yield the sequence of tribonacci numbers, reduced
    by the given modulus.
    """
    while True:
        yield a
        a, b, c = b, c, (a+b+c) % modulus


def lrs(cs, xs):
    """
    Yield the terms of the linear recursive sequence
    with initial terms xs and coefficients cs
    """
    xs = xs.copy()
    rcs = list(reversed(cs))
    while True:
        yield xs[0]
        *(xs[:-1]), xs[-1] = *(xs[1:]), sum(c*x for (c, x) in zip(rcs, xs))


def fastfib(n, modulus, coeffs=[1, 1], init=[0, 1]):
    """
    By default, return the n'th fibonacci number,
    where the indexing is the following:

    is F0 = 0, F1 = 1, F2 = 1, F3 = 2 ...

    Use this when you just want a single fibonacci number
    with large n.

    >>> [fastfib(i,1000) for i in range(13)]
    [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]

    More generally, this will produce the n^th term in
    a linear recursive sequence specified by the given
    coefficients and starting values.

    coeffs = [a1, a2, ...]
    init   = [F0, F1, ...]

    Fn = a1*F(n-1) + a2*F(n-2) + ...
    """
    coeffs = coeffs[::-1]
    order = len(coeffs)
    Id = np.identity(order-1, dtype=object)
    Z = np.zeros((order-1, 1), dtype=object)
    C = np.array([coeffs], dtype=object)
    F = np.block([[Z, Id], [C]])

    M = _matrix_modpow(F, n, modulus)

    return M.dot(np.array(init).transpose())[0] % modulus


def fastfib_decimate(n, modulus, coeffs=[1, 1], init=[0, 1], m=1, i=0):
    """
    Return the n'th element of the sequence

    F(i), F(m+i), F(2*m+i), ...

    where F is specified by the given coefficients and initial values.
    """
    coeffs = coeffs[::-1]
    order = len(coeffs)
    Id = np.identity(order-1, dtype=object)
    Z = np.zeros((order-1, 1), dtype=object)
    C = np.array([coeffs], dtype=object)
    F = np.block([[Z, Id], [C]])

    init = [fastfib(i+j, modulus, coeffs[::-1], init) for j in range(order)]
    F = _matrix_modpow(F, m, modulus)
    M = _matrix_modpow(F, n, modulus)

    return M.dot(np.array(init).transpose())[0] % modulus


def fastfib_sum1(n, m):
    """return the sum of the first n fibonacci numbers, mod m"""
    return fastfib(n+1, m) - 1


def fastfib_sum2(n, m):
    """return the sum of the squares of the first n fibonacci numbers, mod m"""
    return (fastfib(n-1, m) * fastfib(n, m)) % m


def fastfib_sum3(n, m):
    """return the sum of the cubes of the first n fibonacci numbers, mod m"""
    if n<=1:
        return 0
    else:
        return (
            (fastfib(3*n-1, m)+(-1)**n*6*fastfib(n-2, m)+5)*pow(10, -1, m)
        ) % m


def find_recursion(xs, deg):
    """
    Find a linear recursion of degree deg that is satisfied
    by the sequence xs.
    """
    amat = [xs[i:deg+i] for i in range(deg)]
    A = Matrix(amat)
    bmat = [xs[i+1:deg+i+1] for i in range(deg)]
    B = Matrix(bmat)
    return list((B*A.inv())[-1, :])[::-1]


def find_recursion2(xs, deg):
    """
    Find a linear recursion of degree deg that is satisfied
    by the sequence xs.
    """
    a = np.array([xs[i:deg+i] for i in range(deg)], dtype=float)
    b = np.array(xs[deg:deg+deg], dtype=float)
    vs = list(map(round, np.linalg.solve(a, b)))
    for i in range(4*deg-1):
        x = xs[i+deg]
        y = sum(v*x for (v, x) in zip(vs, xs[i:i+deg]))
        if x != y:
            raise ValueError(f'{deg=} {x=} {y=} Failed, higher degree?')
    return list(reversed(vs))


def cumulative_lrs(cs, xs):
    """
    if we have a linear recursive sequence with initial values

    x0, x1, ...

    satisfying

    xn = c1 x(n-1) + c2 x(n-2) + ...

    return the coefficients and initial terms for the sequence
    of cumulative sums.

    e.g. the cumulative sum of Fibonacci numbers starts 0, 1, 2
    and satisfies a(n) = 2*a(n-1)-a(n-3).

    >>> cumulative_lrs([1,1],[0,1])

    ([2, 0, -1], [0, 1, 2])
    """

    nextx = sum(c*x for (c, x) in zip(reversed(cs), xs))
    cumcs = [cs[0]+1] + [c1-c0 for (c0, c1) in
                         pairwise(cs)] + [-cs[-1]]
    cumxs = list(np.cumsum(xs+[nextx]))
    return cumcs, cumxs


def power_lrs(cs, k):
    """
    Return the coefficients of the linear recurrence satisfied
    by the k^{th} powers of a sequence whose coefficients are
    given by cs. (the initial segment of values is not needed)
    """
    if k==0:
        return [1]

    def radical(P):
        return prod(fact for fact, mult in sqf_list(P)[1])

    x, y = symbols('x y')
    P = x**len(cs) - sum(c*x**k for (k, c) in enumerate(reversed(cs)))
    Q = P
    while k-1:
        P = resultant(y**degree(P) * P.subs(x, x/y), Q.subs(x, y), y)
        k -= 1
        P = radical(P)

    power_cs = [expand(P).coeff(x, i) for i in range(degree(P)+1)]
    power_cs = list(reversed([-x for x in power_cs][:-1]))
    return listmap(int, power_cs)


def decimate_lrs(cs, xs, m):
    """
    Return the coefficents of the linear recurrence satisfied
    by the 'every m element' decimations of the linear recurrence
    with the given coefficients cs and initial elements xs.
    """
    d = len(cs)
    #
    # compute all of the terms of the sequence that
    # we will need for the matrix
    #
    f = listtake((d-1)+d*m+1, lrs(cs, xs))
    #
    # now build the d x (d+1) matrix
    #
    M = np.array([[f[d0 + d1*m] for d1 in range(d+1)]
                                for d0 in range(d)], dtype=object)

    new_cs = []
    for i in range(M.shape[1]):
        minor =  sympy.Matrix(np.delete(M, i, axis=1))
        c = (-1)**i * sympy.det(minor) # round(np.linalg.det(minor))
        new_cs.append(c)

    denom = new_cs[-1]
    assert all(new_c % denom == 0 for new_c in new_cs)
    new_cs = np.array(new_cs[:-1], dtype=object)
    new_cs //= denom
    return listmap(int, list(-new_cs)[::-1])


# Here are some examples of using this to get various
# sequences:

#
# The fibonacci numbers
#
# [fastfib(i,10**6,[1,1],[0,1]) for i in range(13)]

#
# Sums of fibonacci numbers
#
# [fastfib(i,10**6,[2,0,-1],[0,1,2]) for i in range(13)]

#
# The squares of fibonacci numbers
#
# [fastfib(i,10**6,[2,2,-1],[0,1,1]) for i in range(13)]

#
# Sums of squares of fibonacci numbers
#
# [fastfib(i,10**6,[2,2,-1],[0,1,2]) for i in range(13)]

#
# Cubes of fibonacci numbers:
#
# [fastfib(i,10**6,[3,6,-3,-1],[0,1,1,8]) for i in range(13)]

#
# Sums of cubes of fibonacci numbers
#
# [fastfib(i,10**6,[4,3,-9,2,1],[0,1,2,10,37]) for i in range(13)]

#
# Sums of cubes of the arithmetic subsequence of fibonacci
# numbers: F11, F35, F59, F83, F107, F131, ...
#
# Here are the coefficents and initial values for the sequence
# of cubes of Fibonacci numbers:
#
# cs = [3, 6,-3,-1]
# xs = [0, 1, 1, 8]
#
# Here we compute the coeffs and initial values
# for the decimated sequence:
#
# xs2 = listtake(100, lrs(cs, xs))[11:84:24]
# cs2 = decimate_lrs(cs, xs, 24)
#
# xs3, cs3 = cumulative_lrs(xs2, cs2)
#

# from math import sqrt, log, ceil

# phi_p = (1.0+sqrt(5.0))/2.0
# phi_m = (1.0-sqrt(5.0))/2.0

# def fib(i):
#     return int((phi_p**(i+1)-phi_m**(i+1))/sqrt(5.0)+0.5)

# def inverse_fib_floor(n):
#     return int(ceil(log(n)/log(phi_p)))

# def zeckendorf(n):
#     a, b = 1, 2
#     fibs = []

#     while a <= n:
#         fibs.append(a)
#         a, b = b, a+b
#     #print(fibs)

#     rep_fibs = []
#     rep_idxs = []
#     i = len(fibs)-1
#     while n != 0:
#         rep_fibs.append(fibs[i])
#         rep_idxs.append(i)
#         n -= fibs[i]
#         #print(i, fibs[i], n)
#         i -= 2
#         while i >= 0 and fibs[i] > n:
#             i -= 1

#     return rep_idxs, rep_fibs


# def fibstr(zs):
#     """
#     If we write

#     x = \sum_{i\in I} F_i

#     as a sum of distinct fibonacci numbers,
#     return the string representation (*'s for
#     the fibonacci numbers included, -'s for
#     those omited.
#     """
#     zstr = ''
#     zprev= 0
#     for z in zs:
#         zstr += '-'*(z-zprev)+'*'
#         zprev=z+1
#     return zstr

