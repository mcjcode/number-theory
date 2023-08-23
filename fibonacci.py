import numpy as np
import sympy
from utilities import fastpow


def _matrix_modpow(a, k, p):
    one = np.identity(a.shape[0], dtype=object)
    mul = lambda x, y: x.dot(y) % p
    return fastpow(a, k, mul, one)


def tribonacci(a=0,b=0,c=1,modulus=2):
    while True:
        yield a
        a, b, c = b, c, (a+b+c) % modulus


def tribonacci2(modulus):
    a, b, c = 0, 0, 1
    while True:
        yield a, b
        d = (a+b+c)%modulus
        e = (b+c+d)%modulus
        a, b, c = c, d, e


def fastfib(n, modulus, coeffs=[1, 1], init=[0,1]):
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
    I = np.identity(order-1, dtype=object)
    Z = np.zeros((order-1,1), dtype=object)
    C = np.array([coeffs], dtype=object)    
    F = np.block([[Z,I],
                  [ C ]])
    return _matrix_modpow(F, n, modulus).dot(np.array(init).transpose())[0] % modulus


def fastfib_decimate(n, modulus, coeffs=[1, 1], init=[0, 1], m=1, i=0):
    """
    Return the n'th element of the sequence
    
    F(i), F(m+i), F(2*m+i), ...
    
    where F is specified by the given coefficients and initial values.
    """
    coeffs = coeffs[::-1]
    order = len(coeffs)
    I = np.identity(order-1, dtype=object)
    Z = np.zeros((order-1,1), dtype=object)
    C = np.array([coeffs], dtype=object)
    F = np.block([[Z,I],
                  [ C ]])
    
    init = [fastfib(i+j, modulus, coeffs[::-1], init) for j in range(len(coeffs))]
    F = _matrix_modpow(F, m, modulus)
    
    return _matrix_modpow(F, n, modulus).dot(np.array(init).transpose())[0] % modulus


def fastfib_sum1(n, m):
    """return the sum of the first n fibonacci numbers, mod m"""
    return fastfib(n+1, m)


def fastfib_sum2(n, m):
    """return the sum of the squares of the first n fibonacci numbers, mod m"""
    return (fastfib(n-1, m) * fastfib(n, m)) % m


def fastfib_sum3(n, m):
    """return the sum of the cubes of the first n fibonacci numbers, mod m"""
    if n<=1:
        return 0
    else:
        return (fastfib(3*n-1,m)+(-1)**n*6*fastfib(n-2,m)+5)//10

    
def find_recursion(xs, deg):
    """
    Find a linear recursion of degree deg that is satisfied
    by the sequence xs.
    """
    amat = [ xs[i  :deg+i  ] for i in range(deg)]
    A = sympy.Matrix(amat)
    bmat = [ xs[i+1:deg+i+1] for i in range(deg)]
    B = sympy.Matrix(bmat)
    return list((B*A.inv())[-1,:])[::-1]


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
    
    ([2, 0, -1], array([0, 1, 2]))
    """
    
    nextx = sum( c*x for (c,x) in zip(cs[::-1],xs) )
    cumcs = [cs[0]+1] + [c1-c0 for (c1,c0) in zip(cs[:-1], cs[1:])] + [-cs[-1]]
    cumxs = np.cumsum(xs+[nextx])
    return cumcs, cumxs



# Here are some examples of using this to get various
# sequences:

#
# The fibonacci numbers
#
[fastfib(i,10**6,[1,1],[0,1]) for i in range(13)]

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


# def fibstr(zz):
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
#     for z in zz:
#         zstr += '-'*(z-zprev)+'*'
#         zprev=z+1
#     return zstr
