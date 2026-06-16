"""
Routines concerning finding near integer nth roots or
testing for the existence of exact integer nth roots.
"""

import math


def sqrtInt(a: int) -> int:
    """
    :param a: a non-negative integer
    :return: the largest integer smaller than the square root of `a`
    """
    if a==0:
        return 0
    n = 1 << (a.bit_length()+1)//2
    n0 = (n + a//n)//2
    while n0 < n:
        n = n0
        n0 = (n + a//n)//2
    return n


def cbrtInt(a: int) -> int:
    """
    :param a: a non-negative integer
    :return: the largest integer smaller than the cube root of `a`
    """
    if a==0:
        return 0
    n = a >> (2*a.bit_length())//3
    q, r = divmod(a, n*n)
    xi = (2*n + (q+(r>0)))//3
    prev = -1
    while xi!=n and xi!=prev:
        n, prev = xi, n
        q, r = divmod(a, n*n)
        xi = (2*n + (q+(r>0)))//3
    if xi**3>a:
        xi -= 1
    return xi


def nrtInt(n: int, a: int) -> int:
    """
    :param n: a non-negative integer
    :param a: a non-negative integer
    :return: the largest integer z such that z**`n` <= `a`
    """
    if a==0:
        return 0

    ans = int(a**(1/n))
    while (ans+1)**n <= a:
        ans += 1
    return ans


_issq_pe = [p**int(math.log(100)/math.log(p)) for p in [2, 3, 5, 7, 11, 13]]
_sqlowbits = {}
for pe in _issq_pe:
    # _pmax = p**e
    bits = [True]*pe
    for x in range(pe):
        bits[(x*x)%pe] = False
    _sqlowbits[pe] = bits


def issq(n: int) -> bool:
    # first check that n is plausibly a square by
    # checking that it is a square mod p**e for
    # small p**e
    if n<0:
        return False
    for pe, bits in _sqlowbits.items():
        if bits[n%pe]:
            return False
    return sqrtInt(n)**2==n
