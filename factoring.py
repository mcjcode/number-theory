import math; from math import gcd, log
from fibonacci import fastfib
from utilities import prod, sqrtInt, issq
from prime_sieve import segmented_sieve


def fermat_factor(n: int, a=1, b=1) -> tuple:
    if (a*b)%2==0:
        raise ValueError('a or b is even')
    
    if gcd(a, b)>1:
        raise ValueError('a and b are not relatively prime')

    if n%2==0:
        return 2, n//2
    
    g = gcd(a*b, n)
    if g>1:
        return g, n//g

    abn = a*b*n
    x = sqrtInt(abn)
    if x*x<abn:
        x += 1
    z = x*x - abn
    while True:
        if issq(z):
            break
        z += 2*x+1
        x += 1

    y = sqrtInt(z)

    assert (x-y)*(x+y) == a*b*n

    u, v = x-y, x+y
    # now clear the factors of a and b from u and v
    
    g = gcd(a, u)
    u //= g
    a //= g

    g = gcd(b, u)
    u //= g
    b //= g

    g = gcd(a, v)
    v //= g
    a //= g

    g = gcd(b, v)
    v //= g
    b //= g

    assert u*v == n
    
    return u, v


def pollard_p_minus_1(n, B, ps=[2, 3, 5]):
    a = 2
    for p in ps:
        if p>B:
            break
        k = int(log(B)/log(p))
        a = pow(a, p**k, n)
        g = gcd(a-1, n)
        if 1<g<n:
            return g
    return False


def williams_p_plus_1(N: int, B: int, A: int) -> int:
    """
    If N has a prime factor p such that p+1 is
    smooth, then p will divide the numbers V_n - 2
    where n is a multiple of p+(D|p) where V is a
    Lucas sequence and D=A**2-4

    :param N: the number we are attempting to factor
    :param B: the smoothness bound (we will check
              the gcd of N and the terms V_{1}, V_{2!},
              ... V_{B!}.
    :param A: a parameter on which the Lucas sequence
              depends.  we want to choose A such that
              (D|p)=1, but we don't know the p yet, so
              we may have to try several values of A.
    :return:  either a positive integer which is a
              proper factor of N, or 0, if no such
              factor was found.

    N = 2**67-1 # = 147573952589676412927
    for A in 9, 17, 23:
        result = williams(N, B=10000, A)
        if result and result!=N:
            print(result)
            break

    will find the factor 193707721 of M67, a Mersene number
    that Mersene reported (incorrectly) as prime, and whose
    factors were presented by Frank Cole in 1903.
    """
    D = A**2-4
    n = 1
    for i in range(1, B+1):
        n *= i
        v = fastfib(n, N, coeffs=[A, -1], init=[2, A])
        g = gcd(N, v-2)
        if g>1 and g!=N:
            return g
    return 0
