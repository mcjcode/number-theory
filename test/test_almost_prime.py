from utilities import sqrtInt, cbrtInt, factorize2
from prime_sieve import segmented_sieve
from lucy import sievecnt

from almost_prime import (
    biprimes,
    triprimes,
    quadprimes,
    almost_primes,
    pi2,
    pqr,
    p2qr,
    p2q2r,
    pkqr,
    pqrs,
)


def almost_primes0(k, N):
    """
    Reference implementation for testing
    """
    retval = 0
    for n in range(1, N+1):
        f = list(factorize2(n))
        if sum(e for _, e in f)==k:
            retval += 1
    return retval


def p2qr0(N):
    retval = 0
    for n in range(1, N+1):
        f = sorted([e for p, e in factorize2(n)])
        if f==[1,1,2]:
            retval += 1
    return retval


def p2q2r0(N):
    retval = 0
    for n in range(1, N+1):
        f = sorted([e for p, e in factorize2(n)])
        if f==[1,2,2]:
            retval += 1
    return retval


def factorization_form_cnt0(form, N):
    retval = 0
    for n in range(1, N+1):
        f = sorted([e for p, e in factorize2(n)])
        if f==form:
            retval += 1
    return retval


def test_almost_primes():
    for n in range(1, 1001):
        for k in range(1, 10):
            assert almost_primes(k, n)==almost_primes0(k, n)

def test_pqr():
    for n in range(1, 1001):
        ps = list(segmented_sieve(sqrtInt(n)))
        V, S0 = sievecnt(n)
        ans1 = factorization_form_cnt0([1, 1, 1], n)
        ans2 = pqr(n, ps, V, S0)
        assert ans1==ans2


def test_pqrs():
    for n in range(1, 1001):
        ps = list(segmented_sieve(sqrtInt(n)))
        V, S0 = sievecnt(n)
        ans1 = factorization_form_cnt0([1, 1, 1, 1], n)
        ans2 = pqrs(n, ps, V, S0)
        assert ans1==ans2

def test_pi2():
    