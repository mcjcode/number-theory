#!/usr/bin/env python
#
"""
Routines for looping over numbers formed from a
given sequence of primes below a fixed bound.
"""

import unittest
from utilities import sqrtInt, prod


def bps_w_rep(n, ps, i=0, only_powerful=False):
    """
    Return all numbers in [1..n] whose prime
    factors are all in the list ps[i:] (repetitions
    allowed)

    If 'only_powerful' is set to True, then
    only return the numbers that are powerful.

    Recursive, so fails if recursion limit
    is not sufficient.
    """
    if len(ps) == i:
        yield 1
    else:
        p = ps[i]
        for pr in bps_w_rep(n, ps, i+1, only_powerful=only_powerful):
            yield pr
        if only_powerful:
            pk = p*p
        else:
            pk = p
        while pk <= n:
            for pr in bps_w_rep(n//pk, ps, i+1, only_powerful=only_powerful):
                yield pk*pr
            pk *= p


def bps_facts_w_rep_recur(n, ps, i=0):
    """
    Return all factorizations of numbers in [1..n]
    whose prime factors are all in the list ps[i:].
    """
    if len(ps) == i or ps[i] > n:
        yield []
    else:
        p = ps[i]
        pk = 1
        e = 0
        while pk <= n:
            for pr in bps_facts_w_rep_recur(n//pk, ps, i+1):
                if e == 0:
                    yield pr
                else:
                    yield [(p, e)] + pr
            pk *= p
            e += 1


def bps_facts_w_rep_powerful_recur(n, ps, i=0):
    """
    Return all factorizations of powerful numbers in [1..n] whose prime
    factors are all in the list ps[i:]
    """
    if len(ps) == i:
        yield []
        return
    p = ps[i]
    if p**2 > n:
        yield []
        return
    pk = 1
    e = 0
    while pk <= n:
        for pr in bps_facts_w_rep_powerful_recur(n//pk, ps, i+1):
            if e == 0:
                yield pr
            else:
                yield [(p, e)] + pr
        if e == 0:
            pk = p*p
            e = 2
        else:
            pk *= p
            e += 1


def bps_facts_w_rep_powerful2(n, ps, i, sofar):
    """
    Set sofar to []
    """
    if len(ps) == i:  # we're out of primes
        yield sofar
        return
    p = ps[i]
    pk = 1
    e = 0
    while pk <= n:
        for pr in bps_facts_w_rep_powerful2(n//pk, ps, i+1, sofar+([(p, i, e, n//pk)] if e else [])):
            yield pr
        pk *= p
        e += 1


def bps_facts_w_rep(n, ps):
    """
    n: the integer upper bound
    ps: the list of primes

    yields: all factorized integers [(p1, e1), ..., (pk, ek)]
        where p1**e1 * ... * pk**ek <= n and all pi are
        from the list ps

    Note: the implementation is iterative, and not recursive,
    and so is suitable for cases where ps is large (i.e. cases
    where the stack would normally grow to the size of ps)
    """
    num_primes = len(ps)
    if n < 1:
        return
    yield []
    val = [[ps[-1], num_primes-1, 1, n//ps[-1]]]
    while True:
        yield [(v[0], v[2]) for v in val]
        #
        # find the largest prime that is
        # 1) as large or larger than the
        #    primes already in the factorization, and
        # 2) less than what is left
        #
        pi = num_primes-1
        p = ps[pi]
        while p > val[-1][3]:
            if p < val[-1][0]:
                break
            pi -= 1
            if pi < 0:
                break
            p = ps[pi]

        if pi < 0:
            return

        #
        # if there is no such prime, we have to back
        # track over the previous prime used.
        #
        if p < val[-1][0]:  # if there is not an add'l prime
            pi = val[-1][1]
            p = ps[pi-1]
            if len(val) == 1:
                val[-1] = [p, pi-1, 1, n//p]
            else:
                del val[-1]
                if p == val[-1][0]:
                    val[-1][2] += 1
                    val[-1][3] //= p
                else:
                    val.append([p, pi-1, 1, val[-1][3]//p])
        elif p == val[-1][0]:  # we can divide by our prime again
            val[-1][2] += 1
            val[-1][3] //= p
        else:
            val.append([p, pi, 1, val[-1][3]//p])


def bps_facts_w_rep_powerful(n, ps):
    """
    n: the integer upper bound
    ps: the list of primes

    yields: all factorized integers [(p1, e1), ..., (pk, ek)]
        where p1**e1 * ... * pk**ek <= n and all pi are
        from the list ps

    Note: the implementation is iterative, and not recursive,
    and so is suitable for cases where ps is large (i.e. cases
    where the stack would normally grow to the size of ps)
    """
    num_primes = len(ps)
    if n < 1:
        return
    yield []
    val = [[ps[-1], num_primes - 1, 2, n // (ps[-1]**2)]]
    while True:
        yield [(v[0], v[2]) for v in val]
        #
        # find the largest prime that is
        # 1) as large or larger than the
        #    primes already in the factorization, and
        # 2) less than what is left
        #
        pi = num_primes - 1
        p = ps[pi]
        last_val = val[-1]
        last_p = last_val[0]
        last_pi = last_val[1]
        last_n = last_val[3]
        sqrt_last_n = sqrtInt(last_n)
        lb = max(last_p, sqrt_last_n)
        if last_p > last_n:  # sqrt_last_n:
            pi = last_pi-1
            if pi >= 0:
                p = ps[pi]
        else:
            while p > lb:
                pi -= 1
                if pi < 0:
                    break
                p = ps[pi]
            if p == last_p and p > last_n:
                pi -= 1
                if pi >= 0:
                    p = ps[pi]
        if pi < 0:
            return

        #
        # if there is no such prime, we have to back
        # track over the previous prime used.
        #
        if p < last_p:  # if there is not an add'l prime
            pi = last_val[1]-1
            p = ps[pi]
            if len(val) == 1:
                last_val[:] = [p, pi, 2, n // (p*p)]
            else:
                del val[-1]
                last_val = val[-1]
                if p == last_val[0]:
                    last_val[2] += 1
                    last_val[3] //= p
                else:
                    val.append([p, pi, 2, last_val[3] // (p*p)])
        elif p == last_p:  # we can divide by our prime again
            last_val[2] += 1
            last_val[3] //= p
        else:
            val.append([p, pi, 2, last_n//(p*p)])


class BPSTest(unittest.TestCase):

    @staticmethod
    def f(facts):
        return prod([p**e for (p, e) in facts])

    def test_bps_facts_w_rep(self):

        lst1 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
        lst2 = set(map(self.f, bps_facts_w_rep(10, [2, 3, 5, 7])))
        self.assertEqual(lst1, lst2)

    def test_bps_facts_w_rep_powerful(self):
        lst1 = {1, 4, 8, 9}
        lst2 = set(map(self.f, bps_facts_w_rep_powerful(10, [2, 3])))
        self.assertEqual(lst1, lst2)
