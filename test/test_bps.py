#!/usr/bin/env python
#
"""
Routines for looping over numbers formed from a
given sequence of primes below a fixed bound.
"""

import unittest
import itertools

from utilities import prod, factorize2

from bps import (
    bps,
    bps_w_rep,
    bps_w_sign,
    bps_facts_w_rep,
    bpsk,
)


class BPSTest(unittest.TestCase):

    @staticmethod
    def f(facts):
        return prod([p**e for (p, e) in facts])

    def test_bps(self):
        #
        # all products of subsets of ps
        # are <=prod(ps).
        #
        ps = [2, 3, 5, 7, 11]
        N = prod(ps)
        assert sum(1 for _ in bps(N, ps))==2**len(ps)
        
    def test_bps_facts_w_rep(self):

        lst1 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
        lst2 = set(map(self.f, bps_facts_w_rep(10, [2, 3, 5, 7])))
        self.assertEqual(lst1, lst2)

    def test_bpsk(self):
        ps = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        for k in range(2, 10):
            cnt1 = 0
            for xs in itertools.combinations(ps, k):
                if prod(xs) <= 100:
                    cnt1 += 1
            cnt2 = len(list(bpsk(ps, 100, k)))
            self.assertEqual(cnt1, cnt2)

    def test_bps_w_sign(self):
        ps = [2, 3, 5, 7, 11, 13, 17, 19]
        N = 1000
        assert len(list(bps_w_sign(ps, N)))==len(list(bps(N, ps)))

    def test_bps_w_sign2(self):
        ps = [2, 3, 5, 7, 11, 13, 17, 19]
        N = 1000
        for s, x in bps_w_sign(ps, N):
            f = list(factorize2(x))
            assert all(e==1 for p, e in f)
            assert s==pow(-1, len(f))
            
    def test_bps_w_rep(self):
        ps = [2, 3, 5, 7]
        N = 10
        assert len(list(bps_w_rep(N, ps)))==N

    def test_bps_w_rep_powerful(self):
        ps = [2, 3, 5, 7]
        N = 1000
        for n in bps_w_rep(N, ps, 0, only_powerful=True):
            for p, e in factorize2(n):
                assert e>=2

    def test_bps_facts_w_rep(self):
        ps = [2, 3, 5, 7]
        N = 1000
        set1 = set(prod(p**e for p, e in f) for f in bps_facts_w_rep(N, ps))
        set2 = set(bps_w_rep(N, ps))
        assert set1==set2

