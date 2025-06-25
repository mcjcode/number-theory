#!/usr/bin/env python
#
"""
Routines for looping over numbers formed from a
given sequence of primes below a fixed bound.
"""

import unittest
import itertools

from utilities import prod
from bps import bps_facts_w_rep, bps_facts_w_rep_powerful, bpsk


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

    def test_bpsk(self):
        ps = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        for k in range(2, 10):
            cnt1 = 0
            for xs in itertools.combinations(ps, k):
                if prod(xs) <= 100:
                    cnt1 += 1
            cnt2 = len(list(bpsk(ps, 100, k)))
            self.assertEqual(cnt1, cnt2)