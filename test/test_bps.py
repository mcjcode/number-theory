#!/usr/bin/env python
#
"""
Routines for looping over numbers formed from a
given sequence of primes below a fixed bound.
"""

import unittest

from utilities import prod
from bps import bps_facts_w_rep, bps_facts_w_rep_powerful


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
