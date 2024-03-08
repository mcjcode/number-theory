#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from lucy import (
    sievecnt,
    sievesum,
    sievecnt_mod3,
    sievecnt_mod4,
    sievesum_mod3,
)


class LucyTest(unittest.TestCase):
    answer_key = {10**1:4,
                  10**2:25,
                  10**3:168,
                  10**4:1229,
                  10**5:9592,
                  10**6:78498}
    
    def test_sievecnt(self):
        for k, v in self.answer_key.items():
            V, S = sievecnt(k)
            self.assertEqual(S[k], v)

    def test_sievecnt_mod4(self):
        for k, v in self.answer_key.items():
            V, S1, S3 = sievecnt_mod4(k)
            self.assertEqual(S1[k]+S3[k]+1, v)

    def test_sievecnt_mod3(self):
        for k, v in self.answer_key.items():
            V, S1, S2 = sievecnt_mod3(k)
            self.assertEqual(S1[k]+S2[k]+1, v)

    def test_sievesum_mod3(self):
        """
        Check that 'sum of the primes <= n congruent to 1(mod3) and 2(mod3)'
        gives the same answer (less 1, for the prime '3') as the just the
        sum of the primes <= n.
        """
        for n in range(4, 101):
            V, Stot = sievesum(n)
            S = sievesum_mod3(n)
            self.assertEqual(S[1][n]+S[2][n]+3, Stot[n])
