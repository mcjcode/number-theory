#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from lucy import (
    sievecnt,
    sievecnt_mod3,
    sievecnt_mod4,
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
            
