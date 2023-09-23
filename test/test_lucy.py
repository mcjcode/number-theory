#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from lucy import (
    sievecnt,
)


class LucyTest(unittest.TestCase):
    def test_sievecnt(self):
        answer_key = {10:4,
                      100:25,
                      1000:168,
                      10000:1229,
                      100000:9592,
                      1000000:78498}
        for k, v in answer_key.items():
            V, S = sievecnt(k)
            self.assertEqual(S[k], v)

