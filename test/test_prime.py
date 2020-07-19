#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from prime import (
    primesum_naive,
    primesum_naive2,
)

sum2mill = 142913828922


class PrimeSumTest(unittest.TestCase):
    def test_primesum_naive(self):
        self.assertEqual(primesum_naive(2*10**6), sum2mill)

    def test_primesum_naive2(self):
        self.assertEqual(primesum_naive2(2*10**6), sum2mill)
