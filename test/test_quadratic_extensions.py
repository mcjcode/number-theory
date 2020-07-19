#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

from prime_sieve import segmented_sieve
from quadratic_extensions import (
    tonelli_shanks,
    legendre,
    legendre_ch,
)


def _ch12(a):
    """
    # the mod-12 character that identifies for which
    # primes p not dividing 12 is 3 a quadratic
    # residue.

    :param a:
    :return: 0/1/-1
    """
    # ----- 0.  1.  2.  3.  4.  5.  6.  7.  8.  9. 10  11
    return [0, +1,  0,  0,  0, -1,  0, -1,  0,  0,  0, +1][a % 12]


class TonelliShanksTest(unittest.TestCase):
    def test_1(self):
        for p in segmented_sieve(1000):
            for a in range(p):
                if legendre(a, p) == +1:
                    r = tonelli_shanks(a, p)
                    self.assertEqual(r*r % p, a)


class LegendreCharacterTest(unittest.TestCase):
    def test_12(self):
        ch = legendre_ch(3)
        for xx in range(12):
            self.assertEqual(ch(xx), _ch12(xx), 'mod 3*4 Legendre character incorrect')
