#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

import numpy as np

from utilities import (
    factorize,
    gcd,
    multiplicities,
)


def phi(nn):
    """
    :param nn: a positive integer
    :return: the number of natural numbers 1 <= kk < nn
             that are relatively prime to nn.
    """
    if nn == 1:
        return 0
    factors = factorize(nn)
    primes, powers = multiplicities(factors)
    return np.prod([pp**(kk-1)*(pp-1) for (pp,kk) in zip(primes,powers)])


class PhiTest(unittest.TestCase):
    def test_phi(self):
        for nn in range(1,100):
            nresid = sum(1 for kk in range(1,nn) if gcd(nn,kk) == 1)
            self.assertEqual(nresid, phi(nn))
