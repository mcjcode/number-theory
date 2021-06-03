#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest

from real_quadratic_fields import (
    pell,
)


class PellTest(unittest.TestCase):

    def test_pell(self):
        for d in range(2, 10000):
            a, b = pell(d**4+1)
            self.assertEqual(a*a-(d**4+1)*b*b, 1, 'pell equn not satisfied for d=%d' % (d**4+1, ))
