#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from finite_field import FiniteField


class FiniteFieldTest(unittest.TestCase):
    def runTest(self):
        for p in [7]:
            for e in [2]:
                ffield = FiniteField(p, 2)
                for elem in list(ffield):
                    self.assertEqual(elem, ffield.one() * elem, msg='1*X == X')
                    self.assertEqual(elem, elem * ffield.one(), msg='X*1 == X')
                    self.assertEqual(ffield.zero() * elem, ffield.zero(), msg='0*X==0')
                    self.assertEqual(elem * ffield.zero(), ffield.zero(), msg='X*0==0')
                    self.assertEqual(elem ** 2, elem * elem, msg='X**2 == X*X')
                    if elem != ffield.zero():
                        self.assertEqual(elem**ffield.order(elem), ffield.one(), msg='X**ord(X)==1')
