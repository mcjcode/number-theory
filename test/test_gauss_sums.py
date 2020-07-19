#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from utilities import squarefree
from gauss_sums import ideal_class_number


class GaussSumTest(unittest.TestCase):
    def test_1(self):
        class_number = ideal_class_number(1)
        self.assertEqual(class_number, 1, 'Class number of Q should be 1.  Was %d.' % class_number)

    def test_2(self):
        class_number = ideal_class_number(-1)
        self.assertEqual(class_number, 1, 'Class number of Q[i] should be 1.  Was %d.' % class_number)

    @staticmethod
    def test_ideal_class_number():
        for d in range(1, 500):
            if squarefree(d):
                _ = ideal_class_number(d)
