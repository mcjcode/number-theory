#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from imaginary_complex_fields import (
    all_reduced_forms,
    proper_reduced_form,
    sum_sq_rep,
)


class QuadraticFormTests(unittest.TestCase):
    def test_one(self):
        for disc in range(-4, -50, -4):
            for form in all_reduced_forms(disc):
                self.assertEqual(proper_reduced_form(*form), form)


class SumOfSquaresAlgoTest(unittest.TestCase):
    def test_one(self):
        ps = [5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137, 149]
        for p in ps:
            a, b = sum_sq_rep(p)
            self.assertEqual(p, a*a+b*b)
