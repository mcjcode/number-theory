#!/usr/bin/env python -i
# -*- coding: utf-8 -*-
"""
General purpose, factorization and modular
arithmetic routines.
"""

import unittest
import math
import numpy as np
from slv import (
    lagrange,
    lagrangeL1,
    )


class SLVTest(unittest.TestCase):


    def test_one(self):
        u = np.array([-1], dtype=object)
        v = np.array([ 3], dtype=object)
        u, v = lagrangeL1(u, v)
        self.assertEqual(u[0], 0)
        self.assertEqual(abs(v[0]), 1)


    def test_FibonacciL1(self):
        a = 1
        b = 2
        for i in range(100):
            print(a, b)
            u = np.array([a], dtype=object)
            v = np.array([b], dtype=object)
            u, v = lagrangeL1(u, v)
            self.assertEqual(u[0], 0)
            self.assertEqual(abs(v[0]), 1)
            a, b = b, a+b
            

    def test_FibonacciL2(self):
        a = 1
        b = 2
        for i in range(100):
            u = np.array([a], dtype=object)
            v = np.array([b], dtype=object)
            u, v = lagrange(u, v)
            self.assertEqual(u[0], 0)
            self.assertEqual(abs(v[0]), 1)            
            a, b = b, a+b


    def test_dim2_L1(self):
        for _ in range(100):
            A = np.array([[1,0],[0,1]], dtype=object)
            for _ in range(20):
                x = np.random.randint(10)
                B = np.array([[1,x],[0,1]], dtype=object)
                A = A.dot(B)
                y = np.random.randint(10)
                B = np.array([[1,0],[y,1]], dtype=object)
                A = A.dot(B)
            self.assertEqual(abs(A[0,0]*A[1,1]-A[1,0]*A[0,1]), 1, f'abs(det({A})) should equal 1')
            u = A[:,0]
            v = A[:,1]
            u, v = lagrangeL1(u, v)
            self.assertEqual(sum(map(abs,u)), 1, f'{A[:,0]} {A[:,1]}. smallest vector {u}, {v} in lattice should have norm 1')


    def test_dim2_L2(self):
        for _ in range(100):
            A = np.array([[1,0],[0,1]], dtype=object)
            for _ in range(20):
                x = np.random.randint(10)
                B = np.array([[1,x],[0,1]], dtype=object)
                A = A.dot(B)
                y = np.random.randint(10)
                B = np.array([[1,0],[y,1]], dtype=object)
                A = A.dot(B)
                
            self.assertEqual(A[0,0]*A[1,1]-A[1,0]*A[0,1], 1, f'det({A}) should equal 1')
            u = A[:,0]
            v = A[:,1]
            u, v = lagrange(u, v)
            self.assertEqual(sum(map(abs,u)), 1, f'smallest vector {u} in lattice should have norm 1.')

            
    def test_minkowski(self):
        for _ in range(100):
            norm = lambda u: math.sqrt(sum(x*x for x in u))
            u = np.random.randint(-100, +100, size=(2,))
            v = np.random.randint(-100, +100, size=(2,))            
            u, v = lagrange(u, v)
            volume = abs(u[0]*v[1] - u[1]*v[0])
            self.assertTrue(norm(u) <= math.sqrt(2*volume))


    def test_minkowskiL1(self):
        for _ in range(100):
            norm = lambda u: sum(map(abs,u))
            u = np.random.randint(-100, +100, size=(2,))
            v = np.random.randint(-100, +100, size=(2,))
            u, v = lagrangeL1(u, v)
            volume = abs(u[0]*v[1] - u[1]*v[0])
            self.assertTrue(norm(u) <= math.sqrt(2*volume), norm(u)/math.sqrt(2*volume))

