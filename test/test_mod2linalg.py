import unittest

import random
import numpy as np
import mod2linalg


class SolveTest(unittest.TestCase):
    def runTest(self):
        for _ in range(100):
            N = random.randint(10, 20)
            A = np.random.randint(0, 2, (N, N), dtype=np.int8)
            if mod2linalg.det_mod_2(A)==1:
                b = np.random.randint(0, 2, (N,), dtype=np.int8)
                v = mod2linalg.solve_linear_equations_mod_2(A, b)
                self.assertTrue(all((A.dot(v)%2)==b))   
            
    
class NullSpace(unittest.TestCase):
    def runTest(self):
        for _ in range(100):
            M = random.randint(10, 20)
            N = M + random.randint(0, 7) - 3
            A = np.random.randint(0, 2, (M, N), dtype=np.int8)
            nullspace = mod2linalg.null_space_mod_2(A)
            print(A)
            print(nullspace)
            print(A.dot(nullspace) % 2)
            for v in A:
                for w in nullspace.transpose():
                    self.assertEqual(v.dot(w)%2, 0)
