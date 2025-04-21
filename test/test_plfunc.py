import unittest

from fractions import Fraction
from plfuncs import PLFunc

class PLFuncTest(unittest.TestCase):

    def test_prune_constant_length(self):
        N = 10
        xs = [Fraction(i, N) for i in range(N+1)]
        ys = [Fraction(1, 1)]*(N+1)
        f = PLFunc(xs, ys)
        g = f.prune()
        self.assertEqual(len(g.xs), 2, 'pruned constant function should only have 2 nodes, has %s' % len(g.xs))

    def test_prune_linear_length(self):
        N = 10
        m = 3
        b = 7
        xs = [Fraction(i, N) for i in range(N+1)]
        ys = [m*x + b for x in xs]
        f = PLFunc(xs, ys)
        g = f.prune()
        self.assertEqual(len(g.xs), 2, 'pruned linear function should only have 2 nodes, has %s' % len(g.xs))

    def test_prune_quadratic_length(self):
        N = 10
        xs = [Fraction(i, N) for i in range(N+1)]
        ys = [((2*x-1)/2)**2 for x in xs]
        f = PLFunc(xs, ys)
        g = f.prune()
        self.assertEqual(len(g.xs), N+1, 'pruned quadratic function should have as many nodes as original, has %s' % len(g.xs))

