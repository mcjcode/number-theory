import unittest
import math

from roots import sqrtInt, cbrtInt, nrtInt, issq, cardano


class RootsTest(unittest.TestCase):

    def test_sqrtInt(self):
        for n in range(1, 1000):
            n2 = n*n
            self.assertEqual(sqrtInt(n2-1), n-1)
            self.assertEqual(sqrtInt(n2), n)
            self.assertEqual(sqrtInt(n2+1), n)

    def test_cbrtInt(self):
        for n in range(1000):
            cbrtn = cbrtInt(n)
            self.assertTrue(cbrtn**3 <= n < (cbrtn+1)**3)

    def test_nrtInt(self):
        for n in range(100000):
            for k in range(4, 7):
                nrtn = nrtInt(k, n)
                assert nrtn**k <= n < (nrtn+1)**k

    def test_issq(self):
        for n in range(1, 10000):
            self.assertTrue(issq(n*n))
            self.assertTrue(not issq(n*n+1))

    def test_cardano(self):
        for x in range(1, 10):
            self.assertAlmostEqual(cardano(0.0, -x), math.cbrt(x), places=10)
