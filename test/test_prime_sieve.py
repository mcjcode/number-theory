import unittest
from utilities import prod
from prime_sieve import (
    factorizations,
    quadratic_sieve,
    segmented_sieve,
)


class PrimeSieveTest(unittest.TestCase):
    def test_factorizations(self):
        """
        Check that the returned list of Counters
        are actually factorizations of their respective
        numbers
        """
        N = 10_000
        facts = factorizations(N)
        for n in range(1, N+1):
            with self.subTest(n=n):
                n2 = prod(p**e for p, e in facts[n].items())
                self.assertEqual(n, n2)


class QuadraticSieveTest(unittest.TestCase):
    def test_factorizations(self):
        N = 10_000
        facts = quadratic_sieve(N)
        for n in range(1, N+1):
            with self.subTest(n=n):
                n2 = prod(p**e for p, e in facts[n])
                self.assertEqual(n**2+1, n2)

def test_segmented_sieve():
    ps = [2,3,5,7]
    for n in range(10):
        assert list(segmented_sieve(n))==[p for p in ps if p<=n]
        