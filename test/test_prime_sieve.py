import unittest
from utilities import prod
from prime_sieve import (
    factorizations,
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
            n2 = prod(p**e for p, e in facts[n].items())
            self.assertEqual(n, n2)
