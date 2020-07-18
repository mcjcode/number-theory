#!/usr/bin/env python
# -*- coding: utf-8 -*-


def farey(n):
    """
    :param n: a positive integer
    :return: yields the farey sequence of order n.

    The farey sequence of order n is all reduced
    fractions between 0 and 1, in reduced terms, with
    denominator <= n, in increasing order.

    It is also a convenient way to obtain a sequence of
    all pairs (a, b) with a <= b <= n with gcd(a, b)=1.

    An O(n) time algorithm with O(1) memory footprint. 
    """
    (a, b, c, d) = (0, 1, 1, n)
    yield a, b
    while c <= n:
        k = (n+b) // d
        a, b, c, d = c, d, k*c-a, k*d-b
        yield a, b
