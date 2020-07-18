#!/usr/bin/env python
#
"""
Routines for counting partitions.
"""

import numpy as np

from utilities import timeit


@timeit
def npartitions(n):
    """
    :param n: the upper bound
    :return: number of unrestricted partitions of k, for k in [0..n]

    Uses Euler's pentagonal theorem.
    """

    def pentagonals():
        k = 1
        while True:
            yield k * (3 * k - 1) // 2
            yield k * (3 * k + 1) // 2
            k += 1

    pents = []
    pent_gen = pentagonals()
    p = next(pent_gen)
    while p <= n:
        pents.append(p)
        p = next(pent_gen)

    arr = [0]*(n+1)
    arr[0] = 1
    for j in range(1, n+1):
        p = 0
        sgn = -1
        cnt = 1
        for k in pents:
            if k > j:
                break
            p -= sgn*arr[j-k]
            cnt = 1-cnt
            if cnt:
                sgn *= -1
        arr[j] = p
    return arr


@timeit
def npartitions_distinct(n):
    """
    :param n: the upper bound
    :return: number of partitions of k into distinct parts, for k in [0..n]

    Note that this is *also* the number of unrestricted partitions of
    k into odd parts.

    cf: http://home.dimacs.rutgers.edu/~asills/talks/PartitionRecurrences.pdf
    """

    parr = npartitions(2*n)
    qarr = [0]*(n+1)
    qarr[0] = 1
    for j in range(1, n+1):
        j2 = 2*j
        k = j2
        sgn = 1
        cnt = 1
        k0 = 0
        p = 0
        while k >= 0:
            p += sgn*parr[k]
            k0 += 1
            k -= k0
            if cnt:
                sgn *= -1
            cnt = 1 - cnt
        qarr[j] = p
    return qarr


def npartitions_distinct_mod(n, modulus):
    """
    Count the number of partitions of k into distinct
    pieces, for all k in [0..n]

    Based on recursion proven by J. Ewell:

    q0(n) - 2*q(n- 1) + 2*q(n- 4)
          - 2*q(n- 9) + 2*q(n-16)
          - 2*q(n-25) + 2*q(n-36)
          - ...

          = (-1)**m  if n=m*(3*m-1)//2 (n is pentagonal
                  0  otherwise
    """
    # build the array of squares
    # 1, 4, 9, 16, ...
    squares = []
    k = 1
    tn = k * k
    while tn <= n:
        squares.append(tn)
        k += 1
        tn = k * k

    # build the array of 2x pentagonal
    # numbers
    pents = []
    m = 1
    pent = m * (3 * m - 1) // 2
    while pent <= n:
        pents.append((pent, (-1) ** (m % 2)))
        pent = m * (3 * m + 1) // 2
        if pent <= n:
            pents.append((pent, (-1) ** (m % 2)))
        else:
            break
        m += 1
        pent = m * (3 * m - 1) // 2
    npents = len(pents)

    arr = np.array([0] * (n + 1), dtype=np.int32)
    arr[0] = 1

    penti = 0
    for n1 in range(1, n + 1):
        while penti < npents and pents[penti][0] < n1:
            penti += 1
        if penti < npents and n1 == pents[penti][0]:
            n1_is_pent = True
            pent_sgn = pents[penti][1]
        else:
            n1_is_pent = False

        k = 0
        accum = 0
        for tr in squares:
            if tr > n1: break
            term = (-1) ** (tr % 2) * arr[n1 - tr]
            accum = (accum - 2 * (-1) ** (tr % 2) * arr[n1 - tr]) % modulus
        if n1_is_pent:
            accum = (accum + pent_sgn) % modulus
            penti += 1

        arr[n1] = accum

    return arr


def npartitions_distinct_odd_mod(n, modulus):
    """
    Count the number of partitions of k into distinct odd
    pieces, for all k in [0..n]

    Based on recursion proven by J. Ewell:

    q0(n) - q0(n-1) - q0(n-3)
          + q0(n-6) + q0(n-10)
          - q0(n-15) - q0(n-21)
          + ...

          = (-1)**m  if n=m*(3*m-1)
                  0  otherwise
    """
    # build the array of triangle numbers
    # 1, 3, 6, 10, ..,
    triangles = []
    tn = 1
    k = 1
    while tn <= n:
        triangles.append(tn)
        k += 1
        tn += k
    # print triangles

    # build the array of double pentagonal
    # numbers
    pent2s = []
    m = 1
    pent2 = m * (3 * m - 1)
    while pent2 <= n:
        pent2s.append((pent2, (-1) ** (m % 2)))
        pent2 = m * (3 * m + 1)
        if pent2 <= n:
            pent2s.append((pent2, (-1) ** (m % 2)))
        else:
            break
        m += 1
        pent2 = m * (3 * m - 1)

    arr = np.array([0] * (n + 1), dtype=np.int32)
    arr[0] = 1
    npent2s = len(pent2s)

    pent2i = 0
    for n1 in range(1, n + 1):

        while pent2i < npent2s and pent2s[pent2i][0] < n1:
            pent2i += 1
        if pent2i < npent2s and n1 == pent2s[pent2i][0]:
            n1_is_2pent = True
            pent_sgn = pent2s[pent2i][1]
        else:
            n1_is_2pent = False

        k = 0
        accum = 0
        for tr in triangles:
            if tr > n1: break
            accum -= (-1) ** (tr % 2) * arr[n1 - tr]

        if n1_is_2pent:
            accum += pent_sgn
            pent2i += 1

        arr[n1] = accum

    return arr
