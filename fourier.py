#
# fourier.py
#

import numpy as np
from utilities import primitive_root
from factoring import factorize


def fourierZ6():
    """
    return a diagonalizing basis for the
    the regular representation of (Z/6Z,*).  I.e.
    the columns of the returned matrix represent
    the elements

    v0 =  g0
    v1 = -g0 + (g2+g4)/2
    v2 =        g4-g2
    v3 = -g0 +            g3
    v4 = (g1+g2-g4-g5)/2
    v5 =  g0 - g3 + (g1-g2-g4+g5)/2

    These have the property that

    vi*vj = 0  (i!=j)
    vi*vi = vi

    """
    #                v0   v1   v2   v3   v4   v5
    return np.array([[1,  -1,   0,  -1,   0,   1],   # g0
                     [0,   0,   0,   0, 1/2, 1/2],   # g1
                     [0, 1/2,-1/2,   0, 1/2,-1/2],   # g2
                     [0,   0,   0,   1,   0,  -1],   # g3
                     [0, 1/2, 1/2,   0,-1/2,-1/2],   # g4
                     [0,   0,   0,   0,-1/2, 1/2]])  # g5


def fourierZ10():
    """
    return a diagonalizing basis for the
    the regular representation of (Z/10Z,*)
    """
    vs = []
    #                0  1   2    3   4   5   6   7   8    9
    vs = np.array([[+4, 0,  0 ,  0,  0,  0,  0,  0,  0,   0],   # v0
                   [-4, 0,  0 ,  0,  0,  4,  0,  0,  0,   0],   # v1
                   [-4, 0,  1 ,  0,  1,  0,  1,  0,  1,   0],   # v2
                   [+0, 0, -1j,  0, -1,  0,  1,  0,  1j,  0],   # v3
                   [+0, 0, -1 ,  0,  1,  0,  1,  0, -1,   0],   # v4
                   [+0, 0,  1j,  0, -1,  0,  1,  0, -1j,  0],   # v5
                   [+4, 1, -1 ,  1, -1, -4, -1,  1, -1,   1],   # v6
                   [+0, 1,  1j, 1j,  1,  0, -1,-1j, -1j, -1],   # v7
                   [+0, 1,  1 , -1, -1,  0, -1, -1,  1,   1],   # v8
                   [+0, 1, -1j,-1j,  1,  0, -1, 1j,  1j, -1]])  # v9
    return np.array(vs/4).transpose()


def generalized_hadamard_transform(a, m):
    """
    m is a k-by-k matrix
    a is a 1d np.array with length a positive power of k
    """
    k = m.shape[0]
    powkn = len(a)  # 16
    powk = 1
    while powk < powkn:  # 1, 2, 4, 8
        for j in range(0, powkn, k*powk):
            for i in range(0, powk):
                a[j+i:j+i+k*powk:powk] = m.dot(a[j+i:j+i+k*powk:powk])
        powk *= k


def hadamard_transform(a):
    """
    Do an in-place 'fast' Walsh-Hadamard transform.
    len(a) must be a power of 2
    """
    pow2n = len(a)
    pow2i = 1
    while pow2i < pow2n:
        for j in range(0, pow2n, 2*pow2i):
            for k1 in range(j, j+pow2i):
                k2 = k1+pow2i
                a[k1], a[k2] = a[k1]+a[k2], a[k1]-a[k2]
        pow2i *= 2


def or_transform(a):
    pow2n = len(a)
    pow2i = 1
    while pow2i < pow2n:
        for j in range(0, pow2n, 2*pow2i):
            for k1 in range(j, j+pow2i):
                k2 = k1+pow2i
                a[k2] = a[k1]+a[k2]
        pow2i *= 2


def inverse_or_transform(a):
    pow2n = len(a)
    pow2i = 1
    while pow2i < pow2n:
        for j in range(0, pow2n, 2*pow2i):
            for k1 in range(j, j+pow2i):
                k2 = k1+pow2i
                a[k2] = -a[k1]+a[k2]
        pow2i *= 2


def convolve(u, v, p):
    n = u.shape[0]
    w = np.zeros((n,), dtype=object)
    for i in range(n):
        for j in range(n):
            w[i] = (w[i] + u[j]*v[(n+i-j)%n]) % p
    return w


def dft_reference(xs, p, inverse=False):
    """
    Just a reference implementation.  Correct but not a 'fft'
    """
    sgn = -1 if inverse else +1

    n = len(xs)
    if (p-1)%n:
        raise ValueError(f'{n} does not divide the order of (Z/{p}Z)*')

    a = primitive_root(p)
    z = pow(a, (p-1)//n, p)

    #
    # z is now an nth root of unity (mod p)
    #
    ys = np.array([0]*n, dtype=object)

    for i in range(n):
        zsgni = pow(z, sgn*i, p)
        coef = 1
        for j in range(n):
            ys[i] += coef*xs[j]
            coef = (coef * zsgni) % p

    return ys % p


def _dft(xs, qs, p, a, inverse=False):
    """
    A safe version of the discrete fourier transform (where
    we already know that the length of n divides p-1)
    """
    n = len(xs)
    sgn = -1 if inverse else +1

    if len(qs)==1:
        return np.array([sum(pow(a, sgn*j*i, p)*xs[i] for i in range(n))
                         for j in range(n)], dtype=object) % p

    q = qs[0]
    z = pow(a, q, p)

    xs2 = np.zeros((n,), dtype=object)
    for i in range(q):
        xs2[i::q] = _dft(xs[i::q], qs[1:], p, z, inverse=inverse)

    xs3 = np.zeros((n,), dtype=object)
    for j in range(n):
        sj = sgn*j
        sj2 = j%(n//q)
        xs3[j] = sum(pow(a, sj*v, p) * xs2[v+sj2*q] for v in range(q)) % p

    return xs3


def dft(xs, p, inverse=False):
    n = len(xs)
    if (p-1)%n:
        raise ValueError(f'Length of xs ({n}) does not divide order of (Z/{p}Z)^*')
    qs = sum(([q]*e for q, e in factorize(n)), [])
    a = pow(primitive_root(p), (p-1)//n, p)

    return _dft(xs, qs, p, a, inverse=inverse)
