#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

from math import sqrt, floor
import numpy as np
from base_complex import infj
from utilities import gcd, bezout


def _rectangle_n_points(n):
    for i in range(2*n):
        yield -n+i, n
        yield n, n-i
        yield n-i, -n
        yield -n, -n+i


def plot_rectangle_n_points():
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(5, 5))
    ax = fig.gca()
    for x, y in _rectangle_n_points(5):
        plt.plot(x, y, '.')
    ax.set_xlim([-7, +7])
    ax.set_ylim([-7, +7])

    plt.show()


def _eisenstein_bound(eps, k, z):
    """
    Return the number of rectangles of lattice points to compute phi_0, k(z) within eps.
    """

    x = z.real
    y = z.imag

    h = min(y, y/sqrt(x**2+y**2))

    nmax = int(floor((8/(h**(2*k)*(2*k-2)*eps))**(1./(2*k-2))))

    return nmax


def unrestricted_eisenstein(k, z):
    """
    Return the value G_k(z) of the unrestricted Eisenstein series of weight k at z.

    G_k(z) = Sum_{(c, d), c>=0, d>0 if c=0} (1/(cz+d)^(2k)).
    """
    if z.imag <= 0:
        raise ValueError("z not in the upper half plane.")

    nmax = _eisenstein_bound(1e-08, k, z)

    retval = 0.0
    for n in range(nmax, 0, -1):
        for (c, d) in _rectangle_n_points(n):
            if c >= 0 and (c > 0 or d > 0):
                v = 1./(c*z+d)**(2*k)
                retval += v
    return retval


def eisenstein(k, z):
    """
    Return the value phi_0, k(z) of the Eisenstein series of weight k at z.
    """

    if z.imag <= 0:
        raise ValueError("z not in the upper half plane.")

    nmax = _eisenstein_bound(1e-08, k, z)

    retval = 0.0
    for n in range(nmax, 0, -1):
        for (c, d) in _rectangle_n_points(n):
            if c >= 0 and (c > 0 or d > 0) and gcd(c, d) == 1:
                v = 1/(c*z+d)**(2*k)
                retval += v
    return retval


def poincare(k, nu, z):
    """
    Return the value of the Poincare series of weight k and character nu at z.
    """

    if z.imag <= 0:
        raise ValueError("z not in the upper half plane.")

    nmax = _eisenstein_bound(1e-08, k, z)

    retval = 0.0
    for n in range(nmax, 0, -1):
        for (c, d) in _rectangle_n_points(n):
            if c >= 0 and (c > 0 or d > 0) and gcd(c, d) == 1:
                b, a = bezout(c, d)
                Tz = (a*z-b)/(c*z+d)
                v = np.exp(2.0*np.pi*nu*1.0j*Tz)/(c*z+d)**(2*k)
                retval += v
    return retval


rho = np.exp(2.0*np.pi*1.0j/6.0)


def _coset_reps(n):
    for a in range(1, n+1):
        if n % a == 0:
            d = n // a
            for b in range(d):
                yield a, b, c, d


def transf(a, b, c, d):
    def f(z):
        if z == infj:
            if c == 0:
                return infj
            else:
                return (1.0*a)/c
        if c*z+d == 0:
            return infj
        else:
            return (a*z+b)/(c*z+d)
    return f


def hecke_operator(n, k, f):
    """
    Return the application of the nth hecke operator to a function.
    """
    def f(z):
        return sum((c*z+d)**(2.*k)*transf(a, b, c, d)(z) for (a, b, c, d) in _coset_reps(n))
    return f
