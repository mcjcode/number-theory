# -*- coding: utf-8 -*-

import numpy as np

from utilities import (
    primitive_root,
    isprime,
    modpow,
    gcd,
    squarefree,
    )
from quadratic_extensions import legendre_ch

from real_quadratic_fields import (
    fundamental_unit,
    discriminant,
    )


def gs_numerical(ch, m, k):
    r"""
    :param ch: a mod m character
    :param m: the modulus
    :param k: an integer
    :return: the Gauss sum :math:`\displaystyle\sum_a(ch(a)\omega^{ak})`

    :math:`\omega = e^{2\pi/m}`

    """
    w = np.exp(2.0j*np.pi/m)
    return sum(ch(a)*w**(k*a) for a in range(1, m))


def L_one_chi(ch, m):
    r"""
    Compute the value of the Dirichlet L-series L(s, ch) at s=1.
    """
    z = np.exp(2.0j*np.pi/m)
    return -(1.0/m) * sum(gs_numerical(ch, m, k)*np.log(1-z**(-k)) for k in range(1, m) if gcd(k, m) == 1)


def all_modm_chars(m):
    r"""
    Returns a list of all mod m characters,
    starting with the trivial character.  Only
    works with prime m for now.
    """
    if not isprime(m):
        raise ValueError("%d is not a prime." % (m, ))
    a = primitive_root(m)
    z = np.exp(2.0j*np.pi/(m-1))

    vals = [None]*(m-1)
    pows = [0] + [pow(a, k, m) for k in range(m-1)]
    for k in range(m-1):
        pows2 = [0] + [(z**k)**i for i in range(m-1)]
        vals[k] = dict(zip(pows, pows2))
    return map(lambda val: (lambda i: val[i%m]), vals)

# Calculate the value of kappa for the Q[z_5], the cyclotomic
# field of 5th roots of unity.  This has 4 complex embeddings,
# no real embeddings, (s=2, r=0, n=r+s=2), (1+z_5) is a fundamental
# unit, the regulator = (1/2)log(||(1+z_5)||), the discriminant
# of the Q[z_5] is 5^(5-2) (i.e. disc(Q[z_p])=p^(p-2)), and there
# are 5 roots of unity in the field.
#
# >>> (2.0*np.pi)**2.0*(0.5)*np.log((1+z)*(1+z).conj())/(5.0*np.sqrt(abs(5**3)))
# (0.33983727824052351+0j)
#
# now compare to the product of the L(1, ch) as ch runs over
# all non-trivial mod 5 characters:
#
# >>> np.product([L_one_chi(ch, 5) for ch in all_modm_chars(5)[1:]])
# (0.33983727824052373+0j)


def kappa(m):
    r"""
    :param m: a square-free integer
    :return: kappa = lim_{t->infty} #{J | ||J|| < t}/t for Q[√m].

    We have that
    :math: abs(L(chi, 1)) = h*kappa
    where h is the class number of Q[√m],
    chi is the 'legendre character' and L is the associated
    Dirichlet L function.
    """
    if not squarefree(m):
        raise ValueError('%d is not square free.' % (m, ))

    if m == 1:
        return 1

    disc = discriminant(m)
    if m > 0:
        u = fundamental_unit(m).real()
        return 2*np.log(u)/np.sqrt(abs(disc))
    else:
        if m == -1:
            w = 2
        elif m == -3:
            w = 3
        else:
            w = 1
        return np.pi / (w*np.sqrt(abs(disc)))


def ideal_class_number(d):
    r"""
    Return the ideal class number of Q[√d], for squarefree integer d.
    """
    if not squarefree(d):
        raise ValueError("%d is not squarefree." % (d, ))

    if d == 1:
        return 1

    abs_disc = abs(discriminant(d))

    ch = legendre_ch(d)

    z = np.exp(2.0j*np.pi/abs_disc)
    rho = np.abs(1.0/np.sqrt(abs(abs_disc)) * sum(
        ch(k) * np.log(1-z**(-k)) for k in range(1, abs_disc) if gcd(k, abs_disc) == 1)
        )

    k = kappa(d)

    # magically rho should be an integral multiple of kappa.
    # make sure this is true before we return the rounded answer.
    assert(abs(round(rho/k) - rho/k) < 0.001)

    return int(round(rho/k))
