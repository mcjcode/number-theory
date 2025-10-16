# -*- coding: utf-8 -*-


import numpy as np
import itertools
from utilities import factorize, factorize2
from numpy import poly1d


def total_poly(p, n):
    """
    Return the polynomial p(x) that is satisfied
    by all elements of the finite field of degree
    p^n.  (i.e. x^(p^n)-x)
    """
    return poly1d([1, 0])**(p**n) - poly1d([1, 0])


def product_of_irreducibles(d, p):
    """
    Return the product of all of the degree d
    irreducible polynomials.
    """
    prime_factors = list(set(factorize(d)))

    numer = total_poly(p, d)
    denom = poly1d([1])

    for i in range(1, len(prime_factors)+1):
        for subset in itertools.combinations(prime_factors, i):
            fact = total_poly(p, d/np.prod(subset))
            if i % 2 == 0:
                numer *= fact
            else:
                denom *= fact
    return (numer/denom)[0]


class FiniteFieldElement(object):

    def __init__(self, value, f):
        self.value = value
        self.f = f

    def __repr__(self):
        return 'FiniteFieldElement(%s, %s)' % (self.value, self.f)

    def __add__(self, other):
        return FiniteFieldElement((self.value+other.value) % self.f.p, self.f)

    def __iadd__(self, other):
        self.value = (self.value + other.value) % self.f.p
        return self

    def __sub__(self, other):
        return FiniteFieldElement((self.value-other.value) % self.f.p, self.f)

    def __neg__(self):
        return FiniteFieldElement(-self.value % self.f.p, self.f)

    def __pow__(self, xp):
        retval = self.f.one()
        runpow = self
        while xp:
            if xp & 1:
                retval *= runpow
            runpow = runpow * runpow
            xp >>= 1
        return retval

    def __mul__(self, other):
        return FiniteFieldElement(self.f.mul(self.value, other.value), self.f)

    def __imul__(self, other):
        self.value = self.f.mul(self.value, other.value)
        return self

    def __eq__(self, other):
        return (self.f == other.f) and (self.value == other.value).all()

    def __ne__(self, other):
        return not (self == other)


class FiniteField(object):

    def __init__(self, p, n):
        self.p = p
        self.n = n
        self.size = self.p**self.n
        self.rpoly = self.defining_polynomial()

    def polynomial_remainder_modp(self, dvd, dvs):
        pc = poly1d(dvd.c)
        x = poly1d([1, 0])
        while pc.order >= dvs.order:
            pc = pc - dvs*(x**(pc.order-dvs.order))*pc.c[0]
            pc = poly1d([coef % self.p for coef in pc.c])
        return pc

    def defining_polynomial(self):
        indicators = np.zeros((self.p, )*self.n, dtype=int)
        for i in range(1, self.n//2+1):
            for acoefs in itertools.product(*[range(self.p)]*i):
                for bcoefs in itertools.product(*[range(self.p)]*(self.n-i)):
                    ccoefs = (poly1d((1, )+acoefs)*poly1d((1, )+bcoefs)).c
                    ccoefs = tuple(x % self.p for x in ccoefs)
                    # print(ccoefs)
                    indicators[tuple(ccoefs[1:])] = 1
        # print(indicators)
        for k in itertools.product(*[range(self.p)]*self.n):
            if indicators[k] == 0:
                return poly1d((1, )+k)

    def defining_polynomial_old(self):
        #
        # find the product of all the monic irreducible
        # polynomials of degree d over Z/pZ.
        #
        pr = product_of_irreducibles(d=self.n, p=self.p)
        #
        # now loop over all degree d monic polynomials
        # and find a factor.  This will be guaranteed
        # to be irreducible.
        #
        for coefs in itertools.product(*[range(self.p)]*(self.n-1)):
            coefs = [1] + list(coefs) + [1]  # monic polynomial of degree self.n
            trial_poly = poly1d(coefs)
            r = self.polynomial_remainder_modp(pr, trial_poly)
            if all(r.c == 0):
                return trial_poly

    def __eq__(self, other):

        return self.p == other.p and self.n == other.n and self.rpoly == other.rpoly

    def order2(self, a):
        """
        :param a: an element of the finite field
        :return: the multiplicative order of the element in the field

        A brute force reference implementation.
        """
        cnt = 1
        b = a
        one = self.one()
        minus_one = -one
        while not b == one:
            if b == minus_one:
                return 2*cnt
            b = b * a
            cnt += 1
        return cnt

    def order(self, a):
        """
        :param a: an element of the finite field
        :return: the multiplicative order of the element in the finite field

        Seems about 2x as fast as a brute force attempt.
        """
        e = self.p**self.n-1
        one = self.one()
        for (q, f) in factorize2(e):
            e //= q**f
            a1 = a**e
            while a1 != one:
                a1 = a1**q
                e *= q
        return e
    
    def primitive_element(self):
        for elt in self:
            if elt == self.zero():
                continue
            if self.order(elt) == self.size-1:
                return elt

    def __repr__(self):
        return 'FiniteField(%d, %d)' % (self.p, self.n)

    def __str__(self):
        return self.__repr__()

    def mul(self, x, y):
        z = poly1d(list(reversed(x)))*poly1d(list(reversed(y)))
        for i in range(z.order+1):
            z[i] %= self.p
        z = self.polynomial_remainder_modp(z, self.rpoly)
        return np.array(list(map(int,list(reversed(z.c))+[0]*(self.n-len(z.c)))))

    def zero(self):
        return FiniteFieldElement(np.array([0]*self.n, dtype=int), self)

    def one(self):
        return FiniteFieldElement(np.array([1]+[0]*(self.n-1), dtype=int), self)

    def __iter__(self):
        return (FiniteFieldElement(np.array(val, dtype=int), self)
                for val in itertools.product(*[range(self.p)]*self.n))

    def jacobi_sum(self, order):
        """
        Returns the Jacobi sum J(chi, chi) associated with the
        character which maps the smallest primitive
        root to the ord-root of unity
        """
        pn = self.p**self.n
        a = self.primitive_element()
        theta = 2 * np.pi / order
        zeta = np.cos(theta)+1j*np.sin(theta)

        chars = np.zeros((self.p, )*self.n, dtype=type(0j))
        chars[(0, )*self.n] = 0j
        pa = self.one()
        char_val = 1+0j
        for k in range(1, pn):
            chars[tuple(pa.value)] = char_val
            pa = pa * a
            char_val *= zeta

        charsum = 0j
        for aindex in itertools.product(*[range(self.p)]*self.n):
            bindex = list((self.p-ai) % self.p for ai in aindex)
            bindex[0] = (bindex[0]+1) % self.p
            bindex = tuple(bindex)
            charsum += chars[aindex]*chars[bindex]

        c2 = charsum.imag/zeta.imag
        tmp = charsum - c2*zeta
        c1 = tmp.real
        return round(c1), round(c2)


def count_projective_plane_curve_points(f, ffield):
    accum = 0
    zero = ffield.zero()
    one = ffield.one()
    for x in ffield:
        for y in ffield:
            accum += f(x, y, one) == zero
        accum += f(x, one, zero) == zero
    accum += f(one, zero, zero) == zero
    return accum


def projective_plane_curve_points(f, ffield):
    pts = []
    zero = ffield.zero()
    one = ffield.one()
    for x in ffield:
        for y in ffield:
            if f(x, y, one) == zero:
                pts.append((x, y, one))
        if f(x, one, zero) == zero:
            pts.append((x, one, zero))
    if f(one, zero, zero) == zero:
        pts.append((one, zero, zero))
    return pts


def count_projective_curve_points(f, ffield):
    accum = 0
    zero = ffield.zero()
    for x in ffield:
        for y in ffield:
            for z in ffield:
                if x != zero or y != zero or z != zero:
                    accum += (f(x, y, z) == zero)
    accum //= ffield.p**ffield.n-1
    return accum


def count_curve_points_affine(f, ffield):
    accum = 0
    zero = ffield.zero()
    for x in ffield:
        for y in ffield:
            if f(x, y) == zero:
                accum += 1
    return accum


def count_diagonal_cubic_points(ffield):
    accum = 0
    zero = ffield.zero()
    cubes = []
    for x in ffield:
        if x**3 not in cubes:
            cubes.append(x**3)
    # cubes = list(set(x**3 for x in ffield))
    if not len(cubes) == (len(list(ffield))-1)/3 + 1:
        print(cubes)
    for x in cubes:
        cx = (1 if x == zero else 3)
        for y in cubes:
            cy = (1 if y == zero else 3)
            for z in cubes:
                cz = (1 if z == zero else 3)
                if x+y+z == zero:
                    accum += cx*cy*cz
    accum -= 1  # (0, 0, 0) is not a projective solution
    assert accum % (ffield.p**ffield.n-1) == 0
    accum //= ffield.p**ffield.n-1
    return accum
