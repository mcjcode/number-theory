from fractions import Fraction
from utilities import sqrtInt

def continued_fraction(n,d):
    """
    yield the terms of the continued fraction
    representation of n/d
    """
    while n%d :
        q, r = divmod(n,d)
        yield n//d
        n, d = d, r
    yield n

def rational_from_cf(cf):
    """
    Return the rational number (a Fraction)
    corresponding to a continued fraction.
    """
    if len(cf)==0:
        return Fraction(1,1)
    elif len(cf)==1:
        return Fraction(cf[0],1)
    else:
        return Fraction(cf[0],1) + 1 / rational_from_cf(cf[1:])

def sum_sq_rep_benjamin(p):
    """
    A prime p=2 or congruent to 1 (mod 4) has a 
    representation as a sum of 2 squares.
    Return such a representation.

    Here we use the algorithm of Benjamin and Zeilberger (2005)
    in their paper 'Pythagorean Primes and Palindromic Continued
    Practions'
    """
    if p==2:
        return 1,1

    for m in range(sqrtInt(p),p//2+1):
        cf = list(continued_fraction(p,m))
        ln = len(cf)
        if cf==list(reversed(cf)):
            if ln%2==0:
                a = rational_from_cf(cf[:ln//2]).numerator
                b = rational_from_cf(cf[:ln//2-1]).numerator
                return a,b
