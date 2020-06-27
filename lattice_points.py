#!/usr/bin/env python
#
# lattice_points.py
#

from math import floor, ceil
from utilities import sqrtInt

def _gauss_circle1(r):
    """
    Return the number of lattice points in R2
    whose distance to the origin is
    less than or equal to r.

    Reference implementation for checking
    against gauss_circle_norm.
    """
    retval = 0
    rceil = int(ceil(r))
    for x in range(-rceil, rceil+1):
        for y in range(-rceil, rceil+1):
            if x**2 + y**2 <= r**2:
                retval += 1
    return retval

def _gauss_circle2(r):
    return _gauss_circle_norm2(r*r)

def _gauss_circle_norm2(r2):
    """
    Return the number of lattice points in R2
    whose norm is less than or equal to r.

    Reference implementation for checking
    against gauss_circle_norm.
    """
    retval = 0
    for i in range((r2-1)//4+1):
        retval += r2 // (4*i+1) - r2 // (4*i+3)
    retval *= 4
    return retval+1 # don't forget the origin (0, 0)


def gauss_circle(ri, trace=False):
    """
    Return the number of lattice points in R2
    whose distance to the origin is less than
    or equal to r.

    O(r) time complexity algorithm based on
    Jacobi's two square theorem.
    """
    return gauss_circle_norm(r*r, trace)

def gauss_circle_norm(r2, trace=False):
    """
    Return the number of lattice points in R2
    whose norm is less than or equal to r.

    O(sqrt(r2)) time complexity algorithm
    based on Jacobi's two square theorem.
    """
    retval = 0

    quotient = 0
    if trace:
        print('=1(4) divisors')
    for d in range(1, sqrtInt(r2//4)):
        count = (r2//d+3)//4 - (r2//(d+1)+3)//4
        quotient = d
        if trace:
            print(count, quotient, '*')
        retval += count * quotient
        
    ub = (r2//(quotient+1)-1)//4
    prev_quot = quotient

    for i in range(ub, -1, -1):
        quotient = r2 // (4*i+1)
        if quotient == prev_quot:
            continue
        if trace:
            print(1, quotient)
        retval += quotient 

    if trace:
        print('=3(4) divisors')     
    quotient = 0   
    for d in range(1, sqrtInt(r2//4)):
        count = (r2//d+1)//4 - (r2//(d+1)+1)//4
        quotient = d
        if trace: 
            print(count, quotient, '*')
        retval -= count * quotient
    ub = (r2//(quotient+1)-1)//4
    prev_quot = quotient

    for i in range(ub, -1, -1):
        quotient = r2 // (4*i+3)
        if quotient == prev_quot:
            continue
        if trace:
            print(1, quotient)
        retval -= quotient

    retval *= 4

    return retval + 1 # don't forget the origin        
