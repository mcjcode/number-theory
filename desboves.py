#!/usr/bin/env python

from fractions import Fraction

__all__ = ['desbovesp']

def gcd(a,b) :
    while b != 0 :
        a, b = b, a-(a/b)*b
    return a

def desboves_tangent(x,y,z) :
    """
    Return the third point of intersection of the tangent
    to x^3 + y^3 = a*z^3 at the point (x,y,z).
    """
    
    x,y,z = x*(x**3+2*y**3), -y*(y**3+2*x**3), -z*(y**3-x**3)

    d = gcd(gcd(x,y),z)
    return x/d, y/d, z/d

# C9WM7
# XPB Q3X GM3 R6K T8 PW

# First/Last: Asellus Borealis
# ThirdStarOfGhost@gmail.com
# pswd: angel0cancer
#

def desboves_secant_x(x1,y1,x2,y2) :
    dx = x2-x1
    dy = y2-y1
    return -3*dy**2*(y1*x2-y2*x1)/(dx**3+dy**3) - (x1+x2)

def desboves_secant(x1, y1, x2, y2) :
    x3 = desboves_secant_x(x1, y1, x2, y2)
    y3 = desboves_secant_x(y1, x1, y2, x2)
    return x3, y3

def pequal(p,q) :
    """
    Are p and q equivalent points of
    projective space?
    """
    n = len(p)
    for i in xrange(n-1) :
        for j in xrange(i,n) :
            if p[i]*q[j] != p[j]*q[i] :
                return False
    return True

def desbovesp(p,q) :
    """
    Calculate the product of the two
    points on the elliptic curve
    x^3 + y^3 = a*z^3.

    This routine does not care about
    the value of a, but p and q must
    satisfy the same cubic equation.
    """
    
    # If either of the points is the
    # identity, then return the other
    # point.
    #
    if pequal(p,[1,-1,0]) :
        return q
    
    if pequal(q,[1,-1,0]) :
        return p

    # If one of the points is the other's
    # reflection in the y=x line, they
    # are additive inverses in the group
    # law.  Return the identity element.
    #
    if pequal([p[1],p[0],p[2]],q) :
        return [1,-1,0]
    
    if pequal(p,q) :
        x, y, z = desboves_tangent(*p)
        return y, x, z

    x3, y3 = desboves_secant(Fraction(p[0],p[2]),
                             Fraction(p[1],p[2]),
                             Fraction(q[0],q[2]),
                             Fraction(q[1],q[2]))

    d = gcd(x3.denominator, y3.denominator)

    y, x, z = (x3.numerator*y3.denominator,
               y3.numerator*x3.denominator,
               x3.denominator*y3.denominator)

    # return the answer with no common factors.
    d = gcd(x,gcd(y,z))
    return x/d, y/d, z/d
