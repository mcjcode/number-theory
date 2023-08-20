#
# slv.py - Shortest Lattice Vector routines
#

import math

def lagrange(u, v):
    """
    An implementation of the two-dimensional Gauss-Lagrange
    algorithm for computing a reduced basis w.r.t. the L2 norm.
    
    :param u: a vector in Rn
    :param v: a vector in Rn
    :return: a reduced basis for the lattice spanned by u and v.
 
    A reduced basis u, v is a basis for which |u| <= |v| <= |v-ku|
    for all integers k.  For a reduced basis u, v, |u| is the minimum
    of the norm on the lattice, and |v| either equals |u| or is the
    next largest value taken on by the norm on the lattice.
    """
    norm = lambda v: v.dot(v)
    if norm(u)==0:
        return u, v
    if norm(u) > norm(v):
        u, v = v, u
    mu = int(round(u.dot(v)/norm(u)))
    while mu:
        v -= mu*u
        if norm(u) > norm(v):
            u, v = v, u
        if norm(u)==0:
            return u, v
        mu = int(round(u.dot(v)/norm(u)))
    return u, v

def mint(u, v):
    """
    Find an integer t that minimizes ||v-tu||_1.
    """
    f = lambda t: sum(abs(v-t*u))

    # if b==0, then the value of a-tb doesn't vary with t
    z = sorted([((a/b), a, b) for a, b in zip(v, u) if b])
    absbs = [abs(b) for (_, _, b) in z]
    slopes = [-sum(absbs)]
    for absb in absbs:
        slopes.append(slopes[-1]+2*absb)

    m = math.inf 
    argm = None
    for i in range(len(slopes)-1):
        if slopes[i]<=0 and slopes[i+1]>=0:
            # we have a local minimum at t=a/b
            t, a, b = z[i]
            q, r = divmod(a, b)
            ts = [q] if r==0 else [math.floor(t), math.ceil(t)]
            for t0 in ts:
                x = f(t0)
                if x<m:
                    m = x
                    argm = t0
    return argm


def lagrangeL1(u, v):
    """
    An implementation of the two-dimensional Gauss-Lagrange
    algorithm for computing a reduced basis w.r.t. the L1 norm.
    
    :param u: a vector in Rn
    :param v: a vector in Rn
    :return: a reduced basis for the lattice spanned by u and v.
 
    A reduced basis u, v is a basis for which |u| <= |v| <= |v-ku|
    for all integers k.  For a reduced basis u, v, |u| is the minimum
    of the norm on the lattice, and |v| either equals |u| or is the
    next largest value taken on by the norm on the lattice.
    """
    norm = lambda w: sum(map(abs,w))
    if norm(u) > norm(v):
        u, v = v, u
    mu = mint(u, v)
    while mu:
        v -= mu*u
        if norm(u) > norm(v):
            u, v = v, u
        mu = mint(u, v)
    return u, v

