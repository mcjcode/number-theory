import math

from utilities import factorize2, sqrtInt, cbrtInt, nrtInt
from prime_sieve import segmented_sieve
from lucy import sievecnt


def biprimes(N):
    retval = 0
    V, S0 = sievecnt(N)
    ps = segmented_sieve(sqrtInt(N))
    for i, p in enumerate(ps):
        retval += S0[N//p] - (i+1) + 1
    return retval


def triprimes(N):
    retval = 0
    V, S0 = sievecnt(N)
    ps = list(segmented_sieve(sqrtInt(N)))
    i = 0
    while ps[i]**3 <= N:
        p1 = ps[i]
        j = i
        while j<len(ps) and ps[j]**2 <= N//p1:
            p2 = ps[j]
            x = S0[N//(p1*p2)] - (j+1) + 1
            retval += x
            j += 1
        i += 1
    return retval


def quadprimes(N):
    retval = 0
    V, S0 = sievecnt(N)
    ps = list(segmented_sieve(sqrtInt(N)))
    i = 0
    while ps[i]**4 <= N:
        p1 = ps[i]
        j = i
        while j<len(ps) and ps[j]**3 <= N//p1:
            p2 = ps[j]
            k = j
            while k<len(ps) and ps[k]**2 <= N//(p1*p2):
                p3 = ps[k]
                x = S0[N//(p1*p2*p3)] - (k+1) + 1
                retval += x
                k += 1
            j += 1
        i += 1
    return retval


def almost_primes(k, N, i=0, j=1, ps=None, V=None, S0=None):
    if j==1:
        ps = list(segmented_sieve(sqrtInt(N)))
        V, S0 = sievecnt(N)
        
    if k==1:
        retval = S0[N] - (i+1) + 1
    else:
        retval = 0
        while i<len(ps) and ps[i]**k <= N:
            p = ps[i]
            retval += almost_primes(k-1, N//p, i, j+1, ps, V, S0)
            i += 1
    return retval
    


def pi2(a, b, N, ps, V, S0):
    g = math.gcd(a, b)
    if g>1:
        return pi2(a//g, b//g, nrtInt(g, N), ps, V, S0)

    if a>b:
        a, b = b, a

    if b<2:
        raise ValueError(f'q exponent {b=} must be at least 2')
        
    retval = 0
    for q in ps:
        if q**b > N:
            break
        if q**(a+b) <= N:
            retval += S0[nrtInt(a, N//q**b)] - 1
        else:
            retval += S0[nrtInt(a, N//q**b)]
    return retval

def pqr(N, ps, V, S0):
    #
    # p*q*r
    #
    retval = 0
    pi = -1
    while (pi:=pi+1) < len(ps) and (p:=ps[pi])<=cbrtInt(N):
        qi = pi
        while (qi:=qi+1) < len(ps) and (q:=ps[qi])<=sqrtInt(N//p):
            retval += S0[N//(p*q)] - S0[q]
    return retval

def p2qr(N, ps, V, S0):
    #
    # p**2*q*r
    #
    retval = 0
    pi = -1
    while (pi:=pi+1) < len(ps) and (p:=ps[pi])<=sqrtInt(N):
        qi = -1
        while (qi:=qi+1) < len(ps) and (q:=ps[qi])<=sqrtInt(N//p**2):
            if q==p:
                continue
            retval += (S0[N//(p**2*q)] - (p**3*q<=N)) - (S0[q] - (p<q))
    return retval

def p2q2r(N, ps, V, S0):
    #
    # p2*q2*r
    #
    retval = 0
    pi = -1
    while (pi:=pi+1) < len(ps) and (p:=ps[pi])<=nrtInt(4, N):
        qi = pi
        while (qi:=qi+1) < len(ps) and (q:=ps[qi])<=sqrtInt(N//p**2):
            retval += S0[N//(p**2*q**2)] - (p<=N/(p**2*q**2)) - (q<=N/(p**2*q**2))
    return retval

def pkqr(N, k, ps, V, S0):
    #
    # p**k*q*r
    #
    retval = 0
    pi = -1
    while (pi:=pi+1) < len(ps) and (p:=ps[pi])<=nrtInt(k, N):
        qi = -1
        while (qi:=qi+1) < len(ps) and (q:=ps[qi])<=sqrtInt(N//p**k):
            if q==p:
                continue
            bd = N//(p**k*q)
            if bd<=q:
                rcount = 0
            else:
                rcount = (S0[bd] - (p<=bd)) - (S0[q] - (p<=q))
            retval += rcount
    return retval


def pqrs(N, ps, V, S0):
    retval = 0
    pi = -1
    while (pi:=pi+1)<len(ps) and (p:=ps[pi])<=nrtInt(4, N):
        qi = pi
        while (qi:=qi+1)<len(ps) and (q:=ps[qi])<=nrtInt(3, N//p):
            ri = qi
            while (ri:=ri+1)<len(ps) and (r:=ps[ri])<=nrtInt(2, N//(p*q)):
                retval += S0[N//(p*q*r)] - S0[r]
    return retval
    
def p2qrs(N, ps, V, S0, log=lambda *x: None):
    x = 0
    pi = -1
    while (pi:=pi+1)<len(ps) and (p:=ps[pi])<sqrtInt(N):
        log(f'{p=}')
        qi = -1
        while (qi:=qi+1)<len(ps) and (q:=ps[qi])**3<(N//p**2):
            if q!=p:
                log(f'  {q=}')
                ri = qi
                while (ri:=ri+1)<len(ps) and (r:=ps[ri])**2<(N//(p**2*q)):
                    if r!=p:
                        bd = N//(p**2*q*r)
                        if bd > r:
                            scount = (S0[bd] - (p<=bd) - (q<=bd) - (r<=bd)) - (S0[r] - (p<=r) - (q<=r) - (r<=r))
                        else:
                            scount = 0
                        log(f'    {r=} {scount=}')
                        x += scount
    return x
