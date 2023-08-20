#!/usr/bin/env python -O
#
# lucy.py
#

from utilities import sqrtInt, cbrtInt

fwidth = 5
column_limit = 20

def printV(V,tag=''):
    fmtstr = '%%%dd' % (fwidth,)
    print(' '.join([fmtstr % (i,) for i in V[-column_limit:]])+(' :%s'%(tag,)))

def printS(V,S,tag=''):
    fmtstr = '%%%dd' % (fwidth,)
    print(' '.join([fmtstr % (S[i],) for i in V[-column_limit:]])+(' :%s'%(tag,)))

def printRuler(V,tag=''):
    print('-'.join(['-'*fwidth for i in V[-column_limit:]])+(' :%s'%(tag,)))

def sievecnt(n):
    """
    Return the count of numbers not excluded
    by a sieve of [1..n] by all primes up to
    and including p.  If p is not given then
    sieve up to sqrt(n). 
    """
    sqrtn = sqrtInt(n)
    V = [n//i for i in range(1,sqrtn+1)]
    V += list(range(V[-1]-1,0,-1))
    S = {i:i-1 for i in V}
    for p in range(2,sqrtn+1):
        if S[p] > S[p-1]:  # p is prime
            sp = S[p-1]  # number of primes smaller than p
            p2 = p*p
            for v in V:
                if v < p2: break
                S[v] -= (S[v//p] - sp)
    return V, S    


def sievecnt_mod4(n):
    """
    Return the number of primes congruent to 1 and 3 (mod 4)
    for less than or equal to n, n//2, n//3, ..., 1
    """
    sqrtn = sqrtInt(n)
    V = [n//i for i in range(1,sqrtn+1)]
    V += list(range(V[-1]-1,0,-1))
    S_1mod4 = {v: (v-1)//4 for v in V} # counting primes = 1(4)
    S_3mod4 = {v: (v+1)//4 for v in V} # counting primes = 3(4)
    for p in range(3,sqrtn+1):
        p2 = p*p
        if S_1mod4[p] > S_1mod4[p-1]: # p is a prime of the form 4k+1
            s1m4pm1 = S_1mod4[p-1]
            s3m4pm1 = S_3mod4[p-1]
            for v in V:
                if v < p2: break
                S_1mod4[v] -= (S_1mod4[v//p]-s1m4pm1)
                S_3mod4[v] -= (S_3mod4[v//p]-s3m4pm1)
        elif S_3mod4[p] > S_3mod4[p-1]: # p is a prime of the form 4k+3
            s1m4pm1 = S_1mod4[p-1]
            s3m4pm1 = S_3mod4[p-1]
            for v in V:
                if v < p2: break 
                S_1mod4[v] -= (S_3mod4[v//p]-s3m4pm1)
                S_3mod4[v] -= (S_1mod4[v//p]-s1m4pm1)

    return V, S_1mod4, S_3mod4


def sievesum(n):
    """
    Return the sum of all numbers not excluded
    by a sieve of [1..n] by all primes up to
    and including p.
    """
    sqrtn = sqrtInt(n)
    V = [n//i for i in range(1,sqrtn+1)]
    V += list(range(V[-1]-1,0,-1))
    S = {i:i*(i+1)//2-1 for i in V}
    for p in range(2,sqrtn+1):
        if S[p] > S[p-1]:  # p is prime
            sp = S[p-1]  # sum of primes smaller than p
            p2 = p*p
            for v in V:
                if v < p2: break
                S[v] -= p*(S[v//p] - sp)
    return V, S

def sievecntsum(n):
    """
    Return both the counts and the sums of the
    primes up to n//k for all positive k.
    """
    sqrtn = sqrtInt(n)
    V = [n//i for i in range(1,sqrtn+1)]
    V += list(range(V[-1]-1,0,-1))
    S0 = {i:i-1          for i in V}
    S1 = {i:i*(i+1)//2-1 for i in V}
    for p in range(2,sqrtn+1):
        if S0[p] > S0[p-1]: # p is prime
            p2 = p*p
            for v in V:
                if v < p2: break
                vmodp = v//p
                D0 =  S0[vmodp] - S0[p-1]
                D1 =  S1[vmodp] - S1[p-1]
                S0[v] -=    D0
                S1[v] -= p*(D1)
    return V, S0, S1

def P23(S0,V,printTables=False):
    """
    V: all positive numbers of the form n//k
       in descending order: [n,n//2,n//3,...,3,2,1]
    S0: gives pi(v) for all v in V
    
    returns: a dictionary S with S[v] = number
    of p*q <= v where p<=q are primes in (cbrtn,n]
    """

    Vr = list(reversed(V))
    if printTables:
        printV(Vr); printS(Vr,S0); printRuler(Vr)

    n = V[0]
    cbrtn = cbrtInt(n)
    cbrtn2 = cbrtn**2
    sqrtn = sqrtInt(n)
    a = S0[cbrtn]

    S = {v:0 for v in V}
    for v in V:
        if v < cbrtn2: break
        cbrtv = cbrtInt(v)
        sqrtv = sqrtInt(v)
        pb = cbrtn+1
        #
        # Effectively, we loop over every prime
        # pb in the range [cbrtn,sqrtv]
        #
        while pb <= sqrtv: ##ub: #sqrtn:
        #for pb in range(cbrtn+1,sqrtn+1): ## cbrtv+1,sqrtv+1):
            #if pb > v//pb: break
            if S0[pb] > S0[pb-1]: # pb is prime
                #print('working on prime=%d'%(pb,))
                #print('%d primes less than %d//%d=%d' % (S0[v//pb],v,pb,v//pb))
                #print('primes in [%d,%d]=%d' % (pb,v//pb,S0[v//pb]-(S0[pb]-1)))
                S[v] += S0[v//pb]#-(S0[pb]-1)
            pb += 1
        S[v] -= S0[sqrtv]*(S0[sqrtv]-1)//2 - a*(a-1)//2

    return S
