#!/usr/bin/env python -O
#
# lucy.py
#

import math
from utilities import sqrtInt, cbrtInt, modinv

fwidth = 5
column_limit = 20
debug = False

def printV(V,tag=''):
    if debug:
        fmtstr = '%%%dd' % (fwidth,)
        print(' '.join([fmtstr % (i,) for i in V[-column_limit:]])+(' :%s'%(tag,)))

def printS(V,S,tag=''):
    if debug:
        fmtstr = '%%%dd' % (fwidth,)
        print(' '.join([fmtstr % (S[i],) for i in V[-column_limit:]])+(' :%s'%(tag,)))

def printRuler(V,tag=''):
    if debug:
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
    printV(V)
    printRuler(V)
    printS(V,S)
    for p in range(2,sqrtn+1):
        if S[p] > S[p-1]:  # p is prime
            sp = S[p-1]  # number of primes smaller than p
            p2 = p*p
            for v in V:
                if v < p2: break
                S[v] -= (S[v//p] - sp)
            printS(V, S, tag=p)
    return V, S    


def sievecnt_mod4(n):
    """
    Return the number of primes congruent to 1 and 3 (mod 4)
    for less than or equal to n, n//2, n//3, ..., 1
    """
    sqrtn = sqrtInt(n)
    V = [n//i for i in range(1,sqrtn+1)]
    V += list(range(V[-1]-1,0,-1))
    Vr = list(reversed(V))
    
    S_1mod4 = {v: (v-1)//4 for v in V} # counting primes = 1(4)
    S_3mod4 = {v: (v+1)//4 for v in V} # counting primes = 3(4)
    printV(Vr)
    printRuler(Vr)
    printS(Vr, S_1mod4)
    printS(Vr, S_3mod4)
    printRuler(V)        
    for p in range(3,sqrtn+1):
        p2 = p*p
        if S_1mod4[p] > S_1mod4[p-1]: # p is a prime of the form 4k+1
            s1m4pm1 = S_1mod4[p-1]
            s3m4pm1 = S_3mod4[p-1]
            for v in V:
                if v < p2: break
                S_1mod4[v] -= (S_1mod4[v//p]-s1m4pm1)
                S_3mod4[v] -= (S_3mod4[v//p]-s3m4pm1)
            printS(Vr, S_1mod4, tag='%s, %s' % (p, 1))
            printS(Vr, S_3mod4, tag='%s, %s' % (p, 3))
            printRuler(V)            
        elif S_3mod4[p] > S_3mod4[p-1]: # p is a prime of the form 4k+3
            s1m4pm1 = S_1mod4[p-1]
            s3m4pm1 = S_3mod4[p-1]
            for v in V:
                if v < p2: break 
                S_1mod4[v] -= (S_3mod4[v//p]-s3m4pm1)
                S_3mod4[v] -= (S_1mod4[v//p]-s1m4pm1)
            printS(Vr, S_1mod4, tag='%s, %s' % (p, 1))
            printS(Vr, S_3mod4, tag='%s, %s' % (p, 3))
            printRuler(V)
    return V, S_1mod4, S_3mod4


def sievecnt_mod3(n):
    """
    Return the number of primes congruent to 1 and 2 (mod 3)
    for less than or equal to n, n//2, n//3, ..., 1
    """
    sqrtn = sqrtInt(n)
    V = [n//i for i in range(1,sqrtn+1)]
    V += list(range(V[-1]-1,0,-1))
    Vr = list(reversed(V))
    S_1mod3 = {v: (v-1)//3 for v in V} # counting primes = 1(3)
    S_2mod3 = {v: (v+1)//3 for v in V} # counting primes = 2(3)
    printV(Vr)
    printRuler(Vr)
    printS(Vr, S_1mod3)
    printS(Vr, S_2mod3)
    printRuler(V)    
    for p in range(2,sqrtn+1):
        p2 = p*p
        if S_1mod3[p] > S_1mod3[p-1]: # p is a prime of the form 3k+1
            s1m3pm1 = S_1mod3[p-1]
            s2m3pm1 = S_2mod3[p-1]
            for v in V:
                if v < p2: break
                S_1mod3[v] -= (S_1mod3[v//p]-s1m3pm1)
                S_2mod3[v] -= (S_2mod3[v//p]-s2m3pm1)
            printS(Vr, S_1mod3, tag='%s, %s' % (p, 1))
            printS(Vr, S_2mod3, tag='%s, %s' % (p, 2))
            printRuler(Vr)                                                
        elif S_2mod3[p] > S_2mod3[p-1]: # p is a prime of the form 3k+2
            s1m3pm1 = S_1mod3[p-1]
            s2m3pm1 = S_2mod3[p-1]
            for v in V:
                if v < p2: break 
                S_1mod3[v] -= (S_2mod3[v//p]-s2m3pm1)
                S_2mod3[v] -= (S_1mod3[v//p]-s1m3pm1)
            printS(Vr, S_1mod3, tag='%s, %s' % (p, 1))
            printS(Vr, S_2mod3, tag='%s, %s' % (p, 2))
            printRuler(Vr)                                
    return V, S_1mod3, S_2mod3


def _f13(n):
    """
    return the sum of all positive numbers <= n
    congruent to 1(mod3)
    """
    #
    #   1  2  3  4  5  6  7  8  9 10 11 12 13
    #   0  0  0  4  4  4 11 11 11 21 21 21 34
    #
    # of the form 3n-2, 3n-1 or 3n then the sum is
    #
    #     1 +   4 +   7 + ... + 3n-2 
    #  3n-2  3n-5 +     +   4 +    1
    #
    #  (3n-1)*n//2 - 1
    #
    if (n+1)%3==0: #n is of the form 3n-1
        n -= 1
    elif n%3==0:
        n -= 2
    return (n+1)*((n+2)//3)//2


def _f23(n):
    """
    return the sum of all positive numbers <= n
    congruent to 1(mod3)
    """
    #
    #   1  2  3  4  5  6  7  8  9 10 11 12 13
    #   0  2  2  2  7  7  7 15 15 15 26 26 26
    #
    #  of the form 3n-1, 3n, 3n+1, then the sum is
    #
    #    2 + 5 + 8 + ... + 3n-1
    # 3n-1 +             +  2
    #
    # n*(3n+1)//2
    #
    if n%3==0:
        n -= 1
    elif (n-1)%3==0:
        n -= 2
    return (n+2)*(n+1)//3//2

    
def sievesum_mod3(n):
    """
    Return the sum of primes congruent to 1 and 2 (mod 3)
    less than or equal to n, n//2, n//3, ..., 1
    """
    sqrtn = sqrtInt(n)
    V = [n//i for i in range(1,sqrtn+1)]
    V += list(range(V[-1]-1, 0, -1))
    Vr = list(reversed(V))
    
    S_1mod3 = {v:_f13(v)-1 for v in V}
    S_2mod3 = {v:_f23(v)   for v in V}
    printV(Vr)
    printRuler(Vr)
    printS(Vr, S_1mod3)
    printS(Vr, S_2mod3)
    printRuler(Vr)
    for p in range(2,sqrtn+1):
        p2 = p*p
        if S_1mod3[p] > S_1mod3[p-1]: # p is prime =1(mod3)
            s1m3pm1 = S_1mod3[p-1]
            s2m3pm1 = S_2mod3[p-1]
            for v in V:
                if v < p2: break
                S_1mod3[v] -= p*(S_1mod3[v//p]-s1m3pm1)
                S_2mod3[v] -= p*(S_2mod3[v//p]-s2m3pm1)
            printS(Vr, S_1mod3, tag='%d, %d'%(p, 1))
            printS(Vr, S_2mod3, tag='%d, %d'%(p, 2))
            printRuler(Vr)                                
        elif S_2mod3[p] > S_2mod3[p-1]: # p is prime =2(mod3)
            s1m3pm1 = S_1mod3[p-1]
            s2m3pm1 = S_2mod3[p-1]
            for v in V:
                if v < p2: break
                S_1mod3[v] -= p*(S_2mod3[v//p]-s2m3pm1)
                S_2mod3[v] -= p*(S_1mod3[v//p]-s1m3pm1)
            printS(Vr, S_1mod3, tag='%d, %d'%(p, 1))
            printS(Vr, S_2mod3, tag='%d, %d'%(p, 2))
            printRuler(Vr)
    return {1:S_1mod3, 2:S_2mod3}


def sievecnt_modp(
        n,
        q,
        f=lambda i,q,x: 1,
        F=lambda i,q,x: (x-i)//q+1):
    """
    n: the upper bound
    q: the modulus
    f: the function we want to sum over primes
    F: the sum of f(i) for i from 1 to n

    return: dictionary S where S[i] (for all i in [1,n-1]
            that are relatively prime to n) is a dictionary
            with keys of the form n//j for all positive
            integers j and S[i][n//j] equals the sum of f
            over all primes not dividing the modulus q less
            than or equal to n//j

    By default, f and F are set to give *counts* of primes.
    If you want to find the *sums* of primes in the range,
    then here is how you want to set f and F:
    
    f=lambda i,q,x: x 
    F=lambda i,q,x: ((x-i)//q+1)*(2*i+((x-i)//q)*q)//2
    """
    residues = [i for i in range(1,q) if math.gcd(i,q)==1]
    sqrtn = sqrtInt(n)
    V = [n//i for i in range(1, sqrtn+1)]
    V += list(range(V[-1]-1, 0, -1))
    Vr = list(reversed(V))
    S = {}

    printV(Vr)
    printRuler(Vr)
    S[1] = {v: F(1,q,v)-1 for v in V}
    printS(Vr, S[1])    
    for i in residues[1:]:
        S[i] = {v: F(i,q,v) for v in V}
        printS(Vr, S[i])
    printRuler(Vr)

    invs = {i:modinv(q, i) for i in residues}
    for p in range(2, sqrtn+1):
        p2 = p*p
        for i in residues:
            invi = invs[i] # modinv(q, i)
            if S[i][p] > S[i][p-1]: # p is prime = i(mod p)
                fiqp = f(i,q,p)
                for v in V:
                    if v < p2: break
                    for j in residues:
                        j2 = j*invi % q
                        S[j][v] -= fiqp*(S[j2][v//p] - S[j2][p-1])
                for i in residues:
                    printS(Vr, S[i], tag='%s, %s' % (p, i))
                printRuler(Vr)
    return S
        
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

def sievesumsq(n):
    """
    Return the sum of the squares of all numbers not excluded
    by a sieve of [1..n] by all primes up to
    and including sqrtInt(n).
    """
    sqrtn = sqrtInt(n)
    V = [n//i for i in range(1,sqrtn+1)]
    V += list(range(V[-1]-1,0,-1))
    S = {i:(2*i**3 + 3*i**2 + i)//6-1 for i in V}
    for p in range(2,sqrtn+1):
        if S[p] > S[p-1]:  # p is prime
            sp = S[p-1]  # sum of q**2 for primes q smaller than p
            p2 = p*p
            for v in V:
                if v < p2: break
                S[v] -= p2*(S[v//p] - sp)
    return V, S


def sievesumcb(n):
    """
    Return the sum of the cubes of all numbers not excluded
    by a sieve of [1..n] by all primes up to
    and including sqrtInt(n).
    """
    sqrtn = sqrtInt(n)
    V = [n//i for i in range(1,sqrtn+1)]
    V += list(range(V[-1]-1,0,-1))
    S = {i:((i+2)*(i**3-i)+2*(i**2+i))//4-1 for i in V}
    for p in range(2,sqrtn+1):
        if S[p] > S[p-1]:  # p is prime
            sp = S[p-1]  # sum of q**3 for primes q smaller than p
            p2 = p*p
            p3 = p2*p
            for v in V:
                if v < p2: break
                S[v] -= p3*(S[v//p] - sp)
    return V, S



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
        printS(Vr, S)
    return S
