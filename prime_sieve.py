#!/usr/bin/env python
#
# prime_sieve.py
#

import math
from utilities import sqrtInt, timeit
import numpy as np

def zeroit(arr, starti, stepi):
    arr[starti::stepi] = 0

def segmented_sieve(n, trace=False):
    """
    Return the list of all prime numbers
    less than or equal to n.  Runtime
    complexity of O(n) (or a bit better)
    with space complexity of O(sqrt(n)).
    """
    seglen = sqrtInt(n+1)

    # always make seglen a multiple of 6
    seglen = max(6, 6*(seglen//6))+6 # add 6 so that the first segment contains sqrt(n)
    
    wheel = [4, 2]
    wheel_len = len(wheel)
 
    if trace:
        print('segment length = %d' % (seglen, ))
    nsegs  = (n+1)//seglen
    if seglen*nsegs < n+1:
        nsegs += 1
    if trace:
        print('number of segments = %d' % (nsegs, ))
    ps = [2, 3]
    for p in ps:
        yield p
    segi = 0
    seg_start = 0
    seg = np.ones((seglen, ), dtype=np.int8)
    seg[:4] = [0, 0, 1, 1]
    
    while segi < nsegs:
        seg_end   = seg_start+seglen-1 
        
        if trace: 
            print('segment number = %d' % (segi))
            print('segment start = %d' % (seg_start, ))
            print('segment end = %d' % (seg_end, ))

        ub = math.sqrt(seg_end)
        for p in ps:
            if p > ub:
                break
            #seg[(p-seg_start)%p::p] = 0
            zeroit(seg, (p-seg_start)%p, p)
            if trace:
                print(p, '%s'%(seg, ))
        
        starti = 1
        endi = min(seglen, ub-seg_start+1)
        stepi = 2

        i = starti
        wheeli = 0
        while i < endi: 
            if seg[i]:
                p = seg_start + i
                yield p
                if segi==0:
                    ps.append(p)
                #seg[p*p-seg_start::p]
                zeroit(seg, p*p-seg_start, p)
                if trace:
                    print(p, '%s'%(seg, ))
            i += wheel[wheeli] #stepi
            wheeli = (wheeli + 1) % wheel_len

        maxi = min(seglen, n-seg_start+1)
        while i < maxi:
            if seg[i]: # seq_start+i is prime
                p = seg_start+i
                yield p
                if segi==0:
                    ps.append(p) 
            #i += 2
            i += wheel[wheeli]
            wheeli = (wheeli + 1) % wheel_len
        
        segi += 1
        seg_start += seglen
        seg[:] = 1

    #return ps

def indicator(n):
    nm1 = n-1
    while True:
        yield True
        i=0
        while i<nm1:
            yield False
            i += 1

def sqfree_parts():
    """
    For each positive integer n,
    yield the list of primes dividing
    n, w/o multiplicity.  Horribly
    inefficient.  Just a proof of
    concept.
    """
    yield [] # 1 has no prime divisors
    indicators = []
    i = 2
    while True:
        ps = []
        for (p, pi) in indicators:
            if next(pi):
                ps.append(p)
        if not ps:
            ps = [i]
            g = indicator(i) 
            _ = next(g)
            indicators.append( (i, g) )
        yield ps
        i += 1

def primeFactors(n):
    """
    Return a length n+1 list pf such that pf[k] is the list of primes dividing k
    """
    retval = [[] for _ in range(n+1)]
    for p in range(2, n+1):
        if not retval[p]:  # p is prime
            for j in range(p, n+1, p):
                retval[j].append(p)
    return retval
