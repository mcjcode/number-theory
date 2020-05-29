def bps_w_rep(n,ps,i=0,only_powerful=False):
    """
    Return all numbers in [1..n] whose prime
    factors are all in the list ps[i:]

    If 'only_powerful' is set to True, then
    only return the numbers that are powerful.
    """
    if len(ps) == i:
        yield 1
    else:
        p = ps[i]
        for pr in bps_w_rep(n,ps,i+1,only_powerful=only_powerful):
            yield pr
        if only_powerful:
            pk = p*p
        else:
            pk = p
        while pk <= n:
            for pr in bps_w_rep(n//pk,ps,i+1,only_powerful=only_powerful):
                yield pk*pr
            pk *= p
    
def bps_facts_w_rep(n,ps,i=0,only_powerful=False):
    """
    Return all factorizations of numbers in [1..n] whose prime
    factors are all in the list ps[i:]

    If 'only_powerful' is set to True, then
    only return the numbers that are powerful.

    Only succeeds if the number of primes in ps[i:]) 
    """
    if len(ps) == i or ps[i] > n or (only_powerful and (ps[i]**2 > n)):
        yield []
    else:
        p = ps[i]
        for pr in bps_facts_w_rep(n,ps,i+1,only_powerful=only_powerful):
            yield pr
        if only_powerful:
            pk = p*p
            e = 2
        else:
            pk = p
            e = 1
        while pk <= n:
            for pr in bps_facts_w_rep(n//pk,ps,i+1,only_powerful=only_powerful):
                yield [(p,e)] + pr
            pk *= p
            e += 1

from heapq import merge
def bps1(n,ps):
    if n == 0 : return
    
    yield 1
    xs = [1]
    pi = 0
    nump = len(ps)
    while pi < nump:
        ys = []
        p = ps[pi]
        if p > n:
            break
        for x in xs:
            if p*x > n :
                break
            yield p*x
            ys.append(p*x)
        xs = list(merge(xs,ys))
        pi += 1

def bps(n,ps):
    if n==0 : return

    yield []
    xs = [(1,[])]
    pi = 0
    nump = len(ps)
    while pi < nump:
        ys = []
        p = ps[pi]
        if p > n:
            break
        for (xprod,x) in xs:
            pxprod = p*xprod
            if pxprod > n:
                break
            yield x+[p]
            ys.append((pxprod,x+[p]))
        xs = list(merge(xs,ys))
        pi += 1

def bps_w_rep(n,ps):
    if n==0 : return
    yield []
    xs = [(1,[])]
    pi = 0
    nump = len(ps)
    while pi < nump:
        ys = []
        p = ps[pi]
        for (xprod,x) in xs:
            pp = p
            e = 1
            while pp <= n:
                pxprod = pp*xprod
                if pxprod > n:
                    break
                yield x+[(p,e)]
                ys.append((pxprod,x+[(p,e)]))
                pp *= p
                e += 1
        xs = list(merge(xs,ys))
        pi += 1

def bps_w_rep_powerful(n,ps):
    if n==0 : return
    yield 1,1,[] # the number the 'squareradical', the factorization
    xs = [(1,1,[])]
    pi = 0
    nump = len(ps)
    while pi < nump:
        ys = []
        p = ps[pi]
        for (xprod,x2rad,x) in xs:
            pp = p**2
            x2rad *= pp
            e = 2
            while pp <= n:
                pxprod = pp*xprod
                if pxprod > n:
                    break
                newfact = x+[(p,e)]
                yield pxprod,x2rad,newfact
                ys.append((pxprod,x2rad,newfact))
                pp *= p
                e += 1
        xs = list(merge(xs,ys))
        pi += 1
