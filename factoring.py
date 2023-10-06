import time
import numpy as np
import random
from collections import Counter
import math; from math import gcd, log, sqrt, exp
from utilities import timeit, prod, sqrtInt, cbrtInt, issq, digits, trial_division, isprime, is_miller_rabin_witness
from prime_sieve import segmented_sieve
from primality_tests import isprobprime
from quadratic_extensions import legendre
from mod2linalg import solve_linear_equations_mod_2, null_space_mod_2, det_mod_2
import random


def factors(n):
    f, _ = trial_division(n)
    retval = [1]
    for p, e in f:
        pi = 1
        tmp = []
        for i in range(e+1):
            for fact in retval:
                tmp.append(fact*pi)
            pi *= p
        retval = tmp
    return retval


def fermat_factor(n: int, a=1, b=1) -> tuple:
    """
    Looking for a factorization of an odd number n=u*v
    where the ratio u:v is around a:b (when a=b=1, this
    just becomes Fermat's original method)
    """
    g = gcd(a, b)
    a //= g
    b //= g

    if (a*b)%2==0:
        raise ValueError('a or b is even')
    
    if n%2==0:
        return 2, n//2

    assert n%2==1
    
    g = gcd(a*b, n)
    if g>1:
        return g, n//g

    abn = a*b*n
    x = sqrtInt(abn)
    if x*x<abn:
        x += 1
    z = x*x - abn
    xmy_prev = x - math.sqrt(z)
    while True:
        if issq(z):
            break
        z += 2*x+1
        x += 1
        xmy = x - math.sqrt(z)
        if xmy_prev - xmy < trial_division_threshold:
            facts = trial_division(n, bound=xmy_prev)
            return next(facts)
        
    y = sqrtInt(z)

    assert (x-y)*(x+y) == a*b*n

    u, v = x-y, x+y
    # now clear the factors of a and b from u and v
    
    g = gcd(a, u)
    u //= g
    a //= g

    g = gcd(b, u)
    u //= g
    b //= g

    g = gcd(a, v)
    v //= g
    a //= g

    g = gcd(b, v)
    v //= g
    b //= g

    assert u*v == n
    
    return u, v


def pollard_p_minus_1(n, B, ps=[2, 3, 5]):
    a = 2
    for p in ps:
        if p>B:
            break
        k = int(log(B)/log(p))
        a = pow(a, p**k, n)
        g = gcd(a-1, n)
        if 1<g<n:
            return g
    return False


def williams_p_plus_1(N: int, B: int, A: int) -> int:
    """
    If N has a prime factor p such that p+1 is
    smooth, then p will divide the numbers V_n - 2
    where n is a multiple of p+(D|p) where V is a
    Lucas sequence and D=A**2-4

    :param N: the number we are attempting to factor
    :param B: the smoothness bound (we will check
              the gcd of N and the terms V_{1}, V_{2!},
              ... V_{B!}.
    :param A: a parameter on which the Lucas sequence
              depends.  we want to choose A such that
              (D|p)=1, but we don't know the p yet, so
              we may have to try several values of A.
    :return:  either a positive integer which is a
              proper factor of N, or 0, if no such
              factor was found.

    N = 2**67-1 # = 147573952589676412927
    for A in 9, 17, 23:
        result = williams(N, B=10000, A)
        if result and result!=N:
            print(result)
            break

    will find the factor 193707721 of M67, a Mersene number
    that Mersene reported (incorrectly) as prime, and whose
    factors were presented by Frank Cole in 1903.
    """
    D = A*A-4
    v = A
    for M in range(2, B+1):
        g = gcd(N, v-2)
        if g>1 and g!=N:
            return g
        #
        # https://rieselprime.de/ziki/P%2B1_factorization_method shows that
        # we can compute V_{nm} from V_m in O(log(n)) time, without
        # computing m or m*n
        #
        x = v
        y = (v*v-2)%N
        bits = list(reversed(list(digits(M, base=2))[:-1]))
        for bit in bits:
            if bit==1:
                x=(x*y-v)%N
                y=(y*y-2)%N
            else:
                y=(x*y-v)%N
                x=(x*x-2)%N
        v=x
    return 0

def lehman(n):
    if n<=21:
        return trial_division(n)
    
    fs, n = trial_division(n, bound=cbrtInt(n))
    if n==1:
        return fs, 1
    
    if isprobprime(n):
        return fs + [(n, 1)], 1
    
    _nsixth = int(n**(1/6))
    
    for k in range(1, cbrtInt(n)+1):
        _4kn = 4*k*n
        lb = sqrtInt(_4kn)
        if lb*lb != _4kn:
            lb += 1
        
        for a in range(lb, lb + _nsixth//(4*sqrtInt(k))+1):
            b2 = a*a - _4kn
            if issq(b2):
                b = sqrtInt(b2)
                g = gcd(a+b, n)
                if isprobprime(g):
                    f1, r1 = [(g, 1)], 1
                else:
                    f1, r1 = trial_division(g)
                    
                if isprobprime(n//g):
                    f2, r2 = [(n//g, 1)], 1
                else:
                    f2, r2 = trial_division(n//g)
                fs = Counter(dict(fs))
                f1 = Counter(dict(f1))
                f2 = Counter(dict(f2))
                f = sorted(list((fs+f1+f2).items()))
                
                return f, r1*r2


def hart_one_line(n):
    if issq(n):
        sqrtn = sqrtInt(n)
        f, r = trial_division(sqrtn)
        if r==1:
            return [(p, 2*e) for p, e in f]

    
    fs, n = trial_division(n, bound=cbrtInt(n))
    if n==1:
        return fs
    else:
        if isprobprime(n):
            return fs+[(n, 1)]
            
        ni = n
        while True:
            x = sqrtInt(ni)  ## hot spot I
            if x*x < ni:
                x += 1
            x2 = (x*x)%n
            if issq(x2):
                t = sqrtInt(x2)  ## hot spot II
                g = gcd(x-t, n)
                if g > 1 and g < n:
                    if isprobprime(g):
                        f1, r1 = [(g,1)], 1
                    else:
                        f1, r1 = trial_division(g)
                    if isprobprime(n//g):
                        f2, r2 = [(n//g, 1)], 1
                    else:
                        f2, r2 = trial_division(n//g)
                        
                    if r1==1 and r2==1:
                        fs = Counter({p:e for p, e in fs})
                        f1 = Counter({p:e for p, e in f1})
                        f2 = Counter({p:e for p, e in f2})
                        return sorted(list((fs+f1+f2).items()))
                    else:
                        break
            ni += n


def pollard_rho(n, fxn = lambda x, n: (x*x+1) % n):
    x = 2
    y = 2
    g = 1

    f = lambda x: fxn(x, n)
    
    while g==1:
        x = f(x)
        y = f(f(y))
        g = gcd(abs(x-y), n)
    
    return g
    
    
def compare_times(start, count, method='pollard_rho'):
    xs = []
    start += start%2==0
    for p in range(start, start+count, 2):
        if not isprobprime(p):
            xs.append(p)


    for x in xs:
        print(x, end=' ')
        
        if method=='pollard_rho':
            t0 = time.time()
            fs, r = trial_division(x, cbrtInt(x))

            if r>1:
                if isprobprime(r):
                    t1 = time.time()
                    print('%10.3f seconds ' % (t1-t0), fs + [(r, 1)])
                    continue
                y = pollard_rho(r)
                if y!=r:
                    if isprobprime(y):
                        f1, r1 = [(y, 1)], 1
                    else:
                        f1, r1 = trial_division(y)

                    if isprobprime(r//y):
                        f2, r2 = [(r//y, 1)], 1
                    else:
                        f2, r2 = trial_division(r//y)

                    if r1==1 and r2==1:
                        fs = Counter({p:e for p, e in fs})
                        f1 = Counter({p:e for p, e in f1})
                        f2 = Counter({p:e for p, e in f2})
                        f = sorted(list((fs+f1+f2).items()))
                        t1 = time.time()
                        print('%10.3f seconds*' % (t1-t0), f)
                else:
                    print('FAIL')
            else:
                t1 = time.time()
                print('%10.3f seconds ' % (t1-t0), fs)
        
        elif method=='hart_one_line':
            t0 = time.time()
            y = hart_one_line(x)
            t1 = time.time()
            print('%10.3f seconds' % (t1-t0), y)
                  
        elif method=='lehman':
            t0 = time.time()
            y, r = lehman(x)
            assert r==1
            t1 = time.time()
            print('%10.3f seconds' % (t1-t0), y)
                    

def smooth_trial_division(n, ps=None):
    ans = []
    for p in ps:
        e = 0
        while n%p == 0:
            e += 1
            n //= p
        if e:
            ans.append((p, e))
    return ans, n


def cont_frac(P, Q, n):
    # Q must divide n-P*P
    A0 = 1
    q0 = sqrtInt(n)
    q, A = q0, q0
    yield P, Q, q, A
    while True:
        assert (n-P*P) % Q == 0
        P = Q*q-P
        Q = (n-P*P)//Q
        q = (P+q0)//Q
                
        A, A0 = (q*A + A0) % n, A
        yield P, Q, q, A
        if q==2*q0:
            break


def dixon(n, trace=False):
    """
    Implementation of Dixon's toy factoring method from
    
    Dixon, J. D. (1981). "Asymptotically fast factorization of integers"
    Math. Comp. 36 (153): 255–260.
    """

    bound = int(exp(sqrt(log(n)*log(log(n))) / 2))
    ps = [p for p in segmented_sieve(bound) if legendre(n, p)==1]
    pi = {p:i for i, p in enumerate(ps)}
    nps = len(ps)
    print(bound)
    print(ps)
    nQs_considered = 0
    vecs, Qs, As, cnt = [], [], [], 0
    while True:
        nQs_considered += 1
        A = random.randint(1, n)
        Q = A*A % n
        if Q==0:
            # A is in [1,n) but n divides A*A.
            return(gcd(A,n))
        f, r = smooth_trial_division(Q, ps)
        #print('factoring', Q, f, r)
        if r==1 and Q>1 and Q!=A**2:
            vec = [0]*nps
            for p, e in f:
                vec[pi[p]] = e%2
            vecs.append(vec)
            Qs.append(Q)
            As.append(A)
            #print(Q, vec, A, f)
            cnt += 1
            if cnt>nps:
                break
    Amat = np.array(vecs, dtype=np.int8)     
    N = null_space_mod_2(Amat.transpose())
    while True:
        col = N.dot(np.random.randint(0, 2, (N.shape[1],)))%2
        Q2s = [Qs[i] for i in range(len(Qs)) if col[i]]
        A2s = [As[i] for i in range(len(As)) if col[i]]
        x = sqrtInt(prod(Q2s))%n
        y = prod(A2s)%n
        g = gcd(x+y, n)

        if n%g == 0 and g!=1 and g!=n:
            g2 = gcd(abs(x-y), n)
            if trace:
                maxQdigits = max(len(str(Q)) for Q in Qs); fmtQstr = '%'+str(maxQdigits)+'d'
                maxAdigits = max(len(str(A)) for A in As); fmtAstr = '%'+str(maxAdigits)+'d'
                maxPdigits = max(len(str(p)) for p in ps); fmtPstr = '%'+str(maxPdigits)+'d'
                def formatVec(v):
                    return ' '.join((fmtPstr%x if x else ' '*maxPdigits) for x in v)
                star = ' *'
                table = ' ' * (2+maxQdigits) + formatVec(ps) + '       \n'
                for i in range(Amat.shape[0]):
                    table += star[col[i]] + ' ' +  (fmtQstr % Qs[i]) + formatVec(list(Amat[i])) + ' ' + (fmtAstr % As[i]) + '\n'
                table += str(x) + '^2 = ' + str(y) + '^2 mod ' + str(n) + '\n'
                print('Q\'s considered = ', nQs_considered)
                print('number of eligible Q/A pairs = ', len(Qs))
                print('number of primes in factor base = ',nps)
                print(table)
                print('factors', g, g2, 'found')
            return g, g2


def cfrac(n, trace=False):
    if isprobprime(n):
        if trace:
            print(n, 'is prime')
        return None
    
    orign = n
    while True:
        if trace:
            print('n =', n)

        if issq(n):
            return sqrtInt(n)
        
        bound = int(exp(sqrt(log(n)*log(log(n))) / 2))
        ps = [p for p in segmented_sieve(bound) if legendre(n, p)==1]
        pi = {p:(i+1) for i, p in enumerate(ps)}

        nps = len(ps) + 1

        cnt, vecs, Qs, As, prevA = 0, [], [], [], 1
        for i, (P, Q, q, A) in enumerate(cont_frac(0, 1, n)):
            f, r = smooth_trial_division(Q, ps)
            if r==1 and Q>1 and Q!=(prevA**2):
                vec = [0]*nps
                vec[0] = i%2
                for p, e in f:
                    vec[pi[p]] = e%2
                vecs.append(vec)
                Qs.append(Q)
                As.append(prevA)
                cnt += 1
                if cnt>nps:
                    break
            prevA = A

        nQs_considered = i

        Amat = np.array(vecs, dtype=np.int8)
        if Amat.shape[0]==0:
            n += orign
            continue
            
        N = null_space_mod_2(Amat.transpose())

        def formatVec(v):
            return ' '.join(('%6d'%x if x else '      ') for x in v)

        for kk in range(1):
            col = N.dot(np.random.randint(0, 2, (N.shape[1],)))%2
            Q2s = [Qs[i] for i in range(len(Qs)) if col[i]]
            A2s = [As[i] for i in range(len(As)) if col[i]]
            x = sqrtInt(prod(Q2s))%n
            y = prod(A2s)%n
            g = gcd(x+y, orign)

            if orign%g == 0 and g!=1 and g!=orign:
                g2 = gcd(abs(x-y), orign)
                if trace:
                    star = ' *'
                    table = '         ' + formatVec([-1] + ps) + '       \n'
                    for i in range(Amat.shape[0]):
                        table += star[col[i]] + ' ' + formatVec([Qs[i]] + list(Amat[i]) + [As[i]]) + '\n'
                    table += str(x) + '^2 = ' + str(y) + '^2 mod ' + str(n) + '\n'
                    print('Q\'s considered = ', nQs_considered)
                    print('number of eligible Q/A pairs = ', len(Qs))
                    print('number of primes in factor base = ',nps)
                    print('number of null-vectors attempted = ',kk+1)
                    print('l = ', n//orign)
                    print(table)
                    print('factors', g, g2, 'found')
                return g, g2
                break
            else:
                n += orign
                continue