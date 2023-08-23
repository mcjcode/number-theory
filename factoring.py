from math import gcd, log

from utilities import sqrtInt, issq

def fermat_factor(n: int, a=1, b=1) -> tuple:
    if (a*b)%2==0:
        raise ValueError('a or b is even')
    
    if gcd(a, b)>1:
        raise ValueError('a and b are not relatively prime')

    if n%2==0:
        return 2, n//2
    
    g = gcd(a*b, n)
    if g>1:
        return g, n//g

    abn = a*b*n
    x = sqrtInt(abn)
    if x*x<abn:
        x += 1
    z = x*x - abn
    while True:
        if issq(z):
            break
        z += 2*x+1
        x += 1

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

