import math
def sievecntsum(n):
    p = int(math.sqrt(n))
    if (p+1)**2 <= n: p += 1
    V = [n//i for i in range(1,p+1)]
    V += list(range(V[-1]-1,0,-1))
    S0 = {i:i-1          for i in V}
    S1 = {i:i*(i+1)//2-1 for i in V}
    SP = {i:i*(i+1)//2-1 for i in V}
    for p in range(2,p+1):
        if S0[p] > S0[p-1]: # p is prime
            p2 = p*p
            for v in V:
                if v < p2: break
                vmodp = v//p
                D0 =  S0[vmodp] - S0[p-1]
                D1 =  S1[vmodp] - S1[p-1]
                S0[v] -=    D0
                S1[v] -= p*(D1)
                SP[v] -= p*(D1-D0)
    return S0[n], S1[n], SP[n]



#print(sievecntsum(10**12) % 10**9)
