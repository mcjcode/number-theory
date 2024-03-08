#!/usr/bin/env python -i
#
# fourier.py
#

def generalized_hadamard_transform(a, m):
    """
    m is a k-by-k matrix
    a is a 1d np.array with length a positive power of k
    """
    k = m.shape[0]
    powkn = len(a) # 16
    powk = 1
    while powk < powkn: # 1, 2, 4, 8
        for j in range(0, powkn, k*powk):
            for i in range(0, powk):
                a[j+i:j+i+k*powk:powk] = m.dot(a[j+i:j+i+k*powk:powk])
            #a[j:j+k*powk:powk] = m.dot(a[j:j+k*powk:powk])
        powk *= k
        
def hadamard_transform(a):
    """
    Do an in-place 'fast' Walsh-Hadamard transform.
    len(a) must be a power of 2
    """
    pow2n = len(a)
    pow2i = 1
    while pow2i < pow2n:
        for j in range(0, pow2n, 2*pow2i):
            for k1 in range(j, j+pow2i):
                k2 = k1+pow2i
                a[k1], a[k2] = a[k1]+a[k2], a[k1]-a[k2]
        pow2i *= 2

def or_transform(a):
    pow2n = len(a)
    pow2i = 1
    while pow2i < pow2n:
        for j in range(0, pow2n, 2*pow2i):
            for k1 in range(j, j+pow2i):
                k2 = k1+pow2i
                a[k2] = a[k1]+a[k2]
        pow2i *= 2

def inverse_or_transform(a):
    pow2n = len(a)
    pow2i = 1
    while pow2i < pow2n:
        for j in range(0, pow2n, 2*pow2i):
            for k1 in range(j, j+pow2i):
                k2 = k1+pow2i
                a[k2] = -a[k1]+a[k2]
        pow2i *= 2
    
