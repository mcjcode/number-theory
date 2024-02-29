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
    
