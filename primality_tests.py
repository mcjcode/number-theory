import random
import math


from utilities import (
    is_miller_rabin_witness,
    issq,
)

def isprobprime(p, prob=1.0/(10**12)):
    if p in [2, 3]:
        return True
    
    for _ in range(int(-math.log2(prob))):
        a = random.randint(2, p-1)
        while issq(a):
            a = random.randint(2, p-1)
        if is_miller_rabin_witness(p, a):
            return False
    return True
                