import numpy as np
from fractions import Fraction

import numpy as np
from fractions import Fraction

def ch(m) :
    """
    Take a symmetric matrix m and calculate an upper
    triangular matrix a s.t. tr(a)*m*a is diagonal.
    """    
    m = m.copy()
    n = m.shape[0]
    a = np.array([[Fraction(0,1)]*n]*n,dtype=Fraction)
    for i in xrange(n) :
        a[i,i] = Fraction(1,1)
    
    print a
    
    for i in xrange(n-1) :
        for j in xrange(i+1,n) :
            r1 = Fraction(m[j,i],m[i,i])
            r2 = Fraction(m[i,j],m[i,i])
            m[j,:] = m[j,:] - m[i,:] * r1
            m[:,j] = m[:,j] - m[:,i] * r2
            a[:,j] = a[:,j] - a[:,i] * r1
    return m, a

def test_ch() :
    aa = map(lambda L:map(Fraction,L),[[5,3,3],[3,3,3],[3,3,8]])
    mm = np.array(aa)
    m, a = ch(mm)
    assert a.transpose().dot(mm).dot(a) = m

             
