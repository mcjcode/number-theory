from utilities import issq, sqrtInt

def cornacchias_algorithm(d, m, rs):
    """
    Return a primitive solution (x, y) to the Diophantine
    equation x^2 + dy^2 = m, where 1<=d<m and d
    and m are coprime. 'rs' is the list of square
    roots of -d (mod m).
    """
    for r1 in rs:
        if (r1**2+d)%m==0:
            r0 = m        
            if 2*r1 >= m:
                r1 = m - r1
            assert r1 < m/2
        
            while r1**2 >= m:
                r0, r1 = r1, r0 % r1
        
            if (m-r1**2)%d==0 and issq((m-r1**2)//d):
                x = r1
                y = sqrtInt((m-r1**2)//d)
                return x, y
    return None