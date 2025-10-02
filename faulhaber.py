import sympy
from sympy.functions.combinatorial.numbers import binomial, bernoulli


def faulhaber_polynomial(p):
    """
    Return a function implementing the p-th
    faulhaber polynomial F_p, where F_p(n)
    equals the sum of the first n pth powers
    of natural numbers
    """    
    n = sympy.symbols('n')
    Fp = sum(binomial(p+1, r)*bernoulli(r)*n**(p+1-r) for r in range(p+1))/(p+1)
    def f(m):
        return int(Fp.subs({n: m}))
    return f


def generate_polynomials(modulus):
    """
    generate the code for the first fifty faulhaber polynomials
    packaged as a single function F(p, n) in straight python
    without invoking sympy evaluation.  exec the returned
    string.
    """
    n = sympy.symbols('n')    
    retval  = 'modulus = %s\n' % modulus
    retval += 'cs= ['
    for p in range(51):
        Fp = sum(sympy.functions.combinatorial.numbers.binomial(p+1, r)*sympy.functions.combinatorial.numbers.bernoulli(r)*n**(p+1-r) for r in range(p+1))/(p+1)*math.factorial(p+1)
        retval += '%s,\n' % [Fp.coeff(n, i)%modulus for i in range(p+2)]
    retval += '    ]\n'
    retval += 'def F(p, n):\n'
    retval += '    c = cs[p]\n'
    retval += '    retval = c[p+1]\n'
    retval += '    for i in range(p, -1, -1):\n'
    retval += '        retval = (retval*n + c[i]) % modulus\n'
    retval += '    return (retval*pow(math.factorial(p+1), -1, modulus))%modulus\n'
    return retval
