import pytest


from real_quadratic_fields import pell


def test_pell():
    for d in range(2, 10000):
        a, b = pell(d**4+1)
        assert a*a-(d**4+1)*b*b==1
            
