#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
from math import sqrt, floor
import numpy as np
from base_complex import infj

from utilities import gcd, euclidean_algorithm

def _rectangle_n_points(nn):
    for ii in range(2*nn) :
        yield (-nn+ii,nn)
        yield (nn,nn-ii)
        yield (nn-ii,-nn)
        yield (-nn,-nn+ii)


def plot_rectangle_n_points():
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(5,5))
    ax = fig.gca()
    for x,y in _rectangle_n_points(5):
        plt.plot(x,y,'.')
    ax.set_xlim([-7,+7])
    ax.set_ylim([-7,+7])
    
    plt.show()

    
def _eisenstein_bound(eps, k, zz):
    """
    Return the number of rectangles of lattice points to compute phi_0,k(z) within eps.
    """

    xx = zz.real
    yy = zz.imag
    
    hh = min(yy, yy/sqrt(xx**2+yy**2))
    
    nmax = int(floor( (8/(hh**(2*k)*(2*k-2)*eps))**(1./(2*k-2)) ))

    return nmax

  
def unrestricted_eisenstein(k,z):
    """
    Return the value G_k(z) of the unrestricted Eisenstein series of weight k at z.
    
    G_k(z) = Sum_{(c,d),c>=0,d>0 if c=0} (1/(cz+d)^(2k)).
    """
    if z.imag <= 0:
        raise ValueError("z not in the upper half plane.")

    Nmax = _eisenstein_bound(1e-08, k, z)
    
    retval = 0.0
    for nn in range(Nmax,0,-1):
        for (cc,dd) in _rectangle_n_points(nn):
            if cc >= 0 and (cc>0 or dd>0):
                vv = 1./(cc*z+dd)**(2*k)
                #print cc, dd, vv
                retval += vv
    return retval

  
def eisenstein(k,z):
    """
    Return the value phi_0,k(z) of the Eisenstein series of weight k at z.
    """

    if z.imag <= 0:
        raise ValueError("z not in the upper half plane.")

    Nmax = _eisenstein_bound(1e-08, k, z)
    
    retval = 0.0
    for nn in range(Nmax,0,-1):
        for (cc,dd) in _rectangle_n_points(nn):
            if cc >= 0 and (cc>0 or dd>0) and gcd(cc,dd)==1 :
                vv = 1/(cc*z+dd)**(2*k)
                #print cc, dd, vv
                retval += vv
    return retval

 
def poincare(k,nu,z):
    """
    Return the value of the Poincare series of weight k and character nu at z.
    """

    if z.imag <= 0:
        raise ValueError("z not in the upper half plane.")
        
    Nmax = _eisenstein_bound(1e-08, k, z)
    
    retval = 0.0
    for nn in range(Nmax,0,-1):
        for (cc,dd) in _rectangle_n_points(nn):
            if cc >= 0 and (cc>0 or dd>0) and gcd(cc,dd)==1 :
                bb, aa = euclidean_algorithm(cc,dd)
                Tz = (aa*z-bb)/(cc*z+dd)
                vv = np.exp(2.0*np.pi*nu*1.0j*Tz)/(cc*z+dd)**(2*k)
                #print cc, dd, vv
                retval += vv
    return retval


rho = np.exp(2.0*np.pi*1.0j/6.0)


class RectangleNPointsTest(unittest.TestCase):
    
    def runTest(self):
        for nn in range(1,10):
            pts = list(_rectangle_n_points(nn))
            npts = len(pts)
            nuniqpts = len(list(set(pts)))
            
            self.assertEqual(npts, 8*nn)
            self.assertEqual(nuniqpts, 8*nn)

class EisensteinTest(unittest.TestCase):
    
    def runTest(self):
        pass

    @staticmethod
    def _random_pt():
        while True:
            x = np.random.rand()-0.5
            y = np.random.rand()
            if x**2+y**2 <= 1 and (x-1)**2+y**2 >= 1 and (x+1)**2 + y**2 > 1:
                break                
        return x + 1.0j*y

    def _test_modularity(self,f):
        z = self._random_pt()
        for kk in [4,5]:
            self.assertAlmostEqual(f(kk, z+1), f(kk,z))
            self.assertAlmostEqual(f(kk,-1/z), z**(2*kk)*f(kk,z))
                
    def test_eisenstein(self):
        """
        Verify that the Eisenstein series are modular forms.
        """
        self._test_modularity(lambda kk, zz: eisenstein(kk,zz))
        self._test_modularity(lambda kk, zz: eisenstein(kk,zz))
        
    def test_unrestricted_eisenstein(self):
        """
        Verify that the unrestricted Eisenstein series are modular forms.
        """
        self._test_modularity(lambda kk, zz: unrestricted_eisenstein(kk,zz))
        self._test_modularity(lambda kk, zz: unrestricted_eisenstein(kk,zz))

        zz = self._random_pt()
        for kk in range(3,7):
            zeta2k = sum( 1./n**(2*kk) for n in range(100000,0,-1) )
            self.assertAlmostEqual(unrestricted_eisenstein(kk, zz), zeta2k*eisenstein(kk,zz),places=6)
    
    def test_poincare(self):
        """
        Verify that the Poincare series are modular forms.
        """
        self._test_modularity(lambda kk, zz: poincare(kk,0,zz))
        self._test_modularity(lambda kk, zz: poincare(kk,1,zz))

def _coset_reps(nn):
    for aa in range(1,nn+1):
        if nn%aa == 0 :
            dd = nn // aa
            for bb in range(dd):
                yield (aa,bb,cc,dd)

def transf(a,b,c,d):
    def f(z):
        if z==infj:
            if c==0:
                return infj
            else:
                return (1.0*a)/c
        if c*z+d==0:
            return infj
        else:
            return (a*z+b)/(c*z+d)
    return f


def hecke_operator(n,k,f):
    """
    Return the application of the nth hecke operator to a function.
    """
    def f(z):
        return sum( (c*z+d)**(2.*k)*transf(a,b,c,d)(z) for (a,b,c,d) in _coset_reps(n) )
    return f

