#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

from math import atan, sin, cos, exp
import numpy as np

import matplotlib

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path

from utilities import gcd, euclidean_algorithm
from base_complex import infj

import unittest

def c2xy(z):
    """
    Given a (finite) complex number, return the associated ordered pair
    """
    return z.real, z.imag


def Sz(z):
    """
    Apply the linear fractional transformation z -> -1/z to a point.
    """
    if z == infj:
        return complex(0.0, 0.0)
    else:
        if z == complex(0.0, 0.0):
            return infj
        else:
            return -1.0/z


def S(corners):
    """
    Apply the linear fractional transformation z -> -1/z to a sequence of points.
    """
    return list(map(Sz, corners))


def Tz(z):
    if z == infj:
        return infj
    else:
        return z + 1.0+0.0j


def T(corners):
    return list(map(Tz, corners))


def Uz(z):
    if z == infj:
        return infj
    else:
        return z - 1.0+0.0j


def U(corners):
    return list(map(Uz, corners))


def z1z2_to_pts(z1, z2, n):
    """
    Return a sequence of n+1 points along the geodesic from z1 to z2.
	  """
    x1, y1 = c2xy(z1)
    x2, y2 = c2xy(z2)

    # this was an equality test, which proved unwise when
    # used with floats.  Taking as a sign that I should be
    # doing this all over Q[\sqrt(3)].
    #
    if abs(x1 - x2) < 0.00000001:
        dy = (y2 - y1) / n
        return [(x1, y1 + ii * dy) for ii in range(n + 1)]
    #
    # compute the point on the real axis
    # equidistant from z1 and z2
    #
    x = (abs(z1) ** 2 - abs(z2) ** 2) / (2.0 * (x1 - x2))
    radius = np.sqrt((x1 - x) ** 2.0 + y1 ** 2.0)

    if y1 == 0 and y2 == 0:
        if x1 > x2:
            theta1 = 0.0
            theta2 = np.pi
        else:
            theta1 = np.pi
            theta2 = 0.0

    if y1 == 0:
        if x1 < x2:
            theta1 = np.pi
        else:
            theta1 = 0.0
    else:
        theta1 = atan(y1 / (x1 - x))
        if theta1 < 0:
            theta1 = np.pi + theta1

    if y2 == 0:
        if x2 < x1:
            theta2 = np.pi
        else:
            theta2 = 0.0
    else:
        theta2 = atan(y2 / (x2 - x))
        if theta2 < 0:
            theta2 = np.pi + theta2

    dtheta = (theta2 - theta1) / n
    return [(x + cos(theta1 + ii * dtheta) * radius, sin(theta1 + ii * dtheta) * radius) for ii in range(n + 1)]


def tile_patch(corners, **kwdargs):
    """
    Given the three corners of a tile (a fundamental
    region for G, the full modular group), return the
    matplotlib Patch artist that can be added to a
    plot.
    
    Other than the corners of the tile, other keyword
    arguments are passed through to Patch object, and
    so its color, transparency, edge properties can
    all be specified here.
    """
    z1, z2, z3 = corners
    k = 20

    if z1 == infj or z2 == infj or z3 == infj:
        if z2 == infj:
            z1, z2, z3 = z2, z3, z1
        elif z3 == infj:
            z1, z2, z3 = z3, z1, z2
        pth = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO]
        x2, y2 = c2xy(z2)
        x3, y3 = c2xy(z3)
        pts = [(x2, y2), (x2, 2.0), (x3, 2.0), (x3, y3)]

        pts += z1z2_to_pts(z3, z2, 2 * k + 1)[1:]
        pth += ([Path.CURVE3] * (2 * k)) + [Path.CLOSEPOLY]
        verts, codes = pts, pth
    else:
        pth1 = [Path.MOVETO] + ([Path.CURVE3] * (2 * k)) + [Path.LINETO]
        pth2 = ([Path.CURVE3] * (2 * k)) + [Path.LINETO]
        pth3 = ([Path.CURVE3] * (2 * k)) + [Path.CLOSEPOLY]
        pts1 = z1z2_to_pts(z1, z2, 2 * k + 1)
        pts2 = z1z2_to_pts(z2, z3, 2 * k + 1)
        pts3 = z1z2_to_pts(z3, z1, 2 * k + 1)
        verts, codes = pts1 + pts2[1:] + pts3[1:], pth1 + pth2 + pth3

    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, **kwdargs)
    return patch

rho = np.exp(2*np.pi*1.0j/6.0)
D = (rho, rho**2, 0.0+0.0j)

def parse_to_mat(ss):
    M = np.array([[1,0],[0,1]],dtype=int)
    for c in ss:
        if c == 'S':
            A = np.array([[0,-1],[1,0]],dtype=int)
        elif c == 'T':
            A = np.array([[1,1],[0,1]],dtype=int)
        elif c == 'U':
            A = np.array([[1,-1],[0,1]],dtype=int)
        elif c == 'I':
            A = np.array([[1,0],[0,1]],dtype=int)
        M = M.dot(A)
    return M

def mat_to_fcn(M):
    def ffz(z):
        if z == infj:
            if M[1,0]==0:
                return infj
            else:
                return float(M[0,0])/float(M[1,0])
        num = M[0,0]*z+M[0,1]
        den = M[1,0]*z+M[1,1]
        if den==0:
            return infj
        else:
            return num/den
    def ff(zs):
        return map(ffz,zs)
    return ff

def parse_word(ss):
    M = parse_to_mat(ss)
    return mat_to_fcn(M)
    

def invert_letter(letter):
    return {'S':'S', 'T':'U', 'U':'T', 'I':'I'}[letter]


def invert_word(ss):
    return ''.join(invert_letter(s) for s in reversed(ss))


def coset_reps(qq):
    """
    Yield coset representatives for Gamma(q) in SL(2,Z).
    """
    
    for cc in range(qq):
        for dd in range(qq):
            if gcd(cc,gcd(dd,qq))==1:
                ii = 0
                if cc==0:
                    cc=qq
                while True:
                    if gcd(cc,dd+ii*qq) == 1:
                        dd += ii*qq
                        break
                    elif gcd(cc,dd-ii*qq) == 1:
                        dd -= ii*qq
                        break
                    ii += 1
                for aa in range(qq):
                    for bb in range(qq):
                        quot, rem = divmod(aa*dd-bb*cc,qq)
                        if rem==1:
                            for (ee,ff) in [(ee,ff) for ee in range(2*qq) for ff in range(2*qq)]:
                                if (aa+ee*qq)*dd-(bb+ff*qq)*cc==1:
                                    yield np.array([[aa+ee*qq,bb+ff*qq],[cc,dd]],dtype=int)
                                    break
                                elif (aa+ee*qq)*dd-(bb-ff*qq)*cc==1:
                                    yield np.array([[aa+ee*qq,bb-ff*qq],[cc,dd]],dtype=int)
                                    break                                    
                                elif (aa-ee*qq)*dd-(bb+ff*qq)*cc==1:
                                    yield np.array([[aa-ee*qq,bb+ff*qq],[cc,dd]],dtype=int)
                                    break                                                                            
                                elif (aa-ee*qq)*dd-(bb-ff*qq)*cc==1:
                                    yield np.array([[aa-ee*qq,bb-ff*qq],[cc,dd]],dtype=int)
                                    break                                    

def name_to_latex(name):
    retval = '$'
    i = 0
    while i < len(name):
        cc = name[i]
        j = 1
        while i+j < len(name):
            if name[i+j]==cc:
                j += 1
            else:
                break
        if cc=='T':
            if j==1:
                retval += 'T'
            else:
                retval += 'T^{%d}' % (j,)
        elif cc=='U':
            retval += 'T^{%d}' % (-j,)
        elif cc=='S':
            if j%2==1:
                retval += 'S'
        else: # cc=='I'
            pass
        i += j
    retval += '$'
    if retval=='$$':
        retval='$I$'
    return retval


def plot_regions(tile_names, center, shift_name, transform_names):
    """
    Plot a fundamental domain and neighbors for a subgroup G of SL(2,Z).
    
    tile_names - coset representatives for G in SL(2,Z)
    center     - complex number for the center of the fundamental region.
    shift_name - the name of an element g of G, g(D) will be
                 the center of the plot
    transform_names - the name of elements of G for neighboring
                 cells to plot
    """
    rho = np.exp(2*np.pi*1.0j/6.0)
    D = (rho, rho**2, 0.0+0.0j)
    
    shift = parse_word(shift_name)
    shift_inv = parse_word(invert_word(shift_name))
    
    D2 = shift(D)
    shift_inv_name = invert_word(shift_name)
    
    #transfs = [lambda tt: shift(parse_word(tname)(shift_inv(tt))) for tname in tile_names]
    transfs = [parse_word( shift_name+tname+shift_inv_name ) for tname in tile_names]
    tiles = [transf(D2) for transf in transfs]

    names = [shift_name+name+shift_inv_name for name in transform_names]
    transforms = map(parse_word, names)
    
    
    center2 = shift([center])[0]
    vert_xcoords = []
    for f, name in zip(transforms, names):
        for tile in map(f, tiles):
            vert_xcoords += [zz.real for zz in tile]
    xmax = np.amax(np.array(vert_xcoords))
    xmin = np.amin(np.array(vert_xcoords))
    
    x_to_y_ratio = (xmax-xmin)/1.5
    
    fig = plt.figure(figsize=(x_to_y_ratio*3,3))
    ax = fig.gca()
                        
    for f, c, name in zip(transforms, 'rgbccgbrgbccgbr', names) :
        for tile in map(f, tiles) :
            ax.add_patch(tile_patch(tile, facecolor=c))
            x, y = c2xy(f([center2])[0])
            latex_name = name_to_latex(name)
            ax.text(x, y, latex_name,color='w',
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=9)

    ax.grid()
    #ax.set_aspect('equal')
    
    
    ax.set_xlim([xmin,xmax])
    #ax.set_xlim([-2.0,+2.0])
    ax.set_ylim([0.0, 1.5])

    plt.show()

def plot_gamma_2():
    """
    Plot the fundamental region for Gamma(2), the principal
    congruence subgroup of level 2, and neighboring regions
    """
    plot_regions(['I', 'T', 'S', 'TS', 'TSTS', 'TST'],
                  np.exp(2*np.pi*1.0j/6.0),
                  'SUUS', 
                  ['I','TT','UU','STTS','TTSTTS','TSTTSU','SUUS'])

def plot_gamma_3():
    plot_regions(['I'],0.75j,'I',['I','T','U','S','TS','US','ST','STS','USUUS','TSTS','TSTTS','TST'])

class CosetRepsUnitTest(unittest.TestCase):
    def coset_reps_test(self):
        for q in [2,3,5,7,11]:
            nn = len(list(coset_reps))
            self.assertEqual(nn, (q**3-q))

