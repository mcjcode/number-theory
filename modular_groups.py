#!/usr/bin/env python -i
# -*- coding: utf-8 -*-

import unittest
import itertools
from math import atan, sin, cos, exp, floor, sqrt
import numpy as np

try:
    import matplotlib
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
except:
    pass

from utilities import (
    gcd,
    euclidean_algorithm,
    factorize,
    isprime,
)
    
from base_complex import infj

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


def halfplane_to_poincare_disk(zz):
    """
    Send zs —> (zs-i)/(-i*zs+1) = (i*zs+1)/(zs+i) = i + 2/(zs+i)
    """
    if zz==-1.0j:
        return infj
    elif zz==infj:
        return 1.0j
    else:
        return (zz-1.0j)/(-1.0j*zz+1)


def poincare_disk_to_halfplane(ww):
    """
    send ws -> (ws+i)/(i*ws+1) = (-i*ws+1)/(ws-i) = -i + 2/(ws-i)
    """
    if ww==1.0j:
        return infj
    elif ww==infj:
        return -1.0j
    else:
        return (ww+1.0j)/(1.0j*ww+1)


def z1z2_to_pts(z1, z2, n, model):
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
        hp_pts = [(x1, y1 + ii * dy) for ii in range(n + 1)]
    else:
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
        hp_pts = [(x + cos(theta1 + ii * dtheta) * radius, 
                  0 + sin(theta1 + ii * dtheta) * radius) for ii in range(n + 1)]

    if model=='halfplane':
        return hp_pts
    else:
        zs = [xi+1.0j*yi for (xi,yi) in hp_pts]
        zs = map(halfplane_to_poincare_disk,zs)
        return [(zz.real,zz.imag) for zz in zs]


def tile_patch(corners, model='halfplane', **kwdargs):
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

    k = 30

    z1, z2, z3 = corners        
    if z1 == infj or z2 == infj or z3 == infj:
        if z2 == infj:
            z1, z2, z3 = z2, z3, z1
        elif z3 == infj:
            z1, z2, z3 = z3, z1, z2
        assert(z1 == infj); assert(z2 != infj); assert(z3 != infj)
        
        if model=='disk':
            
            w2 = halfplane_to_poincare_disk(z2)
            #print 'w2 = ', w2
            s, t = w2.real, w2.imag
            cx, cy = (s**2+(t-1)**2)/(2*s), 1.0
            rr = abs(cx)
            if s<0:
                th1 = -atan(-(t-cy)/(s-cx))
                th2 = 0.0
                while th1 > 0:
                    th1 -= np.pi
            else: # s>0
                th1 = atan((t-cy)/(s-cx))
                th2 = np.pi
                while th1 < np.pi:
                    th1 += np.pi
            #print '(cx,cy)=(%f,%f)' % (cx,cy)
            #print 'th1 = ', th1*180/np.pi
            #print 'th2 = ', th2*180/np.pi
            dth = (th2-th1)/(2*k+1)
            pts = [(cx+rr*cos(th),cy+rr*sin(th)) for th in np.arange(th1,th2+0.5*dth,dth)]
            pth = [Path.MOVETO] + [Path.CURVE3]*(len(pts)-1)
            assert(len(pts)==len(pth))
            
            w3 = halfplane_to_poincare_disk(z3)
            #print 'w3 = ', w3
            s, t = w3.real, w3.imag
            cx, cy = (s**2+(t-1)**2)/(2*s), 1.0
            rr = abs(cx)
            if s<0:
                th1 = 0.0
                th2 = -atan(-(t-cy)/(s-cx))
                while th2 > 0:
                    th2 -= np.pi
            else: # s>0
                th1 = np.pi
                th2 = atan((t-cy)/(s-cx))
                while th2 < np.pi:
                    th2 += np.pi
            #print '(cx,cy)=(%f,%f)' % (cx,cy)
            #print 'th1 = ', th1*180/np.pi
            #print 'th2 = ', th2*180/np.pi
            dth = (th2-th1)/(2*k)
            pts1 = [(cx+rr*cos(th),cy+rr*sin(th)) for th in np.arange(th1,th2+0.5*dth,dth)]
            pts += pts1
            pth += [Path.CURVE3]*len(pts1)
            assert(len(pts)==len(pth))
        else:
            pth = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO]
            x2, y2 = c2xy(z2)
            x3, y3 = c2xy(z3)
            pts = [(x2, y2), (x2, 1.5), (x3, 1.5), (x3, y3)]

        pts += z1z2_to_pts(z3, z2, 2 * k + 1, model)[1:]
        pth += ([Path.CURVE3] * (2 * k)) + [Path.CLOSEPOLY]
        verts, codes = pts, pth
    else:
        pth1 = [Path.MOVETO] + ([Path.CURVE3] * (2 * k)) + [Path.LINETO]
        pth2 = ([Path.CURVE3] * (2 * k)) + [Path.LINETO]
        pth3 = ([Path.CURVE3] * (2 * k)) + [Path.CLOSEPOLY]
        pts1 = z1z2_to_pts(z1, z2, 2 * k + 1, model)
        pts2 = z1z2_to_pts(z2, z3, 2 * k + 1, model)
        pts3 = z1z2_to_pts(z3, z1, 2 * k + 1, model)
        verts, codes = pts1 + pts2[1:] + pts3[1:], pth1 + pth2 + pth3
    
    #print len(verts)
    #print len(codes)
    
    #for vert in verts: print vert
    #for code in codes: print code

    assert (len(verts)==len(codes))        
    path = Path(verts, codes)
    patch = mpatches.PathPatch(path, **kwdargs)
    return patch

                  
#rho = np.exp(2*np.pi*1.0j/6.0)
#D = (rho, rho**2, 0.0+0.0j)


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
        return list(map(ffz,zs))
    return ff


def parse_word(ss):
    M = parse_to_mat(ss)
    return mat_to_fcn(M)
    

def invert_letter(letter):
    return {'S':'S', 'T':'U', 'U':'T', 'I':'I'}[letter]


def invert_word(ss):
    return ''.join(invert_letter(s) for s in reversed(ss))


def _inv22(mm):
    a = mm[0,0];  b = mm[0,1]
    c = mm[1,0];  d = mm[1,1]    
    return np.array([[ d,-b],
                     [-c, a]],dtype=int)


def _eq(m1,m2,qq):
    [a,b,c,d] = (m1.dot(_inv22(m2))).ravel()
    if qq>0:
        return ((a%qq==1 and b%qq==0 and c%qq==0 and d%qq==1) or
                ((a+1)%qq==0 and b%qq==0 and c%qq==0 and (d+1)%qq==0))
    else: #qq=0
        return ((a==1 and b==0 and c==0 and d==1) or
                (a==-1 and b==0 and c==0 and d==-1))


def coset_reps_alt(qq, return_generators=False):
    """
    Yield coset representatives for Gamma(q) in SL(2,Z).
    """
        
    I = np.array([[1,0],[0,1]],dtype=int)
    S = np.array([[0,-1],[1,0]],dtype=int)
    T = np.array([[1,1],[0,1]],dtype=int)
    U = _inv22(T)

    gens = []
    retval = [I]

    i = 0
    while i<len(retval):
        mm = retval[i]
        mminv = _inv22(mm)
        for A in [S,T,U]: #[S,S.dot(T).dot(S),S.dot(U).dot(S)]:
            mat1 = mm.dot(A)
            #if mat1[0,0] < 0:
            #    mat1 = -mat1
            #elif mat1[0,0]==0 and mat1[0,1]>0:
            #    mat1 = -mat1
            found = False
            for mat2 in retval:
                if _eq(mat1,mat2,qq):
                    found = True
                    if not _eq(mat1,mat2,0):
                        found2=False
                        #mat3=mminv.dot(A).dot(mm)
                        #mat3=mat2.dot(_inv22(mat1))
                        mat3=mat1.dot(_inv22(mat2))
                        for mat4 in gens:           
                            if _eq(mat3,mat4,0):
                                found2=True
                        if not found2:
                            gens.append(mat3)
                    break
            if not found:
                retval.append(mat1)
        i+=1
    if return_generators:
        return retval, gens
    else:
        return retval


def coset_reps(qq):
    """
    Yield coset representatives for Gamma(q) in SL(2,Z).
    """
    zq2 = itertools.product(range(qq), range(qq))
    #
    # this pattern matching no longer seems to work in python 3.6
    #
    # f=lambda (cc,dd):gcd(cc,gcd(dd,qq))==1
    #
    # so replace it with this:
    #
    def f(pair):
        cc, dd = pair
        return gcd(cc,gcd(dd,qq))==1
        
    for (cc,dd) in filter(f, zq2):
        ii = 0
        if cc==0:
            if dd!=1:
                cc=qq
        while True:
            if gcd(cc,dd-ii*qq) == 1:
                dd -= ii*qq
                break
            elif gcd(cc,dd+ii*qq) == 1:
                dd += ii*qq
                break
            ii += 1
        # cc and dd are now relatively prime.
        zq2 = itertools.product(range(qq), range(qq))
        for (aa,bb) in zq2:
            quot, rem = divmod(aa*dd-bb*cc,qq)
            if rem==1:
                # now aa*dd-bb*cc = 1 (mod qq)
                for (ee,ff) in [(ee,ff) for ee in range(-2*qq,2*qq) for ff in range(-2*qq,2*qq)]:
                    if (aa+ee*qq)*dd-(bb+ff*qq)*cc==1:
                        yieldval = np.array([[aa+ee*qq,bb+ff*qq],[cc,dd]],dtype=int)
                        # now normalize so that the image of D under this
                        # transformation lies between -1/2 and q-1/2.
                        f = mat_to_fcn(yieldval)
                        mm = int(floor(f([0.5j])[0].real))
                        yieldval = np.array([[1,mm%qq-mm],[0,1]],dtype=int).dot(yieldval)
                        yield yieldval
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
    #transfs = [lambda tt: shift_inv(parse_word(tname)(shift(tt))) for tname in tile_names]
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
            ax.add_patch(tile_patch(tile, model='halfplane', facecolor=c))
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


def plot_mat(tile_mats, shifts=[np.array([[1,0],[0,1]],dtype=int)], model='halfplane'):
    """
    Plot a fundamental domain and neighbors for a subgroup G of SL(2,Z).
    """
    rho = np.exp(2*np.pi*1.0j/6.0)
    D = (rho, rho**2, infj)

    vert_xcoords = []
    vert_ycoords = []
    
    tile_groups =[]
    for g in shifts:
        transfs = [mat_to_fcn(g.dot(tile_mat)) for tile_mat in tile_mats]
        tiles = [f(D) for f in transfs]
        for tile in tiles:
            vert_xcoords += [zz.real for zz in tile]
            vert_ycoords += [(1.25 if zz==infj else zz.imag) for zz in tile]
        tile_groups.append(tiles)

    xmax = np.amax(np.array(vert_xcoords))
    xmin = np.amin(np.array(vert_xcoords))
    x_to_y_ratio = (xmax-xmin)/1.5
    
    if model=='halfplane':
        ymax = 1.5
        ymin = 0.0
        x_to_y_ratio = 16.0/3.0
    else:
        boxrad = 1.25
        xmin = -boxrad
        xmax = +boxrad
        ymin = -boxrad
        ymax = +boxrad
        x_to_y_ratio = 1.0
        
    fig = plt.figure(figsize=(x_to_y_ratio*3,3))
    ax = plt.subplot(111)
    plt.subplots_adjust(left=0, bottom=0, right=1,
                        top=1, wspace=0, hspace=0)
    ax.set_xticks([])           
    ax.set_yticks([])
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.grid(True)
    
    colors=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    nc=len(colors)
    
    labeli = 0
    for (groupi,group) in enumerate(tile_groups):
        if groupi==0:
            color = '#704800'
        else:
            color = '#3e2e8e'
        for tile, f in zip(group,transfs):
            ax.add_patch(tile_patch(tile,
                                    model=model,
                                    facecolor=color, 
                                    linewidth=0.2))
            if model=='halfplane':
                zz = f([1.3333j])[0]
            else:
                zz = halfplane_to_poincare_disk(f([1.3333j])[0])
            xx, yy = zz.real, zz.imag
            #plt.text(xx,yy,'%d'%(labeli,))
            labeli += 1
    
    plt.show()


def plot_gamma_2():
    """
    Plot the fundamental region for Gamma(2), the principal
    congruence subgroup of level 2, and neighboring regions
    """
    plot_regions(['I', 'T', 'S', 'TS', 'TSTS', 'TST'],
                  np.exp(2*np.pi*1.0j/6.0),
                  'I',
                  #'SUUS', 
                  ['I','TT','UU','STTS','TTSTTS','TSTTSU','SUUS'])


def plot_gamma_3():
    plot_regions(['I'],0.75j,'I',['I','T','U','S','TS','US','ST','STS','USUUS','TSTS','TSTTS','TST'])


def plot_gammak(k,model='halfplane'):
    """
    Plot the fundamental domain of Gamma(k).
    
    model - 'halfplane' or 'disk'
            'halfplane' - plot it in the upper halfplane model.
            'disk'      - plot it in the poincare disk.
    """
    reps, shifts =coset_reps_alt(k,True)
    plot_mat(reps,[np.array([[1,0],[0,1]],dtype=int)]+shifts,model=model)


def psl2q_order(qq):
    """
    Return the order of PSL(2,Z/qZ).
    """
    ps = list(set(factorize(qq)))
    if qq==2:
        return 6
    else:
        retval = qq**3
        for pp in ps:
            retval *= pp**2-1
            retval /= pp**2
        return retval/2


def plot_arc(a1,a2):
    z1 = np.exp(a1*1.0j*np.pi/180)
    z2 = np.exp(a2*1.0j*np.pi/180)
    
    x1, y1 = c2xy(z1)
    x2, y2 = c2xy(z2)
    
    det = x1*y2-x2*y1
    
    cx = (y2-y1)/det
    cy = (x1-x2)/det
    
    radius = sqrt((x1-cx)**2+(y1-cy)**2)
    
    while a2<a1:
        a2 += 360
    if a2-a1<180:
        theta1 = a1 - 90
        theta2 = a2 + 90
    else:
        theta1 = a1 + 90
        theta2 = a2 - 90

    while theta2 > theta1:
        theta2 -= 360
    
    reverseQ = theta1-theta2>180
    arc_patch=mpatches.Arc((cx,cy),width=2*radius,
                           height=2*radius,fill=False,theta1=theta2,theta2=theta1)
    path = arc_patch.get_path()
    verts = radius*path.vertices + [cx,cy]
    if reverseQ :
        verts = verts[-1::-1,:]
    path2 = Path(verts, path.codes)
    
    fig = plt.figure(figsize=(3,3),facecolor='grey')
    ax = fig.gca()
    
    ax.add_patch(mpatches.Circle((0,0),1,fill=True,alpha=0.5,color='grey'))
    ax.add_patch(mpatches.Circle((cx,cy),radius,fill=False,color='g',alpha=0.1))
    ax.add_patch(mpatches.Arrow(0,0,x1,y1,width=0.01))
    ax.add_patch(mpatches.Arrow(0,0,x2,y2,width=0.01))
    ax.add_patch(mpatches.Arrow(x1,y1,cx-x1,cy-y1,width=0.01))
    ax.add_patch(mpatches.Arrow(x2,y2,cx-x2,cy-y2,width=0.01))
    ax.add_patch(arc_patch)
    ax.add_patch(mpatches.PathPatch(path2))
    
    plt.xlim([-2,+2])
    plt.ylim([-2,+2])
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    
    plt.show()    
    return path2


def plot_complex(zs,**kwargs):
    """
    Plot complex numbers in the complex plane.
    
    zs     - the np.array of complex numbers to be plotted
    kwargs - passed along to plt.plot
    """
    xs = zs.real
    ys = zs.imag
    plt.plot(xs,ys,**kwargs)

    
def random_complex_gaussian():
    """
    Generate a random complex number.
    
    The real and imaginary parts of the generated number
    are iid normal.
    """
    x = np.random.normal()
    y = np.random.normal()
    return x + 1.0j*y

    
class CosetRepsUnitTest(unittest.TestCase):
    def runTest(self):
        for q in range(2,12):
            expected = psl2q_order(q)
            actual1 = len(list(coset_reps(q))) / (1 if q==2 else 2)
            actual2 = len(coset_reps_alt(q))
            self.assertEqual(expected, actual1)
            self.assertEqual(expected, actual2)


class HalfplaneToDiskUnitTest(unittest.TestCase):
    def runTest(self):
        for i in range(100):
            z1 = random_complex_gaussian()
            z1 = z1.real + 1.0j*abs(z1.imag)
            z2 = poincare_disk_to_halfplane(halfplane_to_poincare_disk(z1))
            self.assertAlmostEqual(z1, z2)
        for i in range(100):
            th = np.random.uniform()*2*np.pi
            rr = np.sqrt(np.random.uniform())
            z1 = rr * (np.cos(th) + 1.0j*np.sin(th))
            z2 = halfplane_to_poincare_disk(poincare_disk_to_halfplane(z1))
            self.assertAlmostEqual(z1, z2)
