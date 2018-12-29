import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from matplotlib.path import Path

from math import atan, sin, cos, exp
import numpy as np

infj = complex(0.0, float('inf'))


def c2xy(z):
    return z.real, z.imag


def Sz(z):
    if z == infj:
        return complex(0.0, 0.0)
    else:
        if z == complex(0.0, 0.0):
            return infj
        else:
            return -1.0/z


def S(corners):
    return list(map(Sz, corners))


def Tz(z):
    if z == infj:
        return infj
    else:
        return z + 1.0+0.0j


def T(corners):
    return list(map(Tz, corners))


def Tinvz(z):
    if z == infj:
        return infj
    else:
        return z - 1.0+0.0j


def Tinv(corners):
    return list(map(Tinvz, corners))


def z1z2_to_pts(z1, z2, n):
    x1, y1 = c2xy(z1)
    x2, y2 = c2xy(z2)

    if x1 == x2:
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


if __name__ == '__main__':
    fig = plt.figure(figsize=(20, 2))
    ax = fig.gca()

    rho = exp(2*np.pi*1.0j/6.0)

    D = (rho, rho**2, 0.0+0.0j)
    tiles = [D, T(D), S(D), T(S(D)), T(S(T(S(D)))), T(S(T(D)))]
    for tile in tiles:
        ax.add_patch(tile_patch(tile, **{'facecolor': 'r'}))

    T2 = lambda tt: T(T(tt))
    tiles2 = map(T2, tiles)
    for tile in tiles2 :
        ax.add_patch(tile_patch(tile, **{'facecolor': 'g'}))

    ST2S = lambda tt: S(T(T(S(tt))))
    tiles3 = map(ST2S, tiles)
    for tile in tiles3:
        ax.add_patch(tile_patch(tile, **{'facecolor': 'b'}))

    TST2STinv = lambda tt: T(S(T(T(S(Tinv(tt))))))
    tiles4 = map(TST2STinv, tiles)
    for tile in tiles4:
        ax.add_patch(tile_patch(tile, **{'facecolor': 'y'}))

    tiles5 = map(TST2STinv, tiles4)
    for tile in tiles5:
        ax.add_patch(tile_patch(tile, **{'facecolor': 'm'}))

    tiles6 = map(ST2S, tiles3)
    for tile in tiles6:
        ax.add_patch(tile_patch(tile, **{'facecolor': 'c'}))

    ax.grid()
    ax.set_xlim([-1.0, +4.0])
    ax.set_ylim([0.0, 0.4])

    plt.show()
