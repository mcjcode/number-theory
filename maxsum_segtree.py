"""
maxsum_segtree

Implementation of a segment tree with O(log(n)) complexity
for updates, and O(1) retrieval of the maximum sum of a
contiguous subarray.
"""

from dataclasses import dataclass


@dataclass(eq=False) #, match_args=False, slots=True)
class Node:
    """
    left - the maximum sum of contiguous subarray
           that starts at the left end of the array
    right - the maximum sum of a contiguous subarray
           that ends at the right end of the array
    tot - the sum of the entire array
    best - the maximum sum of a subarray of this array
    """
    left: int = 0
    right: int = 0
    tot: int = 0
    best: int = 0

def combine(L: Node, R: Node, out: Node) -> None:
    out.left = max(L.left, L.tot,  L.tot + R.left)
    out.right = max(R.right, R.tot, L.right + R.tot)
    out.tot = L.tot + R.tot
    out.best = max(L.best, R.best, L.right + R.left)


def update(t, pos: int, val) -> None:
    pos += len(t)>>1
    t[pos].left = t[pos].right = t[pos].tot = t[pos].best = val
    
    pos >>= 1
    while pos:
        combine(t[2*pos], t[2*pos+1], t[pos])
        pos >>= 1


def batch_update(t) -> None:
    """
    Sweep through the entire segment tree, from the
    bottom (t[-1]) to the top (t[1]).  This is quicker,
    O(n), than making n full updates, which is O(nlogn),
    where n=len(t)
    """
    pos = len(t)-2
    while pos:
        combine(t[pos], t[pos+1], t[pos//2])
        pos -= 2


def setval(t, pos: int, val):
    """
    Set the element of t at position 'pos' to value val,
    without updating any of the parent elements.
    """
    pos += len(t)>>1
    t[pos].left = t[pos].right = t[pos].tot = t[pos].best = val
