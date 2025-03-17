from collections import deque
import numpy as np


def stern_brocot(continue_condition=bool,
                 start=((0,1), (1,1), (1,0))):
    dqe = deque([start])
    while dqe:
        l, x, r = dqe.popleft()
        yield x
        y = l[0]+x[0], l[1]+x[1] 
        if continue_condition(*y):
            dqe.append((l, y, x))
        y = x[0]+r[0], x[1]+r[1]
        if continue_condition(*y):
            dqe.append((x, y, r))


def pythagorean_triple_tree(cond):
    A=np.array([[ 1,-2, 2],
                [ 2,-1, 2],
                [ 2,-2, 3]])
    B=np.array([[ 1, 2, 2],
                [ 2, 1, 2],
                [ 2, 2, 3]])
    C=np.array([[-1, 2, 2],
                [-2, 1, 2],
                [-2, 2, 3]])
    stack = [(3, 4, 5)]
    while stack:
        abc = stack.pop()
        yield abc
        for m in A, B, C:
            new_abc = m@abc
            if cond(*new_abc):
                stack.append(new_abc)

            
