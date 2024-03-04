from collections import deque

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

