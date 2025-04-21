class PLFunc1:
    def __init__(self, xs, ys):
        self.xs = xs
        self.ys = ys

    def __call__(self, x):
        a, b = 0, len(self.xs)-1
        while b-a > 1:
            xz = self.xs[(z:=(a+b)//2)]
            if xz < x:
                a = z
            elif xz > x:
                b = z
            else:
                return self.ys[z]
        p = (x-self.xs[a])/(self.xs[b]-self.xs[a])
        return self.ys[a] + p*(self.ys[b]-self.ys[a])
        
    def __add__(self, other):
        if type(other) in [int, float]:
            return type(self)(self.xs, [y+other for y in  self.ys])
        xs = sorted(set(self.xs + other.xs))
        ys = [self(x)+other(x) for x in xs]
        return type(self)(xs, ys)
        
    def __sub__(self, other):
        if type(other) in [int, float]:
            return type(self)(self.xs, [y-other for y in  self.ys])
        xs = sorted(set(self.xs + other.xs))
        ys = [self(x)+other(x) for x in xs]
        return type(self)(xs, ys)

    def __radd__(self, other):
        return self.__add__(other)

    def plot(self, ax, **kwargs):
        if 'linewidth' not in kwargs:
            kwargs['linewidth']=0.75
        ax.plot(self.xs, self.ys, **kwargs)

    def integral(self):
        it = zip(self.xs, self.ys)
        retval = 0
        x0, y0 = next(it)
        for x1, y1 in it:
            retval += (y0+y1)/2 * (x1-x0)
            x0, y0 = x1, y1
        return retval

    def diff(self):
        it = zip(self.xs, self.ys)
        retval = []
        x0, y0 = next(it)
        for x1, y1 in it:
            retval.append((y1-y0)/(x1-x0))
            x0, y0 = x1, y1
        return retval

    def prune(self):
        it = zip(self.xs, self.ys)
        x0, y0 = next(it)
        newxs, newys = [x0], [y0]
        x1, y1 = next(it)
        m0 = (y1-y0)/(x1-x0)
        for x2, y2 in it:
            m1 = (y2-y1)/(x2-x1)
            if m1 != m0:
                newxs.append(x1)
                newys.append(y1)
            x0, y0 = x1, y1
            x1, y1 = x2, y2
            m0 = m1
        newxs.append(x2)
        newys.append(y2)
        return type(self)(newxs, newys)
        
# class PLFunc2:
    
#     def __init__(self, xys):
#         self.xys = xys

#     def __call__(self, x):
#         it = iter(self.xys)
#         x0, y0 = next(it)
#         for x1, y1 in it:
#             if x <= x1:
#                 p = (x-x0)/(x1-x0)
#                 return y0 + p*(y1-y0)
#             x0, y0 = x1, y1

#     def __add__(self, other):
#         if type(other) in [int, float]:
#             return type(self)([(x, y+other) for x, y in self.xys])
#         xs = sorted(set([x for x, _ in self.xys] + [x for x, _ in other.xys]))
#         ys = [self(x)+other(x) for x in xs]
#         return type(self)(list(zip(xs, ys)))

#     def __sub__(self, other):
#         if type(other) in [int, float]:
#             return type(self)([(x, y-other) for x, y in self.xys])
#         xs = sorted(set([x for x, _ in self.xys] + [x for x, _ in other.xys]))
#         ys = [self(x)+other(x) for x in xs]
#         return type(self)(list(zip(xs, ys)))
        
#     def __radd__(self, other):
#         return self.__add__(other)

#     def plot(self, ax, **kwargs):
#         if 'linewidth' not in kwargs:
#             kwargs['linewidth']=0.75
#         ax.plot(*list(zip(*self.xys)), **kwargs)


def PLMax1(plf1, plf2):
    xs = sorted(set(plf1.xs + plf2.xs))
    x0 = xs[0]
    new_xs = [x0]
    y10, y20 = plf1(xs[0]), plf2(xs[0])
    ys = [max(y10, y20)]
    for i in range(1, len(xs)):
        x1 = xs[i]
        y11 = plf1(x1)
        y21 = plf2(x1)
        if y10 >= y20 and y11 >= y21:
            new_xs.append(x1)
            ys.append(y11)
        elif y20 >= y10 and y21 >= y11:
            new_xs.append(x1)
            ys.append(y21)
        elif y10 >= y20 and y11 <= y21:
            t = (y10 - y20)/((y21-y20)-(y11-y10))
            x = x0 + t*(x1-x0)
            new_xs += [x, x1]
            ys += [plf1(x), y21]
        elif y10 <= y20 and y11 >= y21:
            t = (y10 - y20)/((y21-y20)-(y11-y10))
            x = x0 + t*(x1-x0)
            new_xs += [x, x1]
            ys += [plf1(x), y11]
        y10 = y11
        y20 = y21
        x0 = x1
    return type(plf1)(new_xs, ys)


# def PLMax2(plf1, plf2):
#     xs = sorted(set([x for x, _ in plf1.xys] + [x for x, _ in plf2.xys]))
#     new_xs = [xs[0]]
#     y10, y20 = plf1(xs[0]), plf2(xs[0])
#     ys = [max(y10, y20)]
#     for i in range(1, len(xs)):
#         x0, x1 = xs[i-1:i+1]
#         #y10 = plf1(x0)
#         #y20 = plf2(x0)
#         y11 = plf1(x1)
#         y21 = plf2(x1)
#         if y10 >= y20 and y11 >= y21:
#             new_xs.append(x1)
#             ys.append(y11)
#         elif y20 >= y10 and y21 >= y11:
#             ys.append(plf2(x1))
#             new_xs.append(x1)
#         elif y10 >= y20 and y11 <= y21:
#             t = (y10 - y20)/((y21-y20)-(y11-y10))
#             x = x0 + t*(x1-x0)
#             new_xs.append(x)
#             ys.append(plf1(x))
#             new_xs.append(x1)
#             ys.append(y21)
#         elif y10 <= y20 and y11 >= y21:
#             t = (y10 - y20)/((y21-y20)-(y11-y10))
#             x = x0 + t*(x1-x0)
#             new_xs.append(x)
#             ys.append(plf1(x))
#             new_xs.append(x1)
#             ys.append(y11)
#         y10 = y11
#         y20 = y21

#     return type(plf1)(list(zip(new_xs, ys)))

PLFunc = PLFunc1
PLMax = PLMax1