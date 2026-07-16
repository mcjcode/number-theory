from itertools import accumulate
from utilities import take, listtake


from fibonacci import (
    fibonacci,
    tribonacci,
    fastfib,
    fastfib_decimate,
    fastfib_sum1,
    fastfib_sum2,
    fastfib_sum3,
    find_recursion,
    find_recursion2,
    cumulative_lrs,
    lrs,
    power_lrs,
    decimate_lrs,
)


def test_fibonacci():
    vals = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
    for x, v in zip(fibonacci(modulus=97), vals):
        assert x==v


def test_tribonacci():
    vals = [0, 0, 1, 1, 2, 4, 7, 13, 24, 44, 81]
    for x, v in zip(tribonacci(modulus=97), vals):
        assert x==v


def test_fastfib():
    N = 10000
    m = 10**9+7

    it = fibonacci(modulus=m)
    for _ in range(N):
        next(it)
    assert next(it)==fastfib(N, modulus=m)

    it = tribonacci(modulus=m)
    for _ in range(N):
        next(it)
    assert next(it)==fastfib(N, modulus=m, coeffs=[1, 1, 1], init=[0, 0, 1])


def test_fastfib_decimate():
    N = 10000
    M = 10**9+7
    m = 10
    i = 3
    assert fastfib_decimate(N, modulus=M)==fastfib(N, modulus=M)
    assert fastfib_decimate(N, modulus=M, m=m, i=i)==fastfib(N*m+i, modulus=M)


def test_fastfib_sum1():
    n = 100
    m = 10**9+7
    g = fibonacci(modulus=m)
    assert sum(take(n, g))%m==fastfib_sum1(n, m=m)


def test_fastfib_sum2():
    n = 100
    m = 10**9+7
    g = fibonacci(modulus=m)
    assert sum(map(lambda x: x*x, take(n, g)))%m==fastfib_sum2(n, m=m)


def test_fastfib_sum3():
    n = 17
    m = 10**9+7
    g = fibonacci(modulus=m)
    assert sum(map(lambda x: x**3, take(n, g)))%m==fastfib_sum3(n, m=m)


def test_find_recursion():
    xs = list(take(20, fibonacci(modulus=10**9)))
    cs = find_recursion(xs, 2)
    assert cs == [1, 1]

    xs = list(take(30, tribonacci(modulus=10**12)))
    cs = find_recursion(xs, 3)
    assert cs == [1, 1, 1]


def test_find_recursion2():
    xs = list(take(20, fibonacci(modulus=10**9)))
    cs = find_recursion2(xs, 2)
    assert cs == [1, 1]

    xs = list(take(30, tribonacci(modulus=10**12)))
    cs = find_recursion2(xs, 3)
    assert cs == [1, 1, 1]


def test_cumulative_lrs():

    # fibonacci case
    cs = [1, 1]
    xs = [0, 1]
    cumcs, cumxs = cumulative_lrs(cs, xs)

    assert cumcs==[2, 0, -1]
    assert cumxs==[0, 1, 2]

    # tribonacci case
    cs = [1, 1, 1]
    xs = [0, 0, 1]
    cumcs, cumxs = cumulative_lrs(cs, xs)

    assert cumcs==[2, 0, 0, -1]
    assert cumxs==[0, 0, 1, 2]

    cs = [1, 0, 1]
    xs = [1, 1, 1]
    cumcs, cumxs = cumulative_lrs(cs, xs)

    xs = list(accumulate(list(take(10, lrs(cs, xs)))))
    ys = list(take(10, lrs(cumcs, cumxs)))

    assert xs==ys


def test_power_lrs():

    # fibonacci square case
    cs = [1, 1]
    sq_cs = power_lrs(cs, 2)
    assert sq_cs==[2, 2, -1]

    N = 10
    p = 10**9+7

    # fibonacci cube case
    # I don't know the lrs but I can compute terms
    cs = [1, 1]
    cb_cs = power_lrs(cs, 3)
    xs = list(x**3 for x in take(N, fibonacci(modulus=p)))
    xs2 = xs[:len(cb_cs)]
    for _ in range(N-len(cb_cs)):
        xcs = zip(xs2[-len(cb_cs):], reversed(cb_cs))
        xs2.append(sum(x*c for (x, c) in xcs))

    assert xs2==xs

    # tribonacci square case
    cs = [1, 1, 1]
    sq_cs = power_lrs(cs, 2)
    xs = list(x**2 for x in take(N, tribonacci(modulus=p)))
    xs2 = xs[:len(sq_cs)]
    for _ in range(N-len(sq_cs)):
        xcs = zip(xs2[-len(sq_cs):], reversed(sq_cs))
        xs2.append(sum(x*c for (x, c) in xcs))

    assert xs2==xs


def test_decimate_lrs():

    xs = [0, 1, 1, 8]
    cs = power_lrs([1, 1], 3)

    xs2 = listtake(100, lrs(cs, xs))[11:84:24]
    cs2 = decimate_lrs(cs, xs, 24)

    cs3, xs3 = cumulative_lrs(cs2, xs2)

    g = lrs(cs3, xs3)

    a, b, i = 0, 1, 0
    nout = 0
    cumsum = 0
    while nout < 20:
        if i%24==11:
            cumsum += a**3
            assert cumsum==next(g)
            nout += 1
        a, b = b, a+b
        i += 1
