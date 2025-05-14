def bps_w_rep(n, ps):
    """
    n: the integer upper bound
    ps: the list of primes

    yields: all integers not exceeding n that can be formed
        from products in ps

    Note: the implementation is iterative, and not recursive,
    and so is suitable for cases where ps is large (i.e. cases
    where the stack would normally grow to the size of ps)
    """
    num_primes = len(ps)
    if n < 1:
        return
    
    yield 1
    #yield []
    
    val = [[ps[-1], num_primes-1, 1, n//ps[-1]]]
    nval = ps[-1]

    while True:
        yield nval
        #yield [(v[0], v[2]) for v in val]
        
        #
        # find the largest prime that is
        # 1) as large or larger than the
        #    primes already in the factorization, and
        # 2) less than what is left
        #
        pi = num_primes-1
        p = ps[pi]
        while p > val[-1][3]:
            if p < val[-1][0]:
                break
            pi -= 1
            if pi < 0:
                break
            p = ps[pi]

        if pi < 0:
            return

        #
        # if there is no such prime, we have to back
        # track over the previous prime used.
        #
        if p < val[-1][0]:  # if there is not an add'l prime
            pi = val[-1][1]
            p = ps[pi-1]
            if len(val) == 1:
                nval //= val[-1][0]**val[-1][2]
                val[-1] = [p, pi-1, 1, n//p]
                nval *= p
            else:
                nval //= val[-1][0]**val[-1][2]
                del val[-1]
                if p == val[-1][0]:
                    val[-1][2] += 1
                    val[-1][3] //= p
                    nval *= p
                else:
                    val.append([p, pi-1, 1, val[-1][3]//p])
                    nval *= p
        elif p == val[-1][0]:  # we can divide by our prime again
            val[-1][2] += 1
            val[-1][3] //= p
            nval *= p
        else:
            val.append([p, pi, 1, val[-1][3]//p])
            nval *= p

