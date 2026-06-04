def digits(n, base=10):
    """
    :param n: an integer
    :return: a list ds so that the base^i digit of n is ds[i]
    """
    while n:
        n, d = divmod(n, base)
        yield d


def num_from_digits(ds, base=10):
    """
    :param ds: a list of integers from [0,base)
    :return: the integer whose base^i digit is ds[i]
    """
    n = 0
    a = 1
    for d in ds:
        n += a*d
        a *= base
    return n
