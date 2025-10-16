# -*- coding: utf-8 -*-


def assert_equal(expected, actual, msg):
    if expected == actual:
        print('.', end=' ')
        return True
    else:
        raise Exception('testing "%s": %s expected, %s actual' % (msg, expected, actual))


def assert_exception(f, msg):
    """
    Test whether the code in the function f
    throws an exception
    """
    try:
        f()
    except Exception:
        print('.', end=' ')
        return
    raise Exception(msg)
