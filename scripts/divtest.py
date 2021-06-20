#!/usr/bin/env python3

from collections.abc import Iterable

import pandas as pd
import numpy as np

# https://github.com/ZZZZzzzzac/numfi
from numfi import numfi


def approx_newton_float(x, n):
    #  if x == 0.0:
    #      return 0.0

    # Scale x to be within [0.5, 1.0]
    x_sc = np.copy(x)
    sc = np.ones_like(x)
    for i, xi in enumerate(x_sc):
        while x_sc[i] < 0.5:
            x_sc[i] = x_sc[i] * 2.0
            sc[i] = sc[i] * 2.0

        while x_sc[i] > 1.0:
            x_sc[i] = x_sc[i] / 2.0
            sc[i] = sc[i] / 2.0

    # Initial estimate
    B = 48.0 / 17.0
    A = 32.0 / 17.0
    init = B - A * x_sc

    # Newthon-Raphson iterations
    y0 = init
    y1 = init
    for _ in range(n):
        y1 = 2.0 * y0 - x_sc * y0 * y0
        y0 = y1

    return y1*sc, sc


def fixed(x, s=0, w=16, f=8, rounding='round', overflow='saturate', fixed=True):
    if not isinstance(x, Iterable):
        x = [x]
    return numfi(x, s=s, w=w, f=f, rounding=rounding, overflow=overflow, fixed=fixed)


def scale(x):
    res = np.zeros(x.shape, dtype=np.int8)

    for i, xx in enumerate(x):
        xx = xx.int

        # Essentially, log2 computation, borrowed from:
        # https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
        log2x = (xx > 0xFF) << 3
        xx >>= log2x

        shift = (xx > 0xF ) << 2
        xx >>= shift
        log2x |= shift;

        shift = (xx > 0x3 ) << 1
        xx >>= shift
        log2x |= shift

        log2x |= (xx >> 1)

        res[i] = 7 - log2x[0]

    return res


def approx_newton_fixed(x, n):
    # Scale x to be within [0.5, 1.0]
    x_sc = fixed(np.copy(x))
    sc = scale(x_sc)
    x_sc = np.where(sc >= 0, x_sc << sc, x_sc >> -sc)

    # Initial estimate
    A = fixed(32.0 / 17.0)
    B = fixed(48.0 / 17.0)
    init = B - A * x_sc

    # Newthon-Raphson iterations
    y0 = init
    y1 = init
    for _ in range(n):
        y1 = fixed(2.0) * y0 - x_sc * y0 * y0
        y0 = y1

    # Scale result back
    res = np.where(sc >= 0, y1 << sc, y1 >> -sc)
    return res, sc


if __name__ == '__main__':
    pd.set_option("display.precision", 10)

    n = 256 * 3
    step = 16
    series = np.arange(1, n, step)

    x = series / 256
    gt = 1.0 / x

    n_iter = 2
    print('iterations: {}'.format(n_iter))
    approx, sc = approx_newton_float(x, n_iter)

    xfi = fixed(x)
    approxfi, scfi = approx_newton_fixed(xfi, n_iter)

    df = pd.DataFrame(data=np.array([x, sc, scfi, xfi, xfi.bin_, gt, approx, approxfi]).T,
        columns=['x', 'sc', 'scfi', 'xfi', 'xfib', 'gt', 'flt_newt', 'fi_newt'])
    print(df)

    nums = [
        8.0 / 255.0 / 16.0,
        32.0 / 17.0,
        48.0 / 17.0,
    ]
    fi = fixed(nums, w=16, f=8)
    print(nums)
    print(fi)
    print(fi.bin)
    print(fi.bin_)
    print(fi.hex)
    print(fi.precision)
