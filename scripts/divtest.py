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
    print(A, B)
    init = B - A * x_sc

    # Newthon-Raphson iterations
    y0 = init
    y1 = init
    for _ in range(n):
        y1 = 2.0 * y0 - x_sc * y0 * y0
        y0 = y1

    return y1, sc


def fixed(x):
    if not isinstance(x, Iterable):
        x = [x]
    return numfi(x, s=0, w=8, f=8, rounding='round', overflow='saturate', fixed=True)


def approx_newton_fixed(x, n):
    #  if x == 0.0:
    #      return 0.0

    # Scale x to be within [0.5, 1.0]
    x_sc = fixed(np.copy(x))
    sc = fixed(np.ones_like(x))
    for i, xi in enumerate(x_sc):
        while x_sc[i] < 0.5:
            x_sc[i] = x_sc[i] << 2
            sc[i] = sc[i] << 2

        while x_sc[i] > 1.0:
            x_sc[i] = x_sc[i] >> 2
            sc[i] = sc[i] >> 2

    # Initial estimate
    A = fixed(32.0 / 17.0)
    B = fixed(48.0 / 17.0)
    print(A, B)
    init = B - A * x_sc

    # Newthon-Raphson iterations
    y0 = init
    y1 = init
    for _ in range(n):
        y1 = fixed(2.0) * y0 - x_sc * y0 * y0
        y0 = y1

    return y1, sc


if __name__ == '__main__':
    pd.set_option("display.precision", 10)

    n = 256
    step = 8
    series = np.arange(1, n, step)

    x = series / n

    approx, sc = approx_newton_float(x, 1)

    gt = 1.0 / x
    gt /= sc

    xfi = numfi(x, s=0, w=8, f=8, rounding='round', overflow='saturate', fixed=True)
    approxfi, scfi = approx_newton_fixed(x, 1)

    df = pd.DataFrame(data=np.array([x, sc, xfi, gt, approx, approxfi]).T, columns=['x', 'sc', 'xfi', 'gt', 'flt_newt', 'fi_newt'])
    print(df)

    #  print("x,gt,approx,diff")
    #  for i in range(0, n+1, step):
    #      x = i / n * 3
    #      x_orig = x
    #      if x == 0:
    #          continue

    #      approx, sc = approx_newton_float(x, 2)

    #      gt = 1.0 / x
    #      gt /= sc

    #      print("{:.20f},{:.20f},{:.20f},{:.20f}".format(x_orig, gt, approx, gt-approx))

    num = 1.0
    fi = numfi([num], s=0, w=8, f=8, rounding='round', overflow='saturate', fixed=True)
    print(num, fi, fi.bin, fi.bin_, fi.hex, fi.precision)
