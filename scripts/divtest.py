#!/usr/bin/env python3

import pandas as pd

# https://github.com/ZZZZzzzzac/numfi
from numfi import numfi


def approx_newton_float(x, n):
    if x == 0.0:
        return 0.0

    # Scale x to be within [0.5, 1.0]
    x_sc = x
    sc = 1
    while x_sc < 0.5:
        x_sc = x_sc * 2.0
        sc = sc * 2.0

    while x_sc > 1.0:
        x_sc = x_sc / 2.0
        sc = sc / 2.0

    # Initial estimate
    init = 48.0 / 17.0 - (32.0 / 17.0) * x_sc

    # Newthon-Raphson iterations
    y0 = init
    y1 = init
    for _ in range(n):
        y1 = 2.0 * y0 - x_sc * y0 * y0
        y0 = y1

    return y1, sc


if __name__ == "__main__":
    n = 256
    step = 8
    print("x,gt,approx,diff")
    for i in range(0, n+1, step):
        x = i / n * 3
        x_orig = x
        if x == 0:
            continue

        approx, sc = approx_newton_float(x, 2)

        gt = 1.0 / x
        gt /= sc

        print("{:.20f},{:.20f},{:.20f},{:.20f}".format(x_orig, gt, approx, gt-approx))

    num = 0.68751
    fi = numfi([num], s=0, w=8, f=8, rounding='round', overflow='saturate', fixed=True)
    print(num, fi, fi.bin, fi.bin_, fi.hex, fi.precision)
