#!/usr/bin/env python3

from collections.abc import Iterable

import pandas as pd
import numpy as np

# https://github.com/ZZZZzzzzac/numfi
from numfi import numfi

from matplotlib import pyplot as plt

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
    return res, sc, x_sc


def limits():
    mins = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 0.0]),
    ]

    maxs = [
        np.array([1 / 32, 0.0, 0.0]),
        np.array([1.0, 1.0, 1.0]),
    ]

    pxs = [
        np.array([1.0, 1.0, 1.0]),
    ]

    dcols = [
        'ep[0]',
        'dp',
        'inv_dp',
        'ep_sc[0]',
        'ep_sc2[0]',
        'diff_dp',
        'diff_dp2',
    ]

    print('Precise')
    d = pd.DataFrame(columns=dcols)

    for mn, mx in zip(mins, maxs):
        ep = mx - mn
        dp = ep.dot(ep)
        inv_dp = 1 / dp
        ep_sc = ep * inv_dp
        ep_sc2 = (ep * inv_dp) / 32

        for px in pxs:
            diff = px - mn
            diff_dp = diff.dot(ep_sc)
            diff_dp2 = diff.dot(ep_sc2)
            dnew = pd.DataFrame(
                [[ ep[0], dp, inv_dp, ep_sc[0], ep_sc2[0], diff_dp, diff_dp2 ]],
                columns=dcols,
            )
            d = d.append(dnew, ignore_index=True)
    print(d)

    print('\nFixed point')
    mins = [ numfi(mn, s=0, w=8, f=8, rounding='floor', overflow='saturate', fixed=True) for mn in mins ]
    maxs = [ numfi(mx, s=0, w=8, f=8, rounding='floor', overflow='saturate', fixed=True) for mx in maxs ]
    d = pd.DataFrame(columns=dcols)
    print('Min:', mins)
    print('Max:', maxs)
    print('{:>19s} {:>19s} {:>19s}  {:>19s}  {:>19s}'.format('ep[0]', 'dp', 'inv_dp', 'ep_sc[0]', 'ep_sc2[0]'))

    for mn, mx in zip(mins, maxs):
        ep = (mx - mn).requantize(w=16, overflow='wrap')
        dp = numfi(ep.dot(ep), s=0, w=21, f=10, rounding='floor', overflow='wrap', fixed=True)
        inv_dp = 1 / dp
        ep_sc = ep * inv_dp
        ep_sc2 = (ep * inv_dp) >> 5

        print('{:>19s}  {:>19s}  {:>19s}  {:>19s}  {:>19s}'.format(ep.bin_[0], dp.bin_[0], inv_dp.bin_[0], ep_sc.bin_[0], ep_sc2.bin_[0]))

        _ep = ep[0].precision * ep[0].int
        _ep_sc = ep_sc[0].precision * ep_sc[0].int
        _ep_sc2 = ep_sc2[0].precision * ep_sc2[0].int
        _dp = dp[0].precision * dp[0].int
        _inv_dp = inv_dp[0].precision * inv_dp[0].int

        for px in pxs:
            diff = px - mn
            diff_dp = diff.dot(ep_sc)
            diff_dp2 = diff.dot(ep_sc2)
            dnew = pd.DataFrame(
                [[ _ep[0], _dp, _inv_dp, _ep_sc[0], _ep_sc2[0], diff_dp, diff_dp2 ]],
                columns=dcols,
            )
            d = d.append(dnew, ignore_index=True)

    print()
    print(d)

    print()
    s = np.linspace(0, 31, 32) / 32
    r = np.zeros([32*32*32, 3])
    i = 0
    for x in s:
        for y in s:
            for z in s:
                r[i] = np.array([x, y, z])
                i += 1
    mask = np.logical_not(r.any(axis=1))  # delete zeroes
    r = np.delete(r, mask, axis=0)
    dp = np.square(r).sum(axis=1)
    inv_dp = 1/dp

    sorted_inv_dp = np.sort(inv_dp)
    nbins = 2
    chunks = np.array_split(sorted_inv_dp, nbins)
    for chunk in chunks:
        print('chunk: {} -- {} ({})'.format(np.min(chunk), np.max(chunk), len(chunk)))

    lt = inv_dp[inv_dp <= 1.0]
    gt = inv_dp[inv_dp > 1.0]

    print('less than one: {} -- {} ({})'.format(np.min(lt), np.max(lt), len(lt)))
    print('more than one: {} -- {} ({})'.format(np.min(gt), np.max(gt), len(gt)))

    plt.scatter(inv_dp, np.zeros_like(inv_dp), s=1)
    plt.scatter(dp, np.ones_like(inv_dp), s=1)
    plt.show()


if __name__ == '__main__':
    pd.set_option("display.precision", 10)

    n = 256 * 3
    step = 16
    series = np.arange(1, n, step)

    x = np.arange(1.0/(32.0*32.0), 3.0, 1.0/16.0)
    gt = 1.0 / x

    n_iter = 2
    print('iterations: {}'.format(n_iter))
    approx, sc = approx_newton_float(x, n_iter)

    xfi = fixed(x)
    approxfi, scfi, x_sc = approx_newton_fixed(xfi, n_iter)

    df = pd.DataFrame(data=np.array([x, sc, scfi, xfi, xfi.bin_, x_sc, gt, approx, approxfi, fixed(approxfi).int]).T,
        columns=['x', 'sc', 'scfi', 'xfi', 'xfib', 'x_sc', 'gt', 'flt_newt', 'fi_newt', 'fi_newt_int'])
    print(df)

    A = 32.0 / 17.0
    B = 48.0 / 17.0
    nums = [
        0.4039370120,
        0.4039370120 * 2.0,
        A,
        B,
        B - A * 0.4039370120 * 2.0,
    ]
    fi = fixed(nums, s=0, w=17, f=15)
    print(nums)
    print(fi)
    print(fi.bin)
    print(fi.bin_)
    print(fi.int)
    print(fi.hex)
    print(fi.precision)

    print()
    limits()
