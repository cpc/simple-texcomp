#include <cstdio>

#include "simple_texcomp.hpp"
#include "simple_mathlib.hpp"

using namespace simple;

/* Copy-pasted from simple_astc_int.cpp */
inline uint32_t log2(uint32_t x)
{
    uint32_t res = (x > 0xffff) << 4;
    x >>= res;

    uint32_t shift = (x > 0xff) << 3;
    x >>= shift;
    res |= shift;

    shift = (x > 0xf ) << 2;
    x >>= shift;
    res |= shift;

    shift = (x > 0x3 ) << 1;
    x >>= shift;
    res |= shift;

    res |= (x >> 1);

    return res;
}

/* Copy-pasted from simple_astc_int.cpp */
inline uint32_t approx_inv_u32(uint32_t x)
{
    // Useful constants
    constexpr unsigned int A = 0x0000F0F1;            // 32.0 / 17.0  Q2.15
    constexpr unsigned int B = 0xB4B4B4B5;            // 48.0 / 17.0  Q2.30
    constexpr unsigned int ONE_Q15 = ((1 << 15) - 1); // 1 in Q0.15

    // First, scale the input to be within [0.5, 1.0]
    const int32_t scale = 15 - (int32_t)log2(x);

    const uint32_t shl = (scale < 0) ?      0 : scale;
    const uint32_t shr = (scale < 0) ? -scale :     0;
     // x_sc is 0.5--1.0, so Q0.16 -> Q0.15
    const uint32_t x_sc = (x << shl) >> (shr+1);

    // Then, compute the initial estimate
    // Q4.30 but x_sc is <= 1 so the two MSB are unused => Q2.30
    const uint32_t A_x_sc = A * x_sc;
    // Q2.30 - Q2.30 = Q3.30 -> Q2.30 (won't overflow)
    const uint32_t init = B - A_x_sc;

    // Newthon-Raphson iterations
    // 1st
    uint32_t y0 = init >> 15;             // Q2.15
    uint32_t y00 = init;                  // Q2.30
    uint32_t tmp = (x_sc * y0) >> 15;     // Q0.15 * Q2.15 = Q2.30 -> Q2.15
    // Q2.30 + (Q2.15 * (Q0.15 - Q2.15))
    uint32_t y1 = y00 + y0 * (ONE_Q15 - tmp);

    // 2nd
    y0 = y1 >> 15;
    y00 = y1;
    tmp = (x_sc * y0) >> 15;
    y1 = y00 + y0 * (ONE_Q15 - tmp);

    // 3rd
    y0 = y1 >> 15;
    y00 = y1;// >> 1;
    tmp = (x_sc * y0) >> 15;
    y1 = y00 + y0 * (ONE_Q15 - tmp);

    // The result is scaled down now, we need to scale it back
    y1 >>= 8; // Q10.22
    return (y1 << shl) >> shr;
}

#ifndef NDEBUG

/* Same as mathlib version but with debug prints */
uint32_t approx_inv_u32_debug(uint32_t x)
{
    // First, scale the input to be within [0.5, 1.0]
    int scale = (int)(log2(x));
    scale = 15 - scale;
    printf("  %2d", scale);

    const uint32_t shl = (scale < 0) ?      0 : scale;
    const uint32_t shr = (scale < 0) ? -scale :     0;
    const uint32_t x_sc = (x << shl) >> (shr+1); // x_sc in 0.5--1.0, so Q0.16 -> Q0.15

    printf("  ");
    print_bin_(x_sc, 15, 15);

    // Then, compute the initial estimate
    constexpr uint32_t A = 0x0000F0F1; // 32.0 / 17.0  Q2.15
    constexpr uint32_t B = 0xB4B4B4B5; // 48.0 / 17.0  Q2.30

    const uint32_t A_x_sc = A * x_sc;  // Q4.30 but x_sc is <= 1 so the two MSB are unused => Q2.30
    const uint32_t init = B - A_x_sc;  // Q2.30 - Q2.30 = Q3.30 -> Q2.30 (won't overflow)

    double a = fixed_to_double(A, 15);
    double b = fixed_to_double(B, 30);
    double x_sc_f = fixed_to_double(x_sc, 15);
    double init_f = b - a * x_sc_f;
    double init_fi = fixed_to_double(init, 30);
    printf("  %10.8f  %10.8f  %11.8f", init_f, init_fi, (init_fi - init_f) * 1e6);

    constexpr uint32_t ONE = (1 << 15) - 1;  // 1 in Q0.15

    // Newthon-Raphson iterations
    // 1st
    uint32_t y0 = init >> 15;             // Q2.15
    uint32_t y00 = init;                  // Q2.30
    uint32_t tmp = (x_sc * y0) >> 15;     // Q0.15 * Q2.15 = Q2.30 -> Q2.15
    uint32_t y1 = y00 + y0 * (ONE - tmp); // Q2.30 + (Q2.15 * (Q0.15 - Q2.15))

    double y1_f = init_f + init_f * (1.0 - x_sc_f * init_f);
    double y1_fi = fixed_to_double(y1, 30);
    printf("  %10.8f  %10.8f  %11.8f", y1_f, y1_fi, y1_fi - y1_f);

    // 2nd
    y0 = y1 >> 15;
    y00 = y1;
    tmp = (x_sc * y0) >> 15;
    y1 = y00 + y0 * (ONE - tmp);

    y1_f = y1_f + y1_f * (1.0 - x_sc_f * y1_f);
    y1_fi = fixed_to_double(y1, 30);
    printf("  %10.8f  %10.8f  %11.8f", y1_f, y1_fi, y1_fi - y1_f);

    // 3rd
    y0 = y1 >> 15;
    y00 = y1;// >> 1;
    tmp = (x_sc * y0) >> 15;
    y1 = y00 + y0 * (ONE - tmp);

    y1_f = y1_f + y1_f * (1.0 - x_sc_f * y1_f);
    y1_fi = fixed_to_double(y1, 30);
    printf("  %10.8f  %10.8f  %11.8f", y1_f, y1_fi, y1_fi - y1_f);

    // The result is scaled down now, we need to scale it back
    y1 >>= 8; // Q10.22
    return (y1 << shl) >> shr;
}

void compute_inv(uint32_t x, double xf)
{
    printf("%10.8f  ", xf);
    print_bin_(x, 18, 16);

    uint32_t inv_x = approx_inv_u32_debug(x);
    double gt = 1.0 / xf;

    double inv_x_fi = fixed_to_double(inv_x, 22);
    printf("  %13.8f  %13.8f  %14.8f", gt, inv_x_fi, inv_x_fi - gt);
    printf("\n");
}

#else // not NDEBUG

uint32_t approx_inv_u32_debug(uint32_t x)
{
    printf("Build and run in debug mode\n");
    return 0;
}

void compute_inv(uint32_t x, double xf)
{
    printf("Build and run in debug mode\n");
}

#endif // not NDEBUG

int main()
{
    printf("%10s  %19s  %2s  %16s"
        "  %10s  %10s  %11s  %10s  %10s  %11s"
        "  %10s  %10s  %11s  %10s  %10s  %11s"
        "  %13s  %13s  %14s\n",
        "dp", "dp Q2.16", "sc", "dp_sc",
        "gt init", "init", "dinit 1e-6", "gt #1", "#1", "d#1",
        "gt #2", "#2", "d#2", "gt #3", "#3", "d#3",
        "gt", "1/dp", "err");

    for (uint32_t diff = 0; diff < 256; ++diff)
    {
        if (diff < 8)
        {
            continue;
        }

        double diff_f = (double)(diff) / 256.0;

        double dp_f = diff_f * diff_f;
        uint32_t dp = diff*diff + 0 + 0;

        compute_inv(dp, dp_f);
    }

    for (uint32_t diff = 0; diff < 256; ++diff)
    {
        if (diff < 8)
        {
            continue;
        }

        double diff_f = (double)(diff) / 256.0;

        double dp_f = diff_f*diff_f + diff_f*diff_f + diff_f*diff_f;
        uint32_t dp = diff*diff + diff*diff + diff*diff;

        compute_inv(dp, dp_f);
    }
}
