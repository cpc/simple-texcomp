#include <stdint.h>

// Size of the weight grid the input block will be downsampled to
constexpr unsigned int WGT_X = 8;
constexpr unsigned int WGT_Y = 5;
constexpr unsigned int WGT_CNT = WGT_X * WGT_Y;

#ifdef __TCE__

// Size of the pixel blocks the image will be split into
constexpr unsigned int BLOCK_X = 12;
constexpr unsigned int BLOCK_Y = 12;
constexpr unsigned int BLOCK_PX_CNT = BLOCK_X * BLOCK_Y;

inline void _load_4u8(volatile const uint8_t* inp, uchar4* out)
{
    _TCE_LD32(inp, *out);
}

inline void _min_4u8(uchar4 px, uchar4 mincol, uchar4* out)
{
    _TCE_MINU8X4(px, mincol, *out);
}

inline void _max_4u8(uchar4 px, uchar4 maxcol, uchar4* out)
{
    _TCE_MAXU8X4(px, maxcol, *out);
}

inline void _saturating_sub_4u8(uchar4 lhs, uchar4 rhs, uchar4* out)
{
    _TCE_SATSUBU8X4(lhs, rhs, *out);
}

inline void _saturating_dot_acc_4u8(
    uchar4 lhs,
    uchar4 rhs,
    uint32_t acc,
    uint32_t* out
) {
    _TCE_SATACCDOTU8X4(lhs, rhs, acc, *out);
}

inline void _shr_round_4u8(uchar4 inp, uchar4 amt, uchar4* out)
{
    _TCE_SHRRU8X4(inp, amt, *out);
}

#else  // not __TCE__

#include "simple_mathlib.hpp"
#include "simple_texcomp.hpp"

using namespace simple;

typedef Vec4u8 uchar4;

constexpr unsigned int BLOCK_PX_CNT = astc::BLOCK_X * astc::BLOCK_Y;

#ifndef NDEBUG

inline void print_minmax_u8(const char *pre, uchar4 mincol, uchar4 maxcol)
{
    printf("%s mincol: %3d %3d %3d %3d\n",
        pre, mincol.x, mincol.y, mincol.z, mincol.w
    );
    printf("%s maxcol: %3d %3d %3d %3d\n",
        pre, maxcol.x, maxcol.y, maxcol.z, maxcol.w
    );
}

#endif // NDEBUG

inline void _load_4u8(const uint8_t* inp, uchar4* out)
{
    *out = uchar4 { inp[0], inp[1], inp[2], inp[3] };
}

inline void _min_4u8(uchar4 px, uchar4 mincol, uchar4* out)
{
    *out = min4u8(px, mincol);
}

inline void _max_4u8(uchar4 px, uchar4 maxcol, uchar4* out)
{
    *out = max4u8(px, maxcol);
}

inline void _saturating_sub_4u8(uchar4 lhs, uchar4 rhs, uchar4* out)
{
    *out = lhs.satsub(rhs);
}

inline void _shr_round_4u8(uchar4 inp, uchar4 amt, uchar4* out)
{
    // Assuming constant amt vector
    *out = shr_round_4u8(inp, amt.x);
}

inline void _saturating_dot_acc_4u8(
    uchar4 lhs,
    uchar4 rhs,
    uint32_t acc,
    uint32_t* out
) {
    *out = lhs.sataccdot(rhs, acc);
}

// possible short-circuit, but changes the result
// inline uint32_t approx_inv_u32(uint32_t x)
// {
//     constexpr uint32_t ONE = ((1 << 30) - 1); // 1 in Q0.15
//     return ONE / x;
// }

namespace simple::astc {

#endif  // __TCE__

/** Return a pointer to uchar4 vector as a uint32_t pointer */
inline uint32_t* as_u32(const uchar4* x)
{
    return (uint32_t*)x;
}

/** Quantize an 8-bit integer vector into 5 bits (32 values)
 *
 * Returns the quantized value and modifies the input vector in-place.
 */
inline uchar4 quantize_5b(uchar4* vec)
{
    constexpr uchar4 three = { 3, 3, 3, 3 };

    uchar4 quant;
    _shr_round_4u8(*vec, three, &quant);

    // _TCE_SHRRU8X4(*vec, three, quant);

    // here, rounding >>2 seems to hurt quality
    uchar4 dequant = (quant << 3) | (quant >> 2);
    *vec = dequant;

    return quant;
}

/* Compute log2 (i.e., the position of the highest bit set) */
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

 /* Approximate 1/x in a fixed point arithmetic.
  *
  * Assumes x is in a Q2.16 format (unsigned). Returns Q10.22.
  */
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

/** Encode a block of pixels
 *
 * The block_id_x/y could be substituted with get_global_id() if converted to
 * OpenCL.
 */
#ifdef __TCE__

void encode_block_int(
    volatile const uint8_t* __restrict__ inp_img,
    unsigned int block_id_x,
    unsigned int block_id_y,
    unsigned int img_w,
    volatile uint8_t* __restrict__ out_data
){

#else  // __TCE__

void encode_block_int(
    const uint8_t* __restrict__ inp_img,
    unsigned int block_id_x,
    unsigned int block_id_y,
    unsigned int img_w,
    uint8_t* __restrict__ out_data
){

#endif  // __TCE__

    /*************************************************************************/
    /***** Here starts the actual function without the preprocessor mess *****/
    /*************************************************************************/

    unsigned int xb = block_id_x * BLOCK_X;
    unsigned int yb = block_id_y * BLOCK_Y;

    // Input pixel block as 1D array
    uchar4 block_pixels[BLOCK_PX_CNT];

    // Color endpoints - min and max colors in the current block
    uchar4 mincol = { 255, 255, 255, 255 };
    uchar4 maxcol = { 0, 0, 0, 0 };

    // Read pixel block into a 1D array and select min/max at the same time
    for (unsigned int y = 0; y < BLOCK_Y; ++y)
    {
        const unsigned int off = img_w * (yb + y) + xb;

        for (unsigned int x = 0; x < BLOCK_X; ++x)
        {
            const unsigned int i_inp = (off + x) << 2; // assuming NCH_RGB = 4
            const unsigned int i_out = y*BLOCK_X + x;

            uchar4 px;
            _load_4u8(inp_img+i_inp, &px);

            _min_4u8(px, mincol, &mincol);
            _max_4u8(px, maxcol, &maxcol);

            block_pixels[i_out] = px;
        }
    }

    // print_minmax_u8("   ", mincol, maxcol);

    // Move the endpoints line segment such that mincol is at zero
    uchar4 ep_vec;
    _saturating_sub_4u8(maxcol, mincol, &ep_vec);

    // Inset min/max (shrink the bounding box to reduce the effect of outliers)
    const uchar4 inset = ep_vec >> 4;
    mincol = mincol + inset;
    maxcol = maxcol - inset;

    // print_minmax_u8("ins", mincol, maxcol);

    // Quantize the endpoints
    const uchar4 mincol_quant = quantize_5b(&mincol);
    const uchar4 maxcol_quant = quantize_5b(&maxcol);

    // print_minmax_u8("mq ", mincol, maxcol);
    // print_minmax_u8("mqq", mincol_quant, maxcol_quant);
    // print_minmax_u8("ep ", ep_vec, ep_vec);
    // printf("e32: %#010x", *as_u32(&ep_vec));

    // Pixels quantized as 2b weights
    uint8_t quantized_weights[WGT_CNT] = { 0 };

    if (*as_u32(&mincol_quant) != *as_u32(&maxcol_quant))
    {
        // Pixels mapped to the endpoint line without quantization and
        // downsampling
        uint8_t ideal_weights[BLOCK_PX_CNT];

        // Projection of pixels onto ep_vec to get the ideal weights
        // First, we normalize the endpoint vector and scale it.
        uint32_t ep_dot;     // Q2.16
        _saturating_dot_acc_4u8(ep_vec, ep_vec, 0, &ep_dot);

        // Q10.22, max. value 1024 for 5b quant
        const uint32_t inv_ep_dot = approx_inv_u32(ep_dot);

        printf("ep_dot: %#010x\n", ep_dot);
        printf("inv   : %#010x\n", inv_ep_dot);
    }
}

#ifndef __TCE__
} // namespace simple::astc
#endif
