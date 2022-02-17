#include <cstdint>

// #include <stdlib.h>

#include "platform.hpp"


// Force-support vdotq_u32 intrinsic (only works with -march=armv8.4-a
#define __ARM_FEATURE_DOTPROD
#include <cstring>
#include <arm_neon.h>

#include "simple_mathlib.hpp"
#include "simple_texcomp.hpp"

#include <Tracy.hpp>

namespace simple::astc {

typedef Vec4u8 uchar4;

constexpr unsigned int BLOCK_PX_CNT = BLOCK_X * BLOCK_Y;

// Size of the weight grid the input block will be downsampled to
constexpr unsigned int WGT_X = 8;
constexpr unsigned int WGT_Y = 5;
constexpr unsigned int WGT_CNT = WGT_X * WGT_Y;

// Encoded block size (in bytes)
constexpr unsigned int BLOCK_SIZE_EXP = 4;
constexpr unsigned int BLOCK_SIZE = (1 << BLOCK_SIZE_EXP);

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

inline void _load_4u8(volatile const uint8_t* inp, uchar4* out)
{
    *out = uchar4 { inp[0], inp[1], inp[2], inp[3] };
}

inline void _load_u32(uint8_t* inp, uint32_t* out)
{
    *out = *reinterpret_cast<uint32_t*>(inp);
}

inline void _store_u32(uint8_t* out, uint32_t val)
{
    // *out = val;
    out[0] = (uint8_t)(val & 0xff);
    out[1] = (uint8_t)((val >> 8) & 0xff);
    out[2] = (uint8_t)((val >> 16) & 0xff);
    out[3] = (uint8_t)((val >> 24) & 0xff);
}

inline void _min_4u8(uchar4 px, uchar4 mincol, uchar4* out)
{
    *out = min4u8(px, mincol);
}

inline void _max_4u8(uchar4 px, uchar4 maxcol, uchar4* out)
{
    *out = max4u8(px, maxcol);
}

inline void _max_u32(uint32_t lhs, uint32_t rhs, uint32_t* out)
{
    *out = u32max(lhs, rhs);
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

inline void _shr_round_u32(uint32_t inp, uint8_t amt, uint32_t* out)
{
    // Assuming constant amt vector
    *out = shr_round_u32(inp, amt);
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

/* Unpack the first three elements of uchar4 vector into 32b array */
inline void _unpack_rgb_u32(uchar4 v, uint32_t out[3])
{
    out[0] = static_cast<uint32_t>(v.x);
    out[1] = static_cast<uint32_t>(v.y);
    out[2] = static_cast<uint32_t>(v.z);
}

/* Unpack the first three elements of uchar4 vector into 32b array */
inline void _unpack_rgb_u8(uchar4 v, uint8_t out[3])
{
    out[0] = v.x;
    out[1] = v.y;
    out[2] = v.z;
}

inline void _reflect_u32(uint32_t inp, uint32_t* out)
{
    const uint8_t rev_bytes[4] = {
        (uint8_t)(astc::bitrev8((uint8_t)((inp >> 24) & 0xff))),
        (uint8_t)(astc::bitrev8((uint8_t)((inp >> 16) & 0xff))),
        (uint8_t)(astc::bitrev8((uint8_t)((inp >> 8) & 0xff))),
        (uint8_t)(astc::bitrev8((uint8_t)((inp >> 0) & 0xff))),
    };

    *out = *reinterpret_cast<const uint32_t*>(rev_bytes);
}

inline void _write_bits(int value, int bitcount, int bitoffset, uint8_t* ptr)
{
    astc::write_bits(value, bitcount, bitoffset, ptr);
}

// Precomputed coefficients used in horizontal bilinear interpolation
constexpr uchar4 BILIN_WEIGHTS_X_8_U8[10] = {          //  IDX:
    { 187,  68,   0,   0, },                           //  {  0,  1,  0, },
    {   0, 111, 128,  16, },                           //  {  1,  2,  3, },
    {   0,   0,  42, 142, }, {  71,   0,   0,   0, },  //  {  2,  3,  4, },
    {  90, 135,  30,   0, },                           //  {  4,  5,  6, },
    {   0,  30, 135,  90, },                           //  {  5,  6,  7, },
    {   0,   0,   0,  71, }, { 142,  42,   0,   0, },  //  {  7,  8,  9, },
    {  16, 128, 111,   0, },                           //  {  8,  9, 10, },
    {   0,   0,  68, 187, },                           //  { 10, 11,  0, },
};

// Precomputed coefficients used in vertical bilinear interpolation
// IDX:
// {  0,  1,  2,  0,  0,  0, },
// {  1,  2,  3,  4,  5,  0, },
// {  3,  4,  5,  6,  7,  8, },
// {  6,  7,  8,  9, 10,  0, },
// {  9, 10, 11,  0,  0,  0, },
constexpr uchar4 BILIN_WEIGHTS_Y_5_U8[9] = {
    { 134,  85,  36,   0, },
    {   0,  34,  68,  85, }, {  51,  17,   0,   0, },
    // bumped two middle values to sum up to 254:
    {   0,   0,   0,   8, }, {  42,  77,  77,  42, }, {   8,   0,   0,   0, },
    {   0,   0,  17,  51, }, {  85,  68,  34,   0, },
    {   0,  36,  85, 134, },
};


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
    // y0 = y1 >> 15;
    // y00 = y1;
    // tmp = (x_sc * y0) >> 15;
    // y1 = y00 + y0 * (ONE_Q15 - tmp);

    // 3rd
    // y0 = y1 >> 15;
    // y00 = y1;// >> 1;
    // tmp = (x_sc * y0) >> 15;
    // y1 = y00 + y0 * (ONE_Q15 - tmp);

    // printf("x: %d - scale: %d - shl: %d - shr %d - x_sc %d - A_x_sc %d - init %d - y1 %d - res: %u\n",
    //     x, scale, shl, shr, x_sc, A_x_sc, init, y1, ((y1 << shl) >> shr));

    // The result is scaled down now, we need to scale it back
    y1 >>= 8; // Q10.22
    // TODO: The << *will* overflow (e.g., for x = 50)
    return ((y1 << shl) >> shr) & 0xffffffff;
}

/* Pack the first three elements of input 32b array into uchar4 vector */
inline void pack_rgb_u8(uint32_t a[3], uchar4* out)
{
    out->x = (uint8_t)(a[0]);
    out->y = (uint8_t)(a[1]);
    out->z = (uint8_t)(a[2]);
}

inline void read_block_min_max(
    const uint8_t* __restrict__ inp,
    unsigned int block_id_x,
    unsigned int block_id_y,
    unsigned int img_w,
    uchar4* mincol,
    uchar4* maxcol,
    uint8x16_t block_pixels_x4[36]  // 36 == 12 * 3 == 9 * 4
) {
    ZoneScopedN("minmax");

    unsigned int xb = block_id_x * BLOCK_X;
    unsigned int yb = block_id_y * BLOCK_Y;

    // Color endpoints - min and max colors in the current block
    uint8x16_t min_x4 = {
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255
    };
    uint8x16_t max_x4 = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    for (unsigned int y = 0; y < BLOCK_Y; ++y)
    {
        const unsigned int off = img_w * (yb + y) + xb;

        uint8x16_t r = vld1q_u8(inp + NCH_RGB * (off + 0));  // pixels 0..3
        uint8x16_t g = vld1q_u8(inp + NCH_RGB * (off + 4));  // pixels 4..7
        uint8x16_t b = vld1q_u8(inp + NCH_RGB * (off + 8));  // pixels 8..11

        min_x4 = vminq_u8(min_x4, r);
        min_x4 = vminq_u8(min_x4, g);
        min_x4 = vminq_u8(min_x4, b);

        max_x4 = vmaxq_u8(max_x4, r);
        max_x4 = vmaxq_u8(max_x4, g);
        max_x4 = vmaxq_u8(max_x4, b);

        block_pixels_x4[3*y+0] = r;
        block_pixels_x4[3*y+1] = g;
        block_pixels_x4[3*y+2] = b;
    }

    const uint8x16_t mask_min_r_x4 = {
        0,   255, 255, 255, 0,   255, 255, 255,
        0,   255, 255, 255, 0,   255, 255, 255
    };
    const uint8x16_t mask_min_g_x4 = {
        255, 0,   255, 255, 255, 0,   255, 255,
        255, 0,   255, 255, 255, 0,   255, 255
    };
    const uint8x16_t mask_min_b_x4 = {
        255, 255, 0,   255, 255, 255, 0,   255,
        255, 255, 0,   255, 255, 255, 0,   255
    };

    const uint8x16_t mask_max_r_x4 = {
        255, 0,   0,   0,   255, 0,   0,   0,
        255, 0,   0,   0,   255, 0,   0,   0
    };
    const uint8x16_t mask_max_g_x4 = {
        0,   255, 0,   0,   0,   255, 0,   0,
        0,   255, 0,   0,   0,   255, 0,   0
    };
    const uint8x16_t mask_max_b_x4 = {
        0,   0,   255, 0,   0,   0,   255, 0,
        0,   0,   255, 0,   0,   0,   255, 0
    };

    const uint8x16_t min_r_x4 = vorrq_u8(min_x4, mask_min_r_x4);
    const uint8x16_t min_g_x4 = vorrq_u8(min_x4, mask_min_g_x4);
    const uint8x16_t min_b_x4 = vorrq_u8(min_x4, mask_min_b_x4);

    const uint8x16_t max_r_x4 = vandq_u8(max_x4, mask_max_r_x4);
    const uint8x16_t max_g_x4 = vandq_u8(max_x4, mask_max_g_x4);
    const uint8x16_t max_b_x4 = vandq_u8(max_x4, mask_max_b_x4);

    const uint8_t min_r = vminvq_u8(min_r_x4);  // requires A64
    const uint8_t min_g = vminvq_u8(min_g_x4);  // requires A64
    const uint8_t min_b = vminvq_u8(min_b_x4);  // requires A64

    const uint8_t max_r = vmaxvq_u8(max_r_x4);  // requires A64
    const uint8_t max_g = vmaxvq_u8(max_g_x4);  // requires A64
    const uint8_t max_b = vmaxvq_u8(max_b_x4);  // requires A64

    *mincol = uchar4 { min_r, min_g, min_b, 255 };
    *maxcol = uchar4 { max_r, max_g, max_b, 255 };
}

inline void downsample_12x12_to_8x5_u8_quant(
    const uint8_t inp[BLOCK_PX_CNT + (16 - BLOCK_X)],
    uint8_t out[WGT_CNT]
){
    ZoneScopedN("bilin");

    constexpr unsigned int w_inp = 12;
    constexpr unsigned int h_inp = 12;
    constexpr unsigned int w_out = 8;
    constexpr unsigned int h_out = 5;

    constexpr uint8x8_t BILIN_WEIGHTS_X_0 = {
        187, 111,  42,  90,  30,  71,  16,  68,
    };

    constexpr uint8x8_t BILIN_WEIGHTS_X_1 = {
         68, 128, 142, 135, 135, 142, 128, 187,
    };

    constexpr uint8x8_t BILIN_WEIGHTS_X_2 = {
          0,  16,  71,  30,  90,  42, 111,  0,
    };

    constexpr uint8x8_t IDX_X = { 0, 1, 2, 4, 5, 7, 8, 10 };

    constexpr uint8x8_t BILIN_WEIGHTS_Y_0[3] = {
        { 134, 134, 134, 134, 134, 134, 134, 134, },
        { 85,  85,  85,  85,  85,  85,  85,  85,  },
        { 36,  36,  36,  36,  36,  36,  36,  36,  },
    };

    constexpr uint8x8_t BILIN_WEIGHTS_Y_1[5] = {
        { 34, 34, 34, 34, 34, 34, 34, 34, },
        { 68, 68, 68, 68, 68, 68, 68, 68, },
        { 85, 85, 85, 85, 85, 85, 85, 85, },
        { 51, 51, 51, 51, 51, 51, 51, 51, },
        { 17, 17, 17, 17, 17, 17, 17, 17, },
    };

    constexpr uint8x8_t BILIN_WEIGHTS_Y_2[6] = {
        { 8,  8,  8,  8,  8,  8,  8,  8,  },
        { 42, 42, 42, 42, 42, 42, 42, 42, },
        { 77, 77, 77, 77, 77, 77, 77, 77, },
        { 77, 77, 77, 77, 77, 77, 77, 77, },
        { 42, 42, 42, 42, 42, 42, 42, 42, },
        { 8,  8,  8,  8,  8,  8,  8,  8,  },
    };

    constexpr uint8x8_t BILIN_WEIGHTS_Y_3[5] = {
        { 17, 17, 17, 17, 17, 17, 17, 17, },
        { 51, 51, 51, 51, 51, 51, 51, 51, },
        { 85, 85, 85, 85, 85, 85, 85, 85, },
        { 68, 68, 68, 68, 68, 68, 68, 68, },
        { 34, 34, 34, 34, 34, 34, 34, 34, },
    };

    constexpr uint8x8_t BILIN_WEIGHTS_Y_4[3] = {
        { 36,  36,  36,  36,  36,  36,  36,  36,  },
        { 85,  85,  85,  85,  85,  85,  85,  85,  },
        { 134, 134, 134, 134, 134, 134, 134, 134, },
    };

    uint8x8_t tmp[h_inp];

    // First, interpolate rows.
    for (unsigned int y = 0; y < h_inp; ++y)
    {
        // Align the input values so that 1st, 2nd and 3rd samples of the input
        // row are at the 0th position of row0/1/2
        const uint8x16_t row0_u16 = vld1q_u8(inp + y*w_inp);
        const uint8x16_t row1_u16 = vextq_u8(row0_u16, row0_u16, 1);
        const uint8x16_t row2_u16 = vextq_u8(row0_u16, row0_u16, 2);

        // We calculate only 8 values => pack them into one 8-byte vector
        const uint8x8_t row0 = vqtbl1_u8(row0_u16, IDX_X);  // table select
        const uint8x8_t row1 = vqtbl1_u8(row1_u16, IDX_X);
        const uint8x8_t row2 = vqtbl1_u8(row2_u16, IDX_X);

        // Results storage
        uint16x8_t res_u16 = { 0 };

        // Calculate the dot product by two multiply-adds
        res_u16 = vmlal_u8(res_u16, row0, BILIN_WEIGHTS_X_0);
        res_u16 = vmlal_u8(res_u16, row1, BILIN_WEIGHTS_X_1);
        res_u16 = vmlal_u8(res_u16, row2, BILIN_WEIGHTS_X_2);

        // Shift back from 16-bit to 8-bit precision and store
        res_u16 = vshrq_n_u16(res_u16, 8);
        tmp[y] = vmovn_u16(res_u16);
    }

    // Next, columns (unrolled)
    uint16x8_t res_u16[h_out] = {
        { 0, 0, 0, 0, 0, 0, 0, 0, },
        { 0, 0, 0, 0, 0, 0, 0, 0, },
        { 0, 0, 0, 0, 0, 0, 0, 0, },
        { 0, 0, 0, 0, 0, 0, 0, 0, },
        { 0, 0, 0, 0, 0, 0, 0, 0, },
    };

    res_u16[0] = vmlal_u8(res_u16[0], tmp[0], BILIN_WEIGHTS_Y_0[0]);
    res_u16[0] = vmlal_u8(res_u16[0], tmp[1], BILIN_WEIGHTS_Y_0[1]);
    res_u16[0] = vmlal_u8(res_u16[0], tmp[2], BILIN_WEIGHTS_Y_0[2]);

    res_u16[1] = vmlal_u8(res_u16[1], tmp[1], BILIN_WEIGHTS_Y_1[0]);
    res_u16[1] = vmlal_u8(res_u16[1], tmp[2], BILIN_WEIGHTS_Y_1[1]);
    res_u16[1] = vmlal_u8(res_u16[1], tmp[3], BILIN_WEIGHTS_Y_1[2]);
    res_u16[1] = vmlal_u8(res_u16[1], tmp[4], BILIN_WEIGHTS_Y_1[3]);
    res_u16[1] = vmlal_u8(res_u16[1], tmp[5], BILIN_WEIGHTS_Y_1[4]);

    res_u16[2] = vmlal_u8(res_u16[2], tmp[3], BILIN_WEIGHTS_Y_2[0]);
    res_u16[2] = vmlal_u8(res_u16[2], tmp[4], BILIN_WEIGHTS_Y_2[1]);
    res_u16[2] = vmlal_u8(res_u16[2], tmp[5], BILIN_WEIGHTS_Y_2[2]);
    res_u16[2] = vmlal_u8(res_u16[2], tmp[6], BILIN_WEIGHTS_Y_2[3]);
    res_u16[2] = vmlal_u8(res_u16[2], tmp[7], BILIN_WEIGHTS_Y_2[4]);
    res_u16[2] = vmlal_u8(res_u16[2], tmp[8], BILIN_WEIGHTS_Y_2[5]);

    res_u16[3] = vmlal_u8(res_u16[3], tmp[6],  BILIN_WEIGHTS_Y_3[0]);
    res_u16[3] = vmlal_u8(res_u16[3], tmp[7],  BILIN_WEIGHTS_Y_3[1]);
    res_u16[3] = vmlal_u8(res_u16[3], tmp[8],  BILIN_WEIGHTS_Y_3[2]);
    res_u16[3] = vmlal_u8(res_u16[3], tmp[9],  BILIN_WEIGHTS_Y_3[3]);
    res_u16[3] = vmlal_u8(res_u16[3], tmp[10], BILIN_WEIGHTS_Y_3[4]);

    res_u16[4] = vmlal_u8(res_u16[4], tmp[9],  BILIN_WEIGHTS_Y_4[0]);
    res_u16[4] = vmlal_u8(res_u16[4], tmp[10], BILIN_WEIGHTS_Y_4[1]);
    res_u16[4] = vmlal_u8(res_u16[4], tmp[11], BILIN_WEIGHTS_Y_4[2]);

    constexpr unsigned int SHR_QUANT = 8 + 6; // Quantize to 2b while storing
    for (unsigned int n = 0; n < h_out; ++n)
    {
        res_u16[n] = vshrq_n_u16(res_u16[n], SHR_QUANT);
        const uint8x8_t res = vmovn_u16(res_u16[n]);
        vst1_u8(&out[w_out*n], res);
    }
}

void encode_block_int(
    const uint8_t* __restrict__ inp_img,
    unsigned int block_id_x,
    unsigned int block_id_y,
    unsigned int img_w,
    uint8_t* __restrict__ out_data
){
    ZoneScopedN("enc_blk_astc");

    // Number of blocks in x-direction
    const unsigned int NB_X = img_w / BLOCK_X;

    // Input pixel block as 1D array & min/max pixel values
    uint8x16_t block_pixels_x4[BLOCK_PX_CNT / 4];
    uchar4 mincol, maxcol;

    // Read pixel block into a 1D array and select min/max at the same time
    read_block_min_max(
        inp_img,
        block_id_x,
        block_id_y,
        img_w,
        &mincol,
        &maxcol,
        block_pixels_x4
    );

    // Move the endpoints line segment such that mincol is at zero
    uchar4 ep_vec;
    _saturating_sub_4u8(maxcol, mincol, &ep_vec);

    // Inset min/max (shrink the bounding box to reduce the effect of outliers)
    const uchar4 inset = ep_vec >> 5;
    mincol = mincol + inset;
    maxcol = maxcol - inset;

    // Quantize the endpoints
    const uchar4 mincol_quant = quantize_5b(&mincol);
    const uchar4 maxcol_quant = quantize_5b(&maxcol);

    // Pixels quantized as 2b weights
    uint8_t quantized_weights[WGT_CNT] = { 0 };

    if (*as_u32(&mincol_quant) != *as_u32(&maxcol_quant))
    {
        ZoneScopedN("px_map");

        // Projection of pixels onto ep_vec to get the ideal weights
        // First, we normalize the endpoint vector and scale it.
        uint32_t ep_dot;     // Q2.16
        _saturating_dot_acc_4u8(ep_vec, ep_vec, 0, &ep_dot);

        // Q10.22, max. value 1024 for 5b quant
        const uint32_t inv_ep_dot = approx_inv_u32(ep_dot);

        // printf("ep_dot: %#010x\n", ep_dot);
        // printf("inv   : %#010x\n", inv_ep_dot);

        // The scaled ep_vec can be max. 32 due to 5b quantization
        uint32_t ep_vec32[3];
        _unpack_rgb_u32(ep_vec, ep_vec32);

        uint32_t ep_sc32[3] = {
            ep_vec32[0] * (inv_ep_dot >> 8),  // Q0.8 * Q10.14 = Q10.22
            ep_vec32[1] * (inv_ep_dot >> 8),
            ep_vec32[2] * (inv_ep_dot >> 8),
        };

        // Scaling the endpoint vector to fit 8 bits. The Q representation
        // depends on the magnitude of the division result above: between Q5.3
        // and Q0.8.
        const uint32_t log2ep[3] = {
            log2(ep_sc32[0]),
            log2(ep_sc32[1]),
            log2(ep_sc32[2]),
        };

        uint32_t sc;
        _max_u32(log2ep[0], log2ep[1], &sc);
        _max_u32(log2ep[2], sc, &sc);
        const uint32_t q_int = sc >= 21 ? sc - 21 : 0; // no. integer bits
        const uint32_t shr = 14 + q_int;               // how much to rshift

        ep_sc32[0] >>= shr;
        ep_sc32[1] >>= shr;
        ep_sc32[2] >>= shr;

        uchar4 ep_sc8;
        pack_rgb_u8(ep_sc32, &ep_sc8);

        // One (almost), represented in the variable Q format.
        const uint32_t shr_res = 8 - q_int;
        const uint32_t one = (1 << shr_res) - 1;

        // vector types
        uint8x16_t mincol_x4 = {
            mincol.x, mincol.y, mincol.z, mincol.w,
            mincol.x, mincol.y, mincol.z, mincol.w,
            mincol.x, mincol.y, mincol.z, mincol.w,
            mincol.x, mincol.y, mincol.z, mincol.w,
        };

        uint32x4_t one_x4 = { one, one, one, one };

        uint8x16_t ep_sc8_x4 = {
            ep_sc8.x, ep_sc8.y, ep_sc8.z, ep_sc8.w,
            ep_sc8.x, ep_sc8.y, ep_sc8.z, ep_sc8.w,
            ep_sc8.x, ep_sc8.y, ep_sc8.z, ep_sc8.w,
            ep_sc8.x, ep_sc8.y, ep_sc8.z, ep_sc8.w,
        };

        int32x4_t shl_res_x4 {
            -(int32_t)(shr_res),
            -(int32_t)(shr_res),
            -(int32_t)(shr_res),
            -(int32_t)(shr_res)
        };

        uint8_t ideal_weights[BLOCK_PX_CNT + (16 - BLOCK_X)] = { 0 };

        for (unsigned int i = 0; i < BLOCK_PX_CNT / 4; ++i)
        {
            uint8x16_t diff_x4 = vqsubq_u8(block_pixels_x4[i], mincol_x4);
            // Requires UDOT instruction. Apparently, -march=armv8.4-a is needed
            // as well.
            uint32x4_t dot_x4 = vdotq_u32(one_x4, diff_x4, ep_sc8_x4);

            uint32x4_t res_x4 = vshlq_u32(dot_x4, shl_res_x4);

            ideal_weights[4*i+0] = (uint8_t)(res_x4[0]);
            ideal_weights[4*i+1] = (uint8_t)(res_x4[1]);
            ideal_weights[4*i+2] = (uint8_t)(res_x4[2]);
            ideal_weights[4*i+3] = (uint8_t)(res_x4[3]);
        }

        downsample_12x12_to_8x5_u8_quant(ideal_weights, quantized_weights);
    }

    // Output buffer for quantized weights and output data
	uint8_t wgt_buf[BLOCK_SIZE] = { 0 };

    // weights ISE encoding
    constexpr uint8_t wgt_bits = 2;
    constexpr uint8_t wgt_byte_cnt = 20 * wgt_bits / 8;
    unsigned int j = 0;
    for (unsigned int i = 0; i < WGT_CNT; i += 4)
    {
        uint32_t wgt;
        _load_u32(quantized_weights + i , &wgt);

        uint8_t res = wgt & 0x03;
        res |= (wgt >> 6) & 0x0c;
        res |= (wgt >> 12) & 0x30;
        res |= (wgt >> 18) & 0xc0;

        wgt_buf[j++] = res;
    }

    // Calculate output data address (16 bytes per block)
    out_data = out_data + ((block_id_y*NB_X + block_id_x) << BLOCK_SIZE_EXP);

    // write out weights
    for (unsigned int i = 0; i < BLOCK_SIZE; i += 4)
    {
        uint32_t wgt4 = *(uint32_t*)(&wgt_buf[12 - i]);
        uint32_t reflected;
        _reflect_u32(wgt4, &reflected);
        _store_u32(&out_data[i], reflected);
    }

    // write out mode, partition, CEM
    constexpr uint16_t block_mode = 102;
    _write_bits(block_mode, 11, 0, out_data);
    constexpr unsigned int partition_count = 1;
    _write_bits(partition_count - 1, 2, 11, out_data);
    constexpr unsigned int color_format = 8;
    _write_bits(color_format, 4, 13, out_data);

    // quantized endpoint output data (layout is R0 R1 G0 G1 B0 B1)
    uint8_t mincol_quant_unpacked[3];
    uint8_t maxcol_quant_unpacked[3];

    _unpack_rgb_u8(mincol_quant, (uint8_t*)(mincol_quant_unpacked));
    _unpack_rgb_u8(maxcol_quant, (uint8_t*)(maxcol_quant_unpacked));

    uint8_t endpoints_q[6] = {
        (uint8_t)(mincol_quant_unpacked[0]),
        (uint8_t)(maxcol_quant_unpacked[0]),
        (uint8_t)(mincol_quant_unpacked[1]),
        (uint8_t)(maxcol_quant_unpacked[1]),
        (uint8_t)(mincol_quant_unpacked[2]),
        (uint8_t)(maxcol_quant_unpacked[2]),
    };

    // write out endpoints
    unsigned int off = 17;  // starting bit position for endpoints data
    constexpr uint8_t ep_bits = 5;
    off = 17;  // starting bit position for endpoints data
    for (unsigned int i = 0; i < 6; ++i)
    {
        _write_bits(endpoints_q[i], ep_bits, off, out_data);
        off += ep_bits;
    }
}

} // namespace simple::astc
