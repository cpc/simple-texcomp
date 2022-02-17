#ifndef SIMPLE_ASTC_INT_PLATFORM_H
#define SIMPLE_ASTC_INT_PLATFORM_H

#ifdef __TCE__

#include <tceops.h>
#include <tce_vector.h>

/* Defines compile-time constants in an external repo using TCE */
#include "constants.hpp"

#define ZoneNamedN(x,y,z)

// Size of the pixel blocks the image will be split into
// constexpr unsigned int BLOCK_X = 12;
// constexpr unsigned int BLOCK_Y = 12;
// constexpr unsigned int BLOCK_PX_CNT = BLOCK_X * BLOCK_Y;

inline void _load_4u8(volatile const uint8_t* inp, uchar4* out)
{
    _TCE_LD32(inp, *out);
}

inline void _load_u32(uint8_t* inp, uint32_t* out)
{
    _TCE_LD32(inp, *out);
}

inline void _store_u32(volatile uint8_t* out, uint32_t val)
{
    _TCE_ST32(out, val);
}

inline void _min_4u8(uchar4 px, uchar4 mincol, uchar4* out)
{
    _TCE_MINU8X4(px, mincol, *out);
}

inline void _max_4u8(uchar4 px, uchar4 maxcol, uchar4* out)
{
    _TCE_MAXU8X4(px, maxcol, *out);
}

inline void _min_u32(uint32_t lhs, uint32_t rhs, uint32_t* out)
{
    _TCE_MIN(lhs, rhs, *out);
}

inline void _max_u32(uint32_t lhs, uint32_t rhs, uint32_t* out)
{
    _TCE_MAX(lhs, rhs, *out);
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

inline void _shr_round_u32(uint32_t inp, uint8_t amt, uint32_t* out)
{
    _TCE_SHRRU(inp, amt, *out);
}

/* Unpack the first three elements of uchar4 vector into 32b array */
inline void _unpack_rgb_u32(uchar4 v, uint32_t out[3])
{
    constexpr uchar4 r_mask = { 1, 0, 0, 0 };
    constexpr uchar4 g_mask = { 0, 1, 0, 0 };
    constexpr uchar4 b_mask = { 0, 0, 1, 0 };

    _TCE_SATACCDOTU8X4(r_mask, v, 0, out[0]);
    _TCE_SATACCDOTU8X4(g_mask, v, 0, out[1]);
    _TCE_SATACCDOTU8X4(b_mask, v, 0, out[2]);
}

/* Unpack the first three elements of uchar4 vector into 8b array */
inline void _unpack_rgb_u8(uchar4 v, uint8_t out[3])
{
    constexpr uchar4 r_mask = { 1, 0, 0, 0 };
    constexpr uchar4 g_mask = { 0, 1, 0, 0 };
    constexpr uchar4 b_mask = { 0, 0, 1, 0 };

    _TCE_SATACCDOTU8X4(r_mask, v, 0, out[0]);
    _TCE_SATACCDOTU8X4(g_mask, v, 0, out[1]);
    _TCE_SATACCDOTU8X4(b_mask, v, 0, out[2]);
}

inline void _reflect_u32(uint32_t inp, uint32_t* out)
{
    _TCE_REFLECT(inp, *out);
}

// Taken from astcenc:
// routine to write up to 8 bits
inline void _write_bits(
	int value,
	int bitcount,
	int bitoffset,
	volatile uint8_t* ptr
) {
	int mask = (1 << bitcount) - 1;
	value &= mask;
	// ptr += bitoffset >> 3;
    uint32_t off = (bitoffset >> 3);
	bitoffset &= 7;
	value <<= bitoffset;
	mask <<= bitoffset;
	mask = ~mask;

    uint8_t a = (ptr[off] & mask) | value;
    uint8_t b = (ptr[off+1] & (mask >> 8)) | (value >> 8);

	// ptr[0] &= mask;
	// ptr[0] |= value;
	// ptr[1] &= mask >> 8;
	// ptr[1] |= value >> 8;

    _TCE_ST8(ptr + off, a);
    _TCE_ST8(ptr + off + 1, b);
}

#else  // __TCE__

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

#endif // __TCE__

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

#ifndef __TCE__
} // namespace simple::astc
#endif


#endif // SIMPLE_ASTC_INT_PLATFORM_H
