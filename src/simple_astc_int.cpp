#include <stdint.h>

#ifdef __TCE__

#include <tceops.h>
#include <tce_vector.h>

/* Defines compile-time constants in the Aamu repo */
#include "constants.hpp"

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

#else  // not __TCE__

#include "simple_mathlib.hpp"
#include "simple_texcomp.hpp"

using namespace simple;

typedef Vec4u8 uchar4;

constexpr unsigned int BLOCK_PX_CNT = astc::BLOCK_X * astc::BLOCK_Y;

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

namespace simple::astc {

#endif  // __TCE__

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

/* Pack the first three elements of input 32b array into uchar4 vector */
inline void pack_rgb_u8(uint32_t a[3], uchar4* out)
{
    out->x = (uint8_t)(a[0]);
    out->y = (uint8_t)(a[1]);
    out->z = (uint8_t)(a[2]);
}

/** Downsample 12x12 block into 8x5 using bilinear filtering and quantize it
 *
 * The inp/out grid size is fixed and one channel per sample is assumed.
 * The output is also quantized.
 */
inline void downsample_12x12_to_8x5_u8_quant(
    const uint8_t inp[BLOCK_PX_CNT],
    uint8_t out[WGT_CNT]
){
    constexpr unsigned int w_inp = 12;
    constexpr unsigned int h_inp = 12;
    constexpr unsigned int w_out = 8;
    constexpr unsigned int h_out = 5;

    // size of the bilin. kernel
    constexpr unsigned int pixel_count_x = 3;
    constexpr unsigned int pixel_count_y = 6;

    // Buffer for holding intermediate results (columns stored as rows for easy
    // access in the second loop).
    // volatile uint8_t tmp[w_out*h_inp];
    volatile uint8_t tmp[h_inp*w_out];

    // First, interpolate rows.
    for (unsigned int y = 0; y < h_inp; ++y)
    {
        const unsigned int base_addr_i = y*w_inp;

        uchar4 inp0, inp4, inp8;
        _load_4u8(inp + base_addr_i + 0, &inp0);
        _load_4u8(inp + base_addr_i + 4, &inp4);
        _load_4u8(inp + base_addr_i + 8, &inp8);

        // m = 0
        const uchar4 wgt0 = BILIN_WEIGHTS_X_8_U8[0];

        uint32_t out0;
        _saturating_dot_acc_4u8(inp0, wgt0, 0, &out0);

        unsigned int addr_o = y;
        tmp[addr_o] = (uint8_t)(out0 >> 8);

        // m = 1
        const uchar4 wgt1 = BILIN_WEIGHTS_X_8_U8[1];

        uint32_t out1;
        _saturating_dot_acc_4u8(inp0, wgt1, 0, &out1);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out1 >> 8);

        // m = 2
        const uchar4 wgt2a = BILIN_WEIGHTS_X_8_U8[2];
        const uchar4 wgt2b = BILIN_WEIGHTS_X_8_U8[3];

        uint32_t out2;
        _saturating_dot_acc_4u8(inp0, wgt2a, 0, &out2);
        _saturating_dot_acc_4u8(inp4, wgt2b, out2, &out2);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out2 >> 8);

        // m = 3
        const uchar4 wgt3 = BILIN_WEIGHTS_X_8_U8[4];

        uint32_t out3;
        _saturating_dot_acc_4u8(inp4, wgt3, 0, &out3);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out3 >> 8);

        // m = 4
        const uchar4 wgt4 = BILIN_WEIGHTS_X_8_U8[5];

        uint32_t out4;
        _saturating_dot_acc_4u8(inp4, wgt4, 0, &out4);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out4 >> 8);

        // m = 5
        const uchar4 wgt5a = BILIN_WEIGHTS_X_8_U8[6];
        const uchar4 wgt5b = BILIN_WEIGHTS_X_8_U8[7];

        uint32_t out5;
        _saturating_dot_acc_4u8(inp4, wgt5a, 0, &out5);
        _saturating_dot_acc_4u8(inp8, wgt5b, out5, &out5);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out5 >> 8);

        // m = 6
        const uchar4 wgt6 = BILIN_WEIGHTS_X_8_U8[8];

        uint32_t out6;
        _saturating_dot_acc_4u8(inp8, wgt6, 0, &out6);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out6 >> 8);

        // m = 7
        const uchar4 wgt7 = BILIN_WEIGHTS_X_8_U8[9];

        uint32_t out7;
        _saturating_dot_acc_4u8(inp8, wgt7, 0, &out7);

        addr_o += h_inp;
        tmp[addr_o] = (uint8_t)(out7 >> 8);
    }

    // Next, columns
    constexpr unsigned int SHR_QUANT = 8 + 6; // Quantize to 2b while storing
    for (unsigned int m = 0; m < w_out; ++m)
    {
        const unsigned int base_addr_i = m*h_inp;

        uchar4 tmp0, tmp4, tmp8;
        _load_4u8(tmp + base_addr_i + 0, &tmp0);
        _load_4u8(tmp + base_addr_i + 4, &tmp4);
        _load_4u8(tmp + base_addr_i + 8, &tmp8);

        // n = 0
        const uchar4 wgt0 = BILIN_WEIGHTS_Y_5_U8[0];

        uint32_t out0;
        _saturating_dot_acc_4u8(tmp0, wgt0, 0, &out0);

        unsigned int addr_o = m;
        out[addr_o] = (uint8_t)(out0 >> SHR_QUANT);

        // n = 1
        const uchar4 wgt1a = BILIN_WEIGHTS_Y_5_U8[1];
        const uchar4 wgt1b = BILIN_WEIGHTS_Y_5_U8[2];

        uint32_t out1;
        _saturating_dot_acc_4u8(tmp0, wgt1a, 0, &out1);
        _saturating_dot_acc_4u8(tmp4, wgt1b, out1, &out1);

        addr_o += w_out;
        out[addr_o] = (uint8_t)(out1 >> SHR_QUANT);

        // n = 2
        const uchar4 wgt2a = BILIN_WEIGHTS_Y_5_U8[3];
        const uchar4 wgt2b = BILIN_WEIGHTS_Y_5_U8[4];
        const uchar4 wgt2c = BILIN_WEIGHTS_Y_5_U8[5];

        uint32_t out2;
        _saturating_dot_acc_4u8(tmp0, wgt2a, 0, &out2);
        _saturating_dot_acc_4u8(tmp4, wgt2b, out2, &out2);
        _saturating_dot_acc_4u8(tmp8, wgt2c, out2, &out2);

        addr_o += w_out;
        out[addr_o] = (uint8_t)(out2 >> SHR_QUANT);

        // n = 3
        const uchar4 wgt3a = BILIN_WEIGHTS_Y_5_U8[6];
        const uchar4 wgt3b = BILIN_WEIGHTS_Y_5_U8[7];

        uint32_t out3;
        _saturating_dot_acc_4u8(tmp4, wgt3a, 0, &out3);
        _saturating_dot_acc_4u8(tmp8, wgt3b, out3, &out3);

        addr_o += w_out;
        out[addr_o] = (uint8_t)(out3 >> SHR_QUANT);

        // n = 4
        const uchar4 wgt4 = BILIN_WEIGHTS_Y_5_U8[8];

        uint32_t out4;
        _saturating_dot_acc_4u8(tmp8, wgt4, 0, &out4);

        addr_o += w_out;
        out[addr_o] = (uint8_t)(out4 >> SHR_QUANT); // TODO: rounding
    }
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

#else  // not __TCE__

void encode_block_int(
    const uint8_t* __restrict__ inp_img,
    unsigned int block_id_x,
    unsigned int block_id_y,
    unsigned int img_w,
    uint8_t* __restrict__ out_data
){

    // Number of blocks in x-direction
    const unsigned int NB_X = img_w / BLOCK_X;

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
    const uchar4 inset = ep_vec >> 5;
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

        // printf("ep_sc32[0]: %#010x\n", ep_sc32[0]);
        // printf("ep_sc32[1]: %#010x\n", ep_sc32[1]);
        // printf("ep_sc32[2]: %#010x\n", ep_sc32[2]);
        // printf("ep_sc8    : %#010x\n", *as_u32(&ep_sc8));
        // printf("one       : %#010x\n", one);
        // print_minmax_u8("          :", mincol, maxcol);

        // Pixels mapped to the endpoint line without quantization and
        // downsampling
        uint8_t ideal_weights[BLOCK_PX_CNT];

        for (unsigned int i = 0; i < BLOCK_PX_CNT; ++i)
        {
            uchar4 diff;
            _saturating_sub_4u8(block_pixels[i], mincol, &diff);

            // dot product max: Q5.3 * Q0.8 + 3xADD = Q8.11 -> sat 5.11 (0.11)
            // dot product min: Q0.8 * Q0.8 + 3xADD = Q3.16 -> sat 0.16
            // recover lost precision by adding one
            uint32_t res;
            _saturating_dot_acc_4u8(diff, ep_sc8, one, &res);
            ideal_weights[i] = (uint8_t)(res >> shr_res);  // -> Q0.8

            // printf("iwgt[%3d] : %#04x\n", i, ideal_weights[i]);
        }

        downsample_12x12_to_8x5_u8_quant(
            ideal_weights,
            quantized_weights
        );
    }

    // Output buffers for quantized weights and output data
	uint8_t wgt_buf[BLOCK_SIZE] = { 0 };

    // weights ISE encoding
    constexpr uint8_t wgt_bits = 2;
    constexpr uint8_t wgt_byte_cnt = 10; // 40 * 2 / 8
    unsigned int j = 0;
    for (unsigned int i = 0; i < WGT_CNT; i += 4)
    {
        uint32_t wgt;
        _load_u32(quantized_weights + i , &wgt);

        uint8_t res = wgt & 0x03;
        res |= (wgt >> 6) & 0x0c;
        res |= (wgt >> 12) & 0x30;
        res |= (wgt >> 18) & 0xc0;

        //debug
        //_TCE_ST8(TEST+j, res);

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

#ifndef __TCE__
} // namespace simple::astc
#endif
