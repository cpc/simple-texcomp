#include <stdint.h>

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

inline void _shr_round_4u8(uchar4 inp, uchar4 amt, uchar4* out)
{
    _TCE_SHRRU8X4(inp, amt, *out);
}

#else  // __TCE__

#include "simple_mathlib.hpp"
#include "simple_texcomp.hpp"

using namespace simple;

typedef Vec4u8 uchar4;

constexpr unsigned int BLOCK_PX_CNT = astc::BLOCK_X * astc::BLOCK_Y;

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

namespace simple::astc {

#endif  // __TCE__


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

    // Move the endpoints line segment such that mincol is at zero
    uchar4 ep_vec;
    _saturating_sub_4u8(maxcol, mincol, &ep_vec);

    // Inset min/max (shrink the bounding box to reduce the effect of outliers)
    const uchar4 inset = ep_vec >> 4;
    mincol = mincol + inset;
    maxcol = maxcol - inset;

    // Quantize the endpoints
    const uchar4 mincol_quant = quantize_5b(&mincol);
    const uchar4 maxcol_quant = quantize_5b(&maxcol);
}

#ifndef __TCE__
} // namespace simple::astc
#endif
