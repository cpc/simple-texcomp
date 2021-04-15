// BC1 compression, based on:
// https://www.nvidia.com/object/real-time-ycocg-dxt-compression.html

#include "simple_texcomp.hpp"
#include "simple_mathlib.hpp"

/* */
#define INSET_MARGIN  (F(8.0) / F(255.0)) / F(16.0)

/* Convert a floating-point color into the RGB565 format */
uint32_t f32_to_rgb565(Vec3f *color)
{
    uint32_t r = (uint32_t)std::round(color->x * F(31.0));
    uint32_t g = (uint32_t)std::round(color->y * F(63.0));
    uint32_t b = (uint32_t)std::round(color->z * F(31.0));

    uint32_t out = (r << 11) | (g << 5) | b;

    r = (r << 3) | (r >> 2);
    g = (g << 2) | (g >> 4);
    b = (b << 3) | (b >> 2);

    color->x = (decimal)(r) * (F(1.0) / F(255.0));
    color->y = (decimal)(g) * (F(1.0) / F(255.0));
    color->z = (decimal)(b) * (F(1.0) / F(255.0));

    return out;
}

#if BC1_SELECT_DIAG == 1
/* Optional selection of either current or oposite diagonal - small potential
 * quality improvement at a small runtime cost
 */
void select_diagonal(
    const Vec3f block[16],
    Vec3f *mincol,
    Vec3f *maxcol
){
    Vec3f center = (*mincol + *maxcol) * F(0.5);

    Vec2f cov = { F(0.0), F(0.0) };
    for (int i = 0; i < 16; ++i) {
        Vec3f t = block[i] - center;
        cov.x += t.x * t.z;
        cov.y += t.y * t.z;
    }

    if (cov.x < F(0.0)) {
        decimal tmp = maxcol->x;
        maxcol->x = mincol->x;
        mincol->x = tmp;
    }

    if (cov.y < F(0.0)) {
        decimal tmp = maxcol->y;
        maxcol->y = mincol->y;
        mincol->y = tmp;
    }
}
#endif  // BC1_SELECT_DIAG

/* Shrink the bounding box */
void inset_bbox(Vec3f *mincol, Vec3f *maxcol)
{
    Vec3f inset = (*maxcol - *mincol) * (F(1.0) / F(16.0)) - INSET_MARGIN;
    *mincol = clamp3f(*mincol + inset, F(0.0), F(1.0));
    *maxcol = clamp3f(*maxcol - inset, F(0.0), F(1.0));
}

/* Write two 16b endpoints to a 32b integer (MSB - mincol, LSB - maxcol) */
uint32_t emit_endpoints(Vec3f *mincol, Vec3f *maxcol)
{
    uint32_t mincol_565 = f32_to_rgb565(mincol);
    uint32_t maxcol_565 = f32_to_rgb565(maxcol);

    // Swap diagonals if the other one was used
    if (maxcol_565 < mincol_565)
    {
        Vec3f tmp = *mincol;
        *mincol = *maxcol;
        *maxcol = tmp;
        return mincol_565 | (maxcol_565 << 16);
    }

    return maxcol_565 | (mincol_565 << 16);
}

/* Write 2-bit indices to a 32b integer */
uint32_t emit_indices(
    const Vec3f block[16],
    const Vec3f &mincol,
    const Vec3f &maxcol
){
    // Compute color palette
    Vec3f palette[4] = {
        maxcol,
        mincol,
        maxcol*EP_LERP2 + mincol*EP_LERP1,
        maxcol*EP_LERP1 + mincol*EP_LERP2,
    };

    // Compute indices
    uint32_t indices = 0;
    for (int i = 0; i < 16; ++i) {
        decimal dist[4];
        dist[0] = distsq3f(block[i], palette[0]);
        dist[1] = distsq3f(block[i], palette[1]);
        dist[2] = distsq3f(block[i], palette[2]);
        dist[3] = distsq3f(block[i], palette[3]);

        uint32_t b[4];
        b[0] = dist[0] > dist[3];
        b[1] = dist[1] > dist[2];
        b[2] = dist[0] > dist[2];
        b[3] = dist[1] > dist[3];

        uint32_t b4 = dist[2] > dist[3];
        uint32_t index = (b[0] & b4) | (( (b[1] & b[2]) | (b[0] & b[3]) ) << 1);

        indices |= index << (2 * i);
    }

    return indices;
}

/* Encode a block of 4x4 pixels into the BC1 format */
void encode_block_bc1(
    const uint8_t block_pixels[NCH_RGB*16],
    uint32_t out[2]
){
    // Convert the block into floating point
    Vec3f block32f[16];
    for (int i = 0; i < 16; ++i)
    {
        block32f[i].x = (decimal)block_pixels[NCH_RGB*i] / F(255.0);
        block32f[i].y = (decimal)block_pixels[NCH_RGB*i+1] / F(255.0);
        block32f[i].z = (decimal)block_pixels[NCH_RGB*i+2] / F(255.0);
    }

    // Determine line through color space
    Vec3f mincol, maxcol;
    find_minmaxcolor_bbox(block32f, &mincol, &maxcol);
#if BC1_SELECT_DIAG == 1
    select_diagonal(block32f, &mincol, &maxcol);
#endif  // BC1_SELECT_DIAG
    inset_bbox(&mincol, &maxcol);

    // Write endpoints
    out[0] = emit_endpoints(&mincol, &maxcol);

    // Write indices
    out[1] = emit_indices(block32f, mincol, maxcol);
}

/* Decode an encoded block into an array of 16 pixels */
void decode_block_bc1(
    const uint32_t enc_block[2],
    uint8_t out_pixels[NCH_RGB*16]
){
    Vec3f palette[4];
    palette[0] = rgb565_to_f32(enc_block[0] & 0xffff);  // maxcol
    palette[1] = rgb565_to_f32((enc_block[0] >> 16) & 0xffff);  // mincol
    palette[2] = palette[0]*EP_LERP2 + palette[1]*EP_LERP1;
    palette[3] = palette[0]*EP_LERP1 + palette[1]*EP_LERP2;

    for (int i = 0; i < 16; ++i)
    {
        int id = ( enc_block[1] >> (2*i) ) & 0x3;
        Vec3f res = palette[id];
        out_pixels[NCH_RGB*i] = (uint8_t)(res.x * F(255.0));
        out_pixels[NCH_RGB*i+1] = (uint8_t)(res.y * F(255.0));
        out_pixels[NCH_RGB*i+2] = (uint8_t)(res.z * F(255.0));
    }
}
