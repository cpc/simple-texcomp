// BC1 compression, based on:
// https://www.nvidia.com/object/real-time-ycocg-dxt-compression.html

#ifndef SIMPLE_BC1_H
#define SIMPLE_BC1_H

#include "simple_bcn_common.h"

/* Optional refinement (can improve quality at small runtime cost)
 * 1 - enable, 0 - disable
 */
#define SELECT_DIAG  1

/* Rounding the bounding box inset outwards (= (8.0/255.0)/16.0) */
#define INSET_MARGIN  (8.0f / 255.0f) / 16.0f

/* Convert a floating-point color into the RGB565 format */
static inline uint32_t f32_to_rgb565(Vec3f *color)
{
    uint32_t r = (uint32_t)round(color->x * 31.0f);
    uint32_t g = (uint32_t)round(color->y * 63.0f);
    uint32_t b = (uint32_t)round(color->z * 31.0f);

    uint32_t out = (r << 11) | (g << 5) | b;

    r = (r << 3) | (r >> 2);
    g = (g << 2) | (g >> 4);
    b = (b << 3) | (b >> 2);

    color->x = (double)(r) * (1.0f / 255.0f);
    color->y = (double)(g) * (1.0f / 255.0f);
    color->z = (double)(b) * (1.0f / 255.0f);

    return out;
}

/* Convert a color from RGB565 format into floating point */
static inline Vec3f rgb565_to_f32(uint16_t color)
{
    return Vec3f {
        (double)( (color >> 11) & 0x1f ) / 0x1f,
        (double)( (color >> 5) & 0x3f ) / 0x3f,
        (double)( color & 0x1f ) / 0x1f,
    };
}

#if SELECT_DIAG == 1
/* Optional selection of either current or oposite diagonal - small performance
 * cost and minimal quality improvement. */
static inline void select_diagonal(
    const Vec3f block[16],
    Vec3f *mincol,
    Vec3f *maxcol
){
    Vec3f center = (*mincol + *maxcol) * 0.5f;

    Vec2f cov = { 0.0f, 0.0f };
    for (int i = 0; i < 16; ++i) {
        Vec3f t = block[i] - center;
        cov.x += t.x * t.z;
        cov.y += t.y * t.z;
    }

#ifndef NDEBUG
    printf("## cov: %16.14f %16.14f\n", cov.x, cov.y);
#endif

    if (cov.x < 0.0f) {
        double tmp = maxcol->x;
        maxcol->x = mincol->x;
        mincol->x = tmp;
    }

    if (cov.y < 0.0f) {
        double tmp = maxcol->y;
        maxcol->y = mincol->y;
        mincol->y = tmp;
    }
}
#endif  // SELECT_DIAG

/* Shrink the bounding box (by half the distance of equidistant points) to
 * eliminate influence of outliers */
static inline void inset_bbox(Vec3f *mincol, Vec3f *maxcol)
{
    Vec3f inset = (*maxcol - *mincol) * (1.0f / 16.0f) - INSET_MARGIN;
    *mincol = clamp3f(*mincol + inset, 0.0f, 1.0f);
    *maxcol = clamp3f(*maxcol - inset, 0.0f, 1.0f);
}

/* Write two 16b endpoints to a 32b integer (MSB - mincol, LSB - maxcol) */
static inline uint32_t emit_endpoints(Vec3f *mincol, Vec3f *maxcol)
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
static inline uint32_t emit_indices(
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
        double dist[4];
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
        block32f[i].x = (double)block_pixels[NCH_RGB*i] / 255.0f;
        block32f[i].y = (double)block_pixels[NCH_RGB*i+1] / 255.0f;
        block32f[i].z = (double)block_pixels[NCH_RGB*i+2] / 255.0f;
    }

    // Determine line through color space
    Vec3f mincol, maxcol;
    find_minmaxcolor_bbox(block32f, &mincol, &maxcol);
#if SELECT_DIAG == 1
    select_diagonal(block32f, &mincol, &maxcol);
#endif  // SELECT_DIAG
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
        out_pixels[NCH_RGB*i] = (uint8_t)(res.x * 255.0f);
        out_pixels[NCH_RGB*i+1] = (uint8_t)(res.y * 255.0f);
        out_pixels[NCH_RGB*i+2] = (uint8_t)(res.z * 255.0f);
    }
}

#endif /* SIMPLE_BC1_H */
