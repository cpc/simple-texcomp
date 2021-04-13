// YCoCg-BC3 compression, based on:
// https://www.nvidia.com/object/real-time-ycocg-dxt-compression.html

#include "simple_texcomp.hpp"
#include "simple_mathlib.hpp"

/*
 * = (8.0/255.0)/16.0 for CoCg and (16.0/255.0)/32.0 for Y */
#define INSET_MARGIN_Y  (F(16.0) / F(255.0)) / F(32.0)
#define INSET_MARGIN_COCG  (F(8.0) / F(255.0)) / F(16.0)

/* Offset to center color value at grey level (= 128.0/255.0) */
#define OFFSET  F(128.0) / F(255.0)

/* Helper type for easier access of 16b and 32b words inside a 64b block
 *
 * This is very useful for the Y block, where you need to write 48 bits of
 * indices which does not fit uint. */
union bc_block_t
{
    uint8_t b8[8];
    uint16_t b16[4];
    uint32_t b32[2];
    uint64_t b64;
};

/* Convert RGB pixel to YCoCg */
static inline Vec3f rgb_to_ycocg(const Vec3f &rgb)
{
    return Vec3f {
        rgb.dot(Vec3f {  F(0.25), F(0.50),  F(0.25) }),
        rgb.dot(Vec3f {  F(0.50), F(0.00), -F(0.50) }) + OFFSET,
        rgb.dot(Vec3f { -F(0.25), F(0.50), -F(0.25) }) + OFFSET,
    };
}

/* Convert a f32 2-channel color into RG channels and scale into B channel */
static inline uint32_t f32scale_to_rgb565(Vec2f *color, uint32_t scale)
{
    uint32_t r = (uint32_t)std::round(color->x * F(31.0));
    uint32_t g = (uint32_t)std::round(color->y * F(63.0));

    uint32_t out = (r << 11) | (g << 5) | (scale - 1);

    r = (r << 3) | (r >> 2);
    g = (g << 2) | (g >> 4);

    color->x = (decimal)(r) * (F(1.0) / F(255.0));
    color->y = (decimal)(g) * (F(1.0) / F(255.0));

    return out;
}

/* Select either current or oposite diagonal */
static inline void select_cocg_diagonal(
    const Vec3f block[16],
    Vec2f *min_cocg,
    Vec2f *max_cocg
){
    Vec2f center = (*min_cocg + *max_cocg) * F(0.5);

    decimal cov = F(0.0);
    for (int i = 0; i < 16; ++i) {
        Vec2f t = {
            block[i].y - center.x,
            block[i].z - center.y,
        };
        cov += t.x * t.y;
    }

    if (cov < 0) {
        decimal tmp = max_cocg->y;
        max_cocg->y = min_cocg->y;
        min_cocg->y = tmp;
    }
}

/* Scale values up in case of low dynamic range */
static inline uint32_t get_cocg_scale(
    const Vec2f &min_cocg,
    const Vec2f &max_cocg
){
    Vec2f m0 = abs2f(min_cocg - OFFSET);
    Vec2f m1 = abs2f(max_cocg - OFFSET);

    decimal m = (decimal)std::fmax(std::fmax(m0.x, m0.y), std::fmax(m1.x, m1.y));

    const decimal s0 = F(64.0) / F(255.0);
    const decimal s1 = F(32.0) / F(255.0);

    uint32_t scale = 1;
    if (m < s0) {
        scale = 2;
    }
    if (m < s1) {
        scale = 4;
    }

    return scale;
}

/* Shrink the bounding box around the Co and Cg channels */
void inset_bbox_cocg(
    Vec2f *min_cocg,
    Vec2f *max_cocg
){
    Vec2f inset = (*max_cocg - *min_cocg) * (F(1.0) / F(16.0)) - INSET_MARGIN_COCG;
    *min_cocg = clamp2f(*min_cocg + inset, F(0.0), F(1.0));
    *max_cocg = clamp2f(*max_cocg - inset, F(0.0), F(1.0));
}

/* Write endpoints to the BC1 block */
void emit_endpoints_cocg(
    Vec2f *min_cocg,
    Vec2f *max_cocg,
    uint32_t scale,
    bc_block_t *out_block
){
    // Scale
    *min_cocg = (*min_cocg - OFFSET) * scale + OFFSET;
    *max_cocg = (*max_cocg - OFFSET) * scale + OFFSET;

    // Inset
    inset_bbox_cocg(min_cocg, max_cocg);

    uint32_t min_cocg_565 = f32scale_to_rgb565(min_cocg, scale);
    uint32_t max_cocg_565 = f32scale_to_rgb565(max_cocg, scale);

    // Swap diagonals (shouldn't be necessary but some decoders treat BC3 CoCg
    // block as BC1 also considering punch-through alpha)
    if (max_cocg_565 < min_cocg_565) {
        uint32_t tmp_565 = min_cocg_565;
        min_cocg_565 = max_cocg_565;
        max_cocg_565 = tmp_565;
        Vec2f tmp = *min_cocg;
        *min_cocg = *max_cocg;
        *max_cocg = tmp;
    }

    out_block->b16[0] = max_cocg_565;
    out_block->b16[1] = min_cocg_565;

    // Rescale
    *min_cocg = (*min_cocg - OFFSET) / (decimal)scale + OFFSET;
    *max_cocg = (*max_cocg - OFFSET) / (decimal)scale + OFFSET;
}

/* Write 2-bit indices to the BC1 block */
void emit_indices_cocg(
    const Vec3f block[16],
    const Vec2f &min_cocg,
    const Vec2f &max_cocg,
    bc_block_t *out_block
){
    // Compute color palette ( mix(x,y,a) = x + (y-x)*a )
    Vec2f lerp1 = max_cocg*EP_LERP2 + min_cocg*EP_LERP1;
    Vec2f lerp2 = max_cocg*EP_LERP1 + min_cocg*EP_LERP2;

    // Compute indices
    uint32_t indices = 0;
    for (int i = 0; i < 16; ++i) {
        Vec2f block_cocg = Vec2f { block[i].y, block[i].z };
        decimal dist[4];
        dist[0] = distsq2f(block_cocg, max_cocg);
        dist[1] = distsq2f(block_cocg, min_cocg);
        dist[2] = distsq2f(block_cocg, lerp1);
        dist[3] = distsq2f(block_cocg, lerp2);

        uint32_t b[4];
        b[0] = dist[0] > dist[3];
        b[1] = dist[1] > dist[2];
        b[2] = dist[0] > dist[2];
        b[3] = dist[1] > dist[3];
        uint32_t b4 = dist[2] > dist[3];
        uint32_t index = (b[0] & b4) | (( (b[1] & b[2]) | (b[0] & b[3]) ) << 1);

        indices |= index << (2 * i);
    }

    out_block->b32[1] = indices;
}

/* Shrink the bounding box around the Y channel */
void inset_bbox_y(
    decimal *min_y,
    decimal *max_y
){
    decimal inset = (*max_y - *min_y) / F(32.0) - INSET_MARGIN_Y;
    *min_y = fclamp(*min_y + inset, F(0.0), F(1.0));
    *max_y = fclamp(*max_y - inset, F(0.0), F(1.0));
}

/* Write Y endpoints into the BC4 block, each 8 bits */
void emit_endpoints_y(
    decimal *min_y,
    decimal *max_y,
    bc_block_t *out_block
){
    inset_bbox_y(min_y, max_y);

    out_block->b8[0] = (uint8_t)(std::round(*max_y * F(255.0)));
    out_block->b8[1] = (uint8_t)(std::round(*min_y * F(255.0)));
}

/* Write 3-bit Y indices into the rest of the BC4 block */
void emit_indices_y(
    const Vec3f block[16],
    const decimal &min_y,
    const decimal &max_y,
    bc_block_t *out_block
){
    decimal mid = (max_y - min_y) / (F(2.0) * F(7.0));

    decimal ab1 = min_y + mid;
    decimal ab2 = (F(6.0) * max_y + F(1.0) * min_y) * (F(1.0) / F(7.0)) + mid;
    decimal ab3 = (F(5.0) * max_y + F(2.0) * min_y) * (F(1.0) / F(7.0)) + mid;
    decimal ab4 = (F(4.0) * max_y + F(3.0) * min_y) * (F(1.0) / F(7.0)) + mid;
    decimal ab5 = (F(3.0) * max_y + F(4.0) * min_y) * (F(1.0) / F(7.0)) + mid;
    decimal ab6 = (F(2.0) * max_y + F(5.0) * min_y) * (F(1.0) / F(7.0)) + mid;
    decimal ab7 = (F(1.0) * max_y + F(6.0) * min_y) * (F(1.0) / F(7.0)) + mid;

    bc_block_t indices;
    indices.b64 = 0;
    uint64_t index;
    for (int i = 0; i < 16; ++i)
    {
        decimal a = block[i].x;
        index = 1;
        index += (a <= ab1);
        index += (a <= ab2);
        index += (a <= ab3);
        index += (a <= ab4);
        index += (a <= ab5);
        index += (a <= ab6);
        index += (a <= ab7);
        index &= 7;
        index ^= (2 > index);
        indices.b64 |= index << (3 * i + 16);
    }

    out_block->b16[1] = indices.b16[1];
    out_block->b32[1] = indices.b32[1];
}

/* Encode a block of 4x4 pixels into the YCoCg-BC3 format */
void encode_block_ycocg_bc3(
    const uint8_t block_pixels[NCH_RGB*16],
    uint32_t out[4]
){
    // Convert the block into floating point and YCoCg color space
    Vec3f block32f_ycocg[16];
    for (int i = 0; i < 16; ++i)
    {
        block32f_ycocg[i] = rgb_to_ycocg(
            Vec3f {
                (decimal)block_pixels[NCH_RGB*i] / F(255.0),
                (decimal)block_pixels[NCH_RGB*i+1] / F(255.0),
                (decimal)block_pixels[NCH_RGB*i+2] / F(255.0),
            }
        );
    }

    // Determine line through color space
    Vec3f mincol, maxcol;
    find_minmaxcolor_bbox(block32f_ycocg, &mincol, &maxcol);

    decimal min_y = mincol.x;
    Vec2f min_cocg = { mincol.y, mincol.z };
    decimal max_y = maxcol.x;
    Vec2f max_cocg = { maxcol.y, maxcol.z };

    select_cocg_diagonal(block32f_ycocg, &min_cocg, &max_cocg);

    // Write CoCg into BC1 block
    bc_block_t out_block_cocg;
    uint32_t scale = get_cocg_scale(min_cocg, max_cocg);
    emit_endpoints_cocg(&min_cocg, &max_cocg, scale, &out_block_cocg);
    emit_indices_cocg(block32f_ycocg, min_cocg, max_cocg, &out_block_cocg);

    // Write Y into BC4 block
    bc_block_t out_block_y;
    emit_endpoints_y(&min_y, &max_y, &out_block_y);
    emit_indices_y(block32f_ycocg, min_y, max_y, &out_block_y);

    // Write local blocks into the output block
    out[0] = out_block_y.b32[0];
    out[1] = out_block_y.b32[1];
    out[2] = out_block_cocg.b32[0];
    out[3] = out_block_cocg.b32[1];
}

/* Decode an encoded block into an array of 16 pixels */
void decode_block_ycocg_bc3(
    const uint32_t enc_block[4],
    uint8_t out_pixels[NCH_RGB*16]
){
    // Read Y palette
    bc_block_t block_y = {
        .b32 = { enc_block[0], enc_block[1] }
    };

    decimal max_y = (decimal)block_y.b8[0] / F(255.0);
    decimal min_y = (decimal)block_y.b8[1] / F(255.0);
    decimal palette_y[8];
    palette_y[0] = max_y;
    palette_y[1] = min_y;
    palette_y[2] = (F(6.0) * max_y + F(1.0) * min_y) * (F(1.0) / F(7.0));
    palette_y[3] = (F(5.0) * max_y + F(2.0) * min_y) * (F(1.0) / F(7.0));
    palette_y[4] = (F(4.0) * max_y + F(3.0) * min_y) * (F(1.0) / F(7.0));
    palette_y[5] = (F(3.0) * max_y + F(4.0) * min_y) * (F(1.0) / F(7.0));
    palette_y[6] = (F(2.0) * max_y + F(5.0) * min_y) * (F(1.0) / F(7.0));
    palette_y[7] = (F(1.0) * max_y + F(6.0) * min_y) * (F(1.0) / F(7.0));

    // Read CoCg palette
    bc_block_t block_cocg = {
        .b32 = { enc_block[2], enc_block[3] }
    };

    Vec3f palette_cocgscale[4];
    palette_cocgscale[0] = rgb565_to_f32(block_cocg.b16[0]);  // maxcol
    palette_cocgscale[1] = rgb565_to_f32(block_cocg.b16[1]);  // mincol
    palette_cocgscale[2] =
        palette_cocgscale[0]*EP_LERP2 + palette_cocgscale[1]*EP_LERP1;
    palette_cocgscale[3] =
        palette_cocgscale[0]*EP_LERP1 + palette_cocgscale[1]*EP_LERP2;

    // Interpolate pixels
    for (int i = 0; i < 16; ++i)
    {
        int id_cocg = ( block_cocg.b32[1] >> (2*i) ) & 0x3;
        Vec3f cocgscale = palette_cocgscale[id_cocg];

        decimal inv_scale = F(1.0) / ((F(255.0) / F(8.0)) * cocgscale.z + F(1.0));
        decimal co = (cocgscale.x - OFFSET) * inv_scale;
        decimal cg = (cocgscale.y - OFFSET) * inv_scale;

        int id_y = ( block_y.b64 >> (16 + 3*i) ) & 0x7;
        decimal y = palette_y[id_y];

		// decimal r = y + co - cg;
		// decimal g = y + cg;
		// decimal b = y - co - cg;

		decimal r = fclamp(y + co - cg, F(0.0), F(1.0));
		decimal g = fclamp(y + cg, F(0.0), F(1.0));
		decimal b = fclamp(y - co - cg, F(0.0), F(1.0));

        out_pixels[NCH_RGB*i] = (uint8_t)(r * F(255.0));
        out_pixels[NCH_RGB*i+1] = (uint8_t)(g * F(255.0));
        out_pixels[NCH_RGB*i+2] = (uint8_t)(b * F(255.0));
    }
}
