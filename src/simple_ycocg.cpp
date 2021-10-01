// RGB -> YCoCg conversion without compression

#include "simple_texcomp.hpp"
#include "simple_mathlib.hpp"

namespace simple::ycocg {

/* Offset to center color value at grey level (= 128.0/255.0) */
constexpr decimal OFFSET = F(128.0) / F(255.0);

/* Convert RGB pixel to YCoCg */
inline Vec3f rgb_to_ycocg(const Vec3f &rgb)
{
    return Vec3f {
        rgb.dot(Vec3f {  F(0.25), F(0.50),  F(0.25) }),
        rgb.dot(Vec3f {  F(0.50), F(0.00), -F(0.50) }) + OFFSET,
        rgb.dot(Vec3f { -F(0.25), F(0.50), -F(0.25) }) + OFFSET,
    };
}

/* Encode a block of 4x4 pixels into the YCoCg format (R-Y, G-Co, B-Cg) */
void encode_block(
    const uint8_t block_pixels[NCH_RGB * 16],
    uint32_t out[16]  // 1 byte per channel
){
    for (int i = 0; i < 16; ++i)
    {
        // YCoCg-R transform

        int r = block_pixels[NCH_RGB*i+0];
        int g = block_pixels[NCH_RGB*i+1];
        int b = block_pixels[NCH_RGB*i+2];

		int co  = r - b;
        int tmp = b + (co >> 1);
        int cg  = g - tmp;
        int y   = tmp + (cg >> 1);

        uint8_t y_u8 = (uint8_t)y;  // range 0..255

        // lossy operation, would need one more bit for storing these...
        uint8_t co_u8 = co >> 1;  // range -255..255 -> 0..255
        uint8_t cg_u8 = cg >> 1;  // range -255..255 -> 0..255

        uint32_t out_px = 0;
        out_px |= (uint32_t)(y_u8);
        out_px |= (uint32_t)(co_u8) << 8;
        out_px |= (uint32_t)(cg_u8) << 16;
        out_px |= 0xff << 24;

        out[i] = out_px;
    }
}

/* Decode an encoded block into an array of 16 pixels */
void decode_block(
    const uint32_t enc_block[16],
    uint8_t out_pixels[NCH_RGB*16]
){
    for (int i = 0; i < 16; ++i)
    {
        int y  = (int)(enc_block[i] & 0xff);
        int co = (int)((int8_t)( (enc_block[i] >>  8) & 0xff )) << 1;
        int cg = (int)((int8_t)( (enc_block[i] >> 16) & 0xff )) << 1;

        int tmp = y - (cg >> 1);
        int g   = iclamp(cg + tmp, 0, 255);
        int b   = iclamp(tmp - (co >> 1), 0, 255);
        int r   = iclamp(b + co, 0, 255);

        out_pixels[NCH_RGB*i]   = (uint8_t)(r);
        out_pixels[NCH_RGB*i+1] = (uint8_t)(g);
        out_pixels[NCH_RGB*i+2] = (uint8_t)(b);
        if (NCH_RGB == 4)
        {
            out_pixels[NCH_RGB*i+3] = 255;
        }
    }
}

} // namespace simple::ycocg
