#ifndef SIMPLE_TEXCOMP_HPP
#define SIMPLE_TEXCOMP_HPP

#include <cstdint>

/* Number of channels if a RGB pixel. Used for calculating data sizes */
#define NCH_RGB  3

/******************************************************************************
 * ASTC
 *
 * Decoding function is purposefully missing. Instead, ARM decoder should be
 * used.
 *****************************************************************************/

/* Input ASTC block size */
#define ASTC_BLOCK_X  12
#define ASTC_BLOCK_Y  12

/** Store compressed file into .astc format
 *
 * Taken from https://github.com/ARM-software/astc-encoder
 */
int store_astc_image(
	const uint8_t* data,
	const uint32_t data_len,
    const int dim_x,
    const int dim_y,
	const char* filename
);

void encode_block_astc(
    const uint8_t block_pixels[NCH_RGB*ASTC_BLOCK_X*ASTC_BLOCK_Y],
    uint32_t out[4]
);

// void decode_block_astc(
//     const uint32_t enc_block[4],
//     uint8_t out_pixels[NCH_RGB*16]
// );


/******************************************************************************
 * BC1
 *****************************************************************************/

/* Optional refinement (can improve quality at small runtime cost)
 * 1 - enable, 0 - disable
 */
#define BC1_SELECT_DIAG  1

/* Encode a block of 4x4 pixels into the BC1 format */
void encode_block_bc1(
    const uint8_t block_pixels[NCH_RGB*16],
    uint32_t out[2]
);

/* Decode an encoded block into an array of 16 pixels */
void decode_block_bc1(
    const uint32_t enc_block[2],
    uint8_t out_pixels[NCH_RGB*16]
);


/******************************************************************************
 * YCoCg-BC3
 *****************************************************************************/

/* Encode a block of 4x4 pixels into the YCoCg-BC3 format */
void encode_block_ycocg_bc3(
    const uint8_t block_pixels[NCH_RGB*16],
    uint32_t out[4]
);

/* Decode an encoded block into an array of 16 pixels */
void decode_block_ycocg_bc3(
    const uint32_t enc_block[4],
    uint8_t out_pixels[NCH_RGB*16]
);


#endif // SIMPLE_TEXCOMP_HPP
