#ifndef SIMPLE_ASTC_H
#define SIMPLE_ASTC_H

#include <cstring>
#include <stdint.h>
#include <stddef.h>

#include <fstream>

#define ASTC_BLOCK_X  12
#define ASTC_BLOCK_Y  12
#define NCH_RGB  3

/* ============================================================================
	ASTC compressed file loading
    Adapted from https://github.com/ARM-software/astc-encoder
============================================================================ */
struct astc_header
{
	uint8_t magic[4];
	uint8_t block_x;
	uint8_t block_y;
	uint8_t block_z;
	uint8_t dim_x[3];			// dims = dim[0] + (dim[1] << 8) + (dim[2] << 16)
	uint8_t dim_y[3];			// Sizes are given in texels;
	uint8_t dim_z[3];			// block count is inferred
};

static const uint32_t ASTC_MAGIC_ID = 0x5CA1AB13;

int store_astc_image(
	const uint8_t* data,
	const size_t data_len,
    const int dim_x,
    const int dim_y,
	const char* filename
) {
	astc_header header;
	header.magic[0] =  ASTC_MAGIC_ID        & 0xFF;
	header.magic[1] = (ASTC_MAGIC_ID >>  8) & 0xFF;
	header.magic[2] = (ASTC_MAGIC_ID >> 16) & 0xFF;
	header.magic[3] = (ASTC_MAGIC_ID >> 24) & 0xFF;

	header.block_x = ASTC_BLOCK_X;
	header.block_y = ASTC_BLOCK_Y;
	header.block_z = 1;

	header.dim_x[0] =  dim_x        & 0xFF;
	header.dim_x[1] = (dim_x >>  8) & 0xFF;
	header.dim_x[2] = (dim_x >> 16) & 0xFF;

	header.dim_y[0] =  dim_y        & 0xFF;
	header.dim_y[1] = (dim_y >>  8) & 0xFF;
	header.dim_y[2] = (dim_y >> 16) & 0xFF;

	header.dim_z[0] =  1        & 0xFF;
	header.dim_z[1] = (1 >>  8) & 0xFF;
	header.dim_z[2] = (1 >> 16) & 0xFF;

 	std::ofstream file(filename, std::ios::out | std::ios::binary);
	if (!file)
	{
		printf("ERROR: File open failed '%s'\n", filename);
		return 1;
	}

    // printf("Saving to %s, %zu bytes\n", filename, data_len);
    // char* tmp = (char*)data;
    // for (size_t i = 0; i < data_len; ++i)
    // {
    //     printf("%3d\n", tmp[i]);
    // }

	file.write((char*)&header, sizeof(astc_header));
	file.write((char*)data, data_len);
	return 0;
}

/* ========================================================================= */

void encode_block_astc(
    const uint8_t block_pixels[NCH_RGB*ASTC_BLOCK_X*ASTC_BLOCK_Y],
    uint32_t out[4]
){
    uint8_t encoded_block_bytes[16] = {
        209, 80, 60, 7, 82, 165, 230, 115, 164, 113, 51, 128, 8, 204, 179, 179
    };

    memcpy(out, encoded_block_bytes, 16);
}

// void decode_block_astc(
//     const uint32_t enc_block[4],
//     uint8_t out_pixels[NCH_RGB*16]
// ){
// }

#endif // SIMPLE_ASTC_H
