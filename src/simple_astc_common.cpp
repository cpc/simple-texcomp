#include <cstdint>
#include <fstream>

#include "platform.hpp"
#include "simple_texcomp.hpp"

namespace simple::astc {

/* ============================================================================
	ASTC compressed file handling
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

const uint32_t MAGIC_ID = 0x5CA1AB13;

int store_astc_image(
	const uint8_t* data,
	const uint32_t data_len,
    const int dim_x,
    const int dim_y,
	const char* filename
) {
	astc_header header;
	header.magic[0] =  MAGIC_ID        & 0xFF;
	header.magic[1] = (MAGIC_ID >>  8) & 0xFF;
	header.magic[2] = (MAGIC_ID >> 16) & 0xFF;
	header.magic[3] = (MAGIC_ID >> 24) & 0xFF;

	header.block_x = BLOCK_X;
	header.block_y = BLOCK_Y;
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
		LOGE("ERROR: File open failed '%s'\n", filename);
		return 1;
	}

    // LOGE("Saving to %s, %zu bytes\n", filename, data_len);
    // char* tmp = (char*)data;
    // for (size_t i = 0; i < data_len; ++i)
    // {
    //     LOGE("%3d\n", tmp[i]);
    // }

	file.write((char*)&header, sizeof(astc_header));
	file.write((char*)data, data_len);
	return 0;
}

/* ========================================================================= */

// Taken from astcenc:
// routine to write up to 8 bits
void write_bits(
	int value,
	int bitcount,
	int bitoffset,
	uint8_t* ptr
) {
	int mask = (1 << bitcount) - 1;
	value &= mask;
	ptr += bitoffset >> 3;
	bitoffset &= 7;
	value <<= bitoffset;
	mask <<= bitoffset;
	mask = ~mask;

	ptr[0] &= mask;
	ptr[0] |= value;
	ptr[1] &= mask >> 8;
	ptr[1] |= value >> 8;
}

/* Taken from astcenc */
int bitrev8(int p)
{
	p = ((p & 0xF) << 4) | ((p >> 4) & 0xF);
	p = ((p & 0x33) << 2) | ((p >> 2) & 0x33);
	p = ((p & 0x55) << 1) | ((p >> 1) & 0x55);
	return p;
}

} // namespace simple::astc
