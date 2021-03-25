#include <cstring>
#include <fstream>

#include "simple_texcomp.hpp"

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

static const uint32_t ASTC_MAGIC_ID = 0x5CA1AB13;

int store_astc_image(
	const uint8_t* data,
	const uint32_t data_len,
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

// Taken from astcenc:
// routine to write up to 8 bits
static inline void write_bits(
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
static inline int bitrev8(int p)
{
	p = ((p & 0xF) << 4) | ((p >> 4) & 0xF);
	p = ((p & 0x33) << 2) | ((p >> 2) & 0x33);
	p = ((p & 0x55) << 1) | ((p >> 1) & 0x55);
	return p;
}

static void print_bin(unsigned int num, unsigned int nb)
{
	for (uint b = 0; b < nb; ++b)
	{
		unsigned int x = (num >> (nb-1-b)) & 1;
		printf("%d", x);
	}
}

/** Downsample input block XxY into the size of MxN
 *
 * Returns 1 in case of error, 0 on success
 */
static int bilinear_downsample(
    const float* inp,  // input values
    int nch,           // number of channels (e.g. 3 for RGB or 1 for grayscale)
    int w_inp,         // width of the input block (X)
    int h_inp,         // height of the input block (Y)
    float* out,        // output values
    int w_out,         // width of the output block (M; M <= X)
    int h_out          // height of the output block (N; N <= Y)
){
    if ((w_out > w_inp) || (h_out > h_inp))
    {
        return 1;
    }

    float x = 0.5f;
    float step_x = 1.0f / (float)w_inp;
    float m = 0.5f;
    float step_m = 1.0f / (float)w_out;

    float nsteps = step_m / step_x;
    printf("number of taps: %6.3f\n", (double)nsteps);

    for (int i = 0; i < w_out; ++i)
    {


        m += step_m;
    }

    return 0;
}


void encode_block_astc(
    const uint8_t block_pixels[NCH_RGB*ASTC_BLOCK_X*ASTC_BLOCK_Y],
    uint32_t out[4]
){
    const uint16_t block_mode = 102;
    const uint8_t wgt_grid_w = 8;
    const uint8_t wgt_grid_h = 5;
    const uint8_t wgt_count = wgt_grid_w * wgt_grid_h;

    const uint8_t color_quant_level = 11;  // range 32
    const uint8_t ep_bits = 5;
    const uint8_t wgt_quant_mode = 2;      // range 4
    const uint8_t wgt_bits = 2;

    const int wgt_bitcount = wgt_count * wgt_bits;

    const int X = 12;
    const int Y = 12;
    const int WX = 8;
    const int WY = 5;
    const float wx_step = 1.0f / ((float)WX - 1);
    const float wy_step = 1.0f / ((float)WY - 1);

    int res = bilinear_downsample(
        NULL,
        1,
        X, Y,
        NULL,
        WX, WY
    );
    if (res) { printf("ERROR"); }

    // printf("texels: %dx%d,  weights: %dx%d\n", X, Y, WX, WY);
    // for (int j = 0; j < WY; ++j)
    // {
    //     const float wy = j * wy_step;
    //     const float y = wy * Y;
    //     printf("%4.2f(%5.2f):", wy, y);
    //     for (int i = 0; i < WX; ++i)
    //     {
    //         const float wx = i * wx_step;
    //         const float x = wx * X;
    //         printf("%6.2f(%5.2f)", wx, x);
    //     }
    //     printf("\n");
    // }

    exit(0);

    // Weights after decimation
    const uint8_t weights[wgt_count] = {  // wgt_count == 40
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 1, 0,
        0, 1, 2, 2, 3, 2, 2, 0,
        0, 1, 1, 1, 1, 1, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
    };

	uint8_t wgt_buf[16];
	uint8_t out_buf[16];
    for (int i = 0; i < 16; ++i)
    {
        wgt_buf[i] = 0;
        out_buf[i] = 0;
    }

    // ISE encoding
    int off = 0;
    for (int i = 0; i < wgt_count; ++i)
    {
        write_bits(weights[i], wgt_bits, off, wgt_buf);
        off += wgt_bits;
    }

    // write out weights
    for (int i = 0; i < 16; ++i)
    {
        out_buf[i] = bitrev8(wgt_buf[15-i]);
    }

    // write out mode, partition, CEM
    write_bits(block_mode, 11, 0, out_buf);
    const int partition_count = 1;
    write_bits(partition_count - 1, 2, 11, out_buf);
    const int color_format = 8;
    write_bits(color_format, 4, 13, out_buf);

    // Unquantized endpoints
    const uint8_t mincol[3] = { 0, 0, 0 };
    const uint8_t maxcol[3] = { 255, 255, 255 };
    if (mincol[0] + mincol[1] + mincol[2] > maxcol[0] + maxcol[1] + maxcol[2])
    {
        printf("ERROR: mincol must be smaller than maxcol\n");
    }

    // Quantized endpoints (layout is R0 R1 G0 G1 B0 B1)
    uint8_t endpoints_q[6] = { 0 };
    for (int i = 0; i < 3; ++i)
    {
        endpoints_q[2*i] = mincol[i] >> 3;
        endpoints_q[2*i+1] = maxcol[i] >> 3;
    }

    // write out endpoints
    off = 17;  // starting bit position for endpoints data
    for (int i = 0; i < 6; ++i)
    {
        printf("ep[%d]: %d\n", i, endpoints_q[i]);
        write_bits(endpoints_q[i], ep_bits, off, out_buf);
        off += ep_bits;
    }

    printf("Output block (msb..lsb):\n");
    for (int i = 0; i < 16; ++i)
    {
        printf("i: %2d  wgt: ", i);
        print_bin(wgt_buf[i], 8);
        printf("  out: ");
        print_bin(out_buf[i], 8);
        printf("\n");
    }

    // uint8_t encoded_block_bytes[16] = {
    //     209, 80, 60, 7, 82, 165, 230, 115, 164, 113, 51, 128, 8, 204, 179, 179
    // };

    memcpy(out, out_buf, 16);
}

// void decode_block_astc(
//     const uint32_t enc_block[4],
//     uint8_t out_pixels[NCH_RGB*16]
// ){
// }
