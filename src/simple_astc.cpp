#include <cstring>
#include <fstream>

#include "simple_mathlib.hpp"
#include "simple_texcomp.hpp"

/* Rounding the bounding box inset outwards (= (8.0/255.0)/16.0)
 *
 * TODO: Optimize this for different quantizations
 */
#define INSET_MARGIN  (8.0 / 255.0) / 16.0

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

/* Pre-computed bilinear filter weights */
static bilinear_weights bilin_weights;

int init_astc(
    int block_size_x,
    int block_size_y,
    int weight_grid_x,
    int weight_grid_y
){
    int ret = populate_bilinear_weights(
        block_size_x,
        block_size_y,
        &bilin_weights,
        weight_grid_x,
        weight_grid_y
    );
    return ret;
}

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

/* Find min/max color as a corners of a bounding box of the block */
static void find_minmaxcolor_bbox_astc(
    const Vec3f* block,
    int pixel_count,
    Vec3f *mincol,
    Vec3f *maxcol
){
    *mincol = { 1.0, 1.0, 1.0 };
    *maxcol = { 0.0, 0.0, 0.0 };

    for (int i = 0; i < pixel_count; ++i)
    {
        *mincol = min3f(*mincol, block[i]);
        *maxcol = max3f(*maxcol, block[i]);
    }
}

/* Shrink the bounding box (by half the distance of equidistant points) to
 * eliminate influence of outliers */
static inline void inset_bbox(Vec3f *mincol, Vec3f *maxcol)
{
    Vec3f inset = (*maxcol - *mincol) * (1.0 / 16.0) - INSET_MARGIN;
    *mincol = clamp3f(*mincol + inset, 0.0, 1.0);
    *maxcol = clamp3f(*maxcol - inset, 0.0, 1.0);
}

/** Quantize decimal value into 2 bits (4 values) */
static inline uint8_t quantize_2b(decimal x)
{
    uint8_t y = (uint8_t)round(x * 15.0);
    return (y << 4) | (y >> 1);
}

/** Quantize decimal value into 5 bits (32 values) */
static inline uint8_t quantize_5b(decimal x)
{
    uint8_t y = (uint8_t)round(x * 31.0);
    return (y << 3) | (y >> 2);
}

void encode_block_astc(
    const uint8_t block_pixels[NCH_RGB*ASTC_BLOCK_X*ASTC_BLOCK_Y],
    uint32_t out[4]
){
    const uint16_t block_mode = 102;
    const uint8_t wgt_grid_w = 8;
    const uint8_t wgt_grid_h = 5;
    const uint8_t wgt_count = wgt_grid_w * wgt_grid_h;

    // const uint8_t color_quant_level = 11;  // range 32
    const uint8_t ep_bits = 5;
    // const uint8_t wgt_quant_mode = 2;      // range 4
    const uint8_t wgt_bits = 2;

    const uint8_t block_size_x = 12;
    const uint8_t block_size_y = 12;
    const uint8_t pixel_count = block_size_x * block_size_y;

    constexpr uint8_t MAX_PIXEL_COUNT = ASTC_MAX_BLOCK_DIM*ASTC_MAX_BLOCK_DIM;

    // Convert the block into floating point
    Vec3f block_flt[MAX_PIXEL_COUNT];
    printf("Pixels R:\n");
    for (int i = 0; i < pixel_count; ++i)
    {
        block_flt[i].x = (decimal)block_pixels[NCH_RGB*i] / 255.0;
        block_flt[i].y = (decimal)block_pixels[NCH_RGB*i+1] / 255.0;
        block_flt[i].z = (decimal)block_pixels[NCH_RGB*i+2] / 255.0;
        printf("%6.3f" , (double)block_flt[i].x);
        if ((i % block_size_x) == (block_size_x - 1))
        {
            printf("\n");
        }
    }

    // Determine line through color space
    Vec3f mincol, maxcol;
    find_minmaxcolor_bbox_astc(block_flt, pixel_count, &mincol, &maxcol);
    inset_bbox(&mincol, &maxcol);

    printf("mincol: %5.3f %5.3f %5.3f  maxcol:  %5.3f %5.3f %5.3f\n",
        (double)mincol.x, (double)mincol.y, (double)mincol.z,
        (double)maxcol.x, (double)maxcol.y, (double)maxcol.z);

    // Move the endpoints line segments such that mincol is at zero
    Vec3f ep_vec = maxcol - mincol;
    // Normalize the endpoint vector
    // decimal norm = 1.0 / std::sqrt(ep_vec.dot(ep_vec));
    // Vec3f ep_vec_norm = ep_vec * norm;

    // It works when the norm is squared, why?
    Vec3f ep_vec_scaled = ep_vec / ep_vec.dot(ep_vec);

    printf("ep_vec: %5.3f %5.3f %5.3f  ep_vec_scaled:  %5.3f %5.3f %5.3f\n",
        (double)ep_vec.x, (double)ep_vec.y, (double)ep_vec.z,
        (double)ep_vec_scaled.x, (double)ep_vec_scaled.y, (double)ep_vec_scaled.z);
    // printf("ep_vec_len: %5.3f\n", (1.0 / norm));

    // Project all pixels onto the endpoint vector. For each pixel, the result
    // tells how far it goes into the endpoint vector direction. Small values
    // (-> 0.0) mean closer to mincol, large values (-> 1.0) mean closer to
    // maxcol.
    // In other words, the resulting array is the array of ideal weights,
    // assuming there is no quantization.
    printf("Weights:\n");
    decimal ideal_weights[MAX_PIXEL_COUNT];
    for (int i = 0; i < pixel_count; ++i)
    {
        ideal_weights[i] = fclamp(
            (block_flt[i] - mincol).dot(ep_vec_scaled), 0.0, 1.0
        );
        printf("%6.3f" , (double)ideal_weights[i]);
        if ((i % block_size_x) == (block_size_x - 1))
        {
            printf("\n");
        }
    }

    decimal downsampled_weights[MAX_PIXEL_COUNT];
    bilinear_downsample(
        ideal_weights,
        block_size_x,
        block_size_y,
        &bilin_weights,
        downsampled_weights,
        wgt_grid_w,
        wgt_grid_h
    );

    printf("Decimated weights:\n");
    for (int i = 0; i < wgt_count; ++i)
    {
        printf("%6.3f" , (double)downsampled_weights[i]);
        if ((i % wgt_grid_w) == (wgt_grid_w - 1))
        {
            printf("\n");
        }
    }

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

    // Weights after decimation
    // const uint8_t weights[wgt_count] = {  // wgt_count == 40
    //     0, 0, 0, 0, 0, 0, 0, 0,
    //     0, 1, 1, 1, 1, 1, 1, 0,
    //     0, 1, 2, 2, 3, 2, 2, 0,
    //     0, 1, 1, 1, 1, 1, 1, 0,
    //     0, 0, 0, 0, 0, 0, 0, 0,
    // };

    uint8_t quantized_weights[MAX_PIXEL_COUNT];
    printf("Quantized weights:\n");
    for (int i = 0; i < wgt_count; ++i)
    {
        quantized_weights[i] = quantize_2b(downsampled_weights[i]) >> 6;
        printf("%4d" , quantized_weights[i]);
        if ((i % wgt_grid_w) == (wgt_grid_w - 1))
        {
            printf("\n");
        }
    }

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
        write_bits(quantized_weights[i], wgt_bits, off, wgt_buf);
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

    // TODO: Move ep quant before weights calculation
    // Quantized endpoints
    const uint8_t mincol_int[3] = {
        quantize_5b(mincol.x),
        quantize_5b(mincol.y),
        quantize_5b(mincol.z)
    };
    const uint8_t maxcol_int[3] = {
        quantize_5b(maxcol.x),
        quantize_5b(maxcol.y),
        quantize_5b(maxcol.z)
    };
    if ( (mincol_int[0] + mincol_int[1] + mincol_int[2])
        > (maxcol_int[0] + maxcol_int[1] + maxcol_int[2]) )
    {
        printf("ERROR: mincol must be smaller than maxcol\n");
    }

    // Quantized endpoints (layout is R0 R1 G0 G1 B0 B1)
    uint8_t endpoints_q[6] = { 0 };
    for (int i = 0; i < 3; ++i)
    {
        endpoints_q[2*i] = mincol_int[i] >> 3;
        endpoints_q[2*i+1] = maxcol_int[i] >> 3;
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
