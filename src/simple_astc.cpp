#include <cstring>
#include <cassert>
#include <fstream>

#include "simple_mathlib.hpp"
#include "simple_texcomp.hpp"

/**
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

// Quantization midpoints for better rounding
// taken from: https://gist.github.com/castano/c92c7626f288f9e99e158520b14a61cf
static const decimal quant_midpoints_2b[4] = {
    0.16666667, 0.50000000, 0.83333333, 1.00000000
};

static const decimal quant_midpoints_3b[8] = {
    0.07058824, 0.21372549, 0.35686275, 0.50000000, 0.64313725, 0.78627451, 0.92941176, 1.00000000
};

static const decimal quant_midpoints_5b[32] = {
    0.01568627, 0.04705882, 0.07843137, 0.11176471, 0.14509804, 0.17647059, 0.20784314, 0.24117647,
    0.27450980, 0.30588235, 0.33725490, 0.37058824, 0.40392157, 0.43529412, 0.46666667, 0.50000000,
    0.53333333, 0.56470588, 0.59607843, 0.62941176, 0.66274510, 0.69411765, 0.72549020, 0.75882353,
    0.79215686, 0.82352941, 0.85490196, 0.88823529, 0.92156863, 0.95294118, 0.98431373, 1.00000000
};

// calculate the tables above
void init_tables() {
    for (int i = 0; i < 3; i++) {
        decimal f0 = decimal(((i+0) << 6) | ((i+0) << 4) | ((i+0) << 2) | (i+0)) / 255.0;
        decimal f1 = decimal(((i+1) << 6) | ((i+1) << 4) | ((i+1) << 2) | (i+1)) / 255.0;
        printf("%.8f, ", (f0 + f1) * 0.5);
    }
    printf("%.8f\n\n", 1.0);

    for (int i = 0; i < 7; i++) {
        decimal f0 = decimal(((i+0) << 5) | ((i+0) << 2) | ((i+0) >> 1)) / 255.0;
        decimal f1 = decimal(((i+1) << 5) | ((i+1) << 2) | ((i+1) >> 1)) / 255.0;
        printf("%.8f, ", (f0 + f1) * 0.5);
    }
    printf("%.8f\n\n", 1.0);

    for (int i = 0; i < 31; i++) {
        decimal f0 = decimal(((i+0) << 3) | ((i+0) >> 2)) / 255.0;
        decimal f1 = decimal(((i+1) << 3) | ((i+1) >> 2)) / 255.0;
        printf("%.8f, ", (f0 + f1) * 0.5);
        if ((i % 8) == 7) {
            printf("\n");
        }
    }
    printf("%.8f\n", 1.0);
}

/* Pre-computed bilinear filter weights */
static bilinear_weights bilin_weights;

int init_astc(
    int block_size_x,
    int block_size_y,
    int weight_grid_x,
    int weight_grid_y
){
    // init_tables();
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

static void print_minmax(const char *pre, Vec3f mincol, Vec3f maxcol)
{
    Vec3f mincol2 = mincol * 255.0;
    Vec3f maxcol2 = maxcol * 255.0;

    printf("%s mincol: %5.3f %5.3f %5.3f  %3.0f %3.0f %3.0f\n",
        pre,
        mincol.x, mincol.y, mincol.z,
        mincol2.x, mincol2.y, mincol2.z
    );
    printf("    maxcol: %5.3f %5.3f %5.3f  %3.0f %3.0f %3.0f\n",
        maxcol.x, maxcol.y, maxcol.z,
        maxcol2.x, maxcol2.y, maxcol2.z
    );
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
        assert(!std::isnan(block[i].x));
        assert(!std::isnan(block[i].y));
        assert(!std::isnan(block[i].z));
        assert((block[i].x >= 0.0) && (block[i].x <= 1.0));
        assert((block[i].y >= 0.0) && (block[i].y <= 1.0));
        assert((block[i].z >= 0.0) && (block[i].z <= 1.0));

        *mincol = min3f(*mincol, block[i]);
        *maxcol = max3f(*maxcol, block[i]);
    }
}

/* Shrink the bounding box */
static inline void inset_bbox(Vec3f *mincol, Vec3f *maxcol)
{
    Vec3f inset = (*maxcol - *mincol) * (1.0 / 16.0) - INSET_MARGIN;
    *mincol = clamp3f(*mincol + inset, 0.0, 1.0);
    *maxcol = clamp3f(*maxcol - inset, 0.0, 1.0);
}

#if ASTC_SELECT_DIAG == 1
/* Optional selection of either current or oposite diagonal - small potential
 * quality improvement at a small runtime cost
 */
static inline bool select_diagonal(
    const Vec3f *block,
    uint8_t pixel_count,
    Vec3f *mincol,
    Vec3f *maxcol
){
    bool swapped = false;
    Vec3f center = (*mincol + *maxcol) * 0.5;

    Vec2f cov = { 0.0, 0.0 };
    for (int i = 0; i < pixel_count; ++i) {
        Vec3f t = block[i] - center;
        cov.x += t.x * t.z;
        cov.y += t.y * t.z;
    }

    // printf("    cov_x: %.5f  cov_y: %.5f\n", cov.x, cov.y);

    if (cov.x < 0.0) {
        decimal tmp = maxcol->x;
        maxcol->x = mincol->x;
        mincol->x = tmp;
        swapped = true;
    }

    if (cov.y < 0.0) {
        decimal tmp = maxcol->y;
        maxcol->y = mincol->y;
        mincol->y = tmp;
        swapped = true;
    }

    return swapped;
}
#endif

/* TODO: Fix this */
/** Quantize decimal value into 2 bits (4 values)
 *
 * The input decimal does not need to be updated in our use case
 */
static inline uint8_t quantize_2b(decimal x)
{
    assert((x >= 0.0) && (x <= 1.0));
    assert(!std::isnan(x));

    // There is a couple of ways to do this. The correct one would be:
    //
    // uint8_t quant = (uint8_t)(x * 3.0);
    // quant += (x > quant_midpoints_2b[quant]);
    //
    // or without ideal rounding (about 0.0016 dB loss):
    //
    // uint8_t quant = u8clamp((uint8_t)(x * 3.0 + 0.5), 0, 3);
    //
    // However, I found that for some reason, quantizing first to 3 bits, then
    // bit-shifting to 2 bits yealds better quality.
    // The ideal rounding yields about .001 dB improvement, and only when it is
    // performed after the first rounding (see the +0.5), otherwise has no
    // effect.

    uint8_t quant = u8clamp((uint8_t)(x * 7.0 + 0.5), 0, 7);
    quant += (x > quant_midpoints_3b[quant]);
    quant >>= 1;

    assert((quant >= 0) && (quant <= 3));
    return quant;
}

/** Quantize a decimal value into 5 bits (32 values)
 *
 * Returns the quantized input decimal as 8-bit integer and also modifies the
 * input decimal to the new value
 */
static inline Vec3i quantize_5b(Vec3f *vec)
{
    assert(!std::isnan(vec->x));
    assert(!std::isnan(vec->y));
    assert(!std::isnan(vec->z));
    assert((vec->x >= 0.0) && (vec->x <= 1.0));
    assert((vec->y >= 0.0) && (vec->y <= 1.0));
    assert((vec->z >= 0.0) && (vec->z <= 1.0));

    // Ideal rounding yields about 0.002 dB improvement. Otherwise, non-ideal
    // rounding:
    // int quant_x = iclamp((int)(vec->x * 31.0 + 0.5), 0, 31);

    int quant_x = (int)(vec->x * 31.0);
    int quant_y = (int)(vec->y * 31.0);
    int quant_z = (int)(vec->z * 31.0);

    quant_x += (vec->x > quant_midpoints_5b[quant_x]);
    quant_y += (vec->y > quant_midpoints_5b[quant_y]);
    quant_z += (vec->z > quant_midpoints_5b[quant_z]);

    int dequant_x = (quant_x << 3) | (quant_x >> 2);
    int dequant_y = (quant_y << 3) | (quant_y >> 2);
    int dequant_z = (quant_z << 3) | (quant_z >> 2);

    vec->x = (decimal)(dequant_x) * (1.0 / 255.0);
    vec->y = (decimal)(dequant_y) * (1.0 / 255.0);
    vec->z = (decimal)(dequant_z) * (1.0 / 255.0);

    assert((vec->x >= 0.0) && (vec->x <= 1.0));
    assert((vec->y >= 0.0) && (vec->y <= 1.0));
    assert((vec->z >= 0.0) && (vec->z <= 1.0));
    return Vec3i { quant_x, quant_y, quant_z };
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
    Vec3f sum = { 0.0, 0.0, 0.0 };
    Vec3f sq_sum = { 0.0, 0.0, 0.0 };
    // TODO: this loop can be merged with find_minmaxcolor_bbox_astc()
    for (int i = 0; i < pixel_count; ++i)
    {
        block_flt[i].x = (decimal)block_pixels[NCH_RGB*i] / 255.0;
        block_flt[i].y = (decimal)block_pixels[NCH_RGB*i+1] / 255.0;
        block_flt[i].z = (decimal)block_pixels[NCH_RGB*i+2] / 255.0;

        sum = sum + block_flt[i];
        Vec3f sq = {
            block_flt[i].x * block_flt[i].x,
            block_flt[i].y * block_flt[i].y,
            block_flt[i].z * block_flt[i].z,
        };
        sq_sum = sq_sum + sq;
    }
    Vec3f avg = sum / (decimal)(pixel_count);
    Vec3f var =
        (sq_sum / pixel_count) - Vec3f { avg.x*avg.x, avg.y*avg.y, avg.z*avg.z };
    Vec3f std = {
        std::sqrt(var.x),
        std::sqrt(var.y),
        std::sqrt(var.z),
    };

    // printf("\n");
    // printf("avg: %5.3f %5.3f %5.3f  std: %5.3f %5.3f %5.3f\n",
    //     avg.x, avg.y, avg.z, std.x, std.y, std.z);

    // Determine line through color space
    Vec3f mincol, maxcol;
    find_minmaxcolor_bbox_astc(block_flt, pixel_count, &mincol, &maxcol);
    // print_minmax("   ", mincol, maxcol);
#if ASTC_SELECT_DIAG == 1
    bool swapped = select_diagonal(block_flt, pixel_count, &mincol, &maxcol);
    // if (swapped)
    // {
    //     print_minmax("swp", mincol, maxcol);
    // }
#endif
    inset_bbox(&mincol, &maxcol);
    // print_minmax("ins", mincol, maxcol);

    // Quantize endpoints
    Vec3i mincol_int = quantize_5b(&mincol);
    Vec3i maxcol_int = quantize_5b(&maxcol);
    // print_minmax("qnt", mincol, maxcol);
    // uint8_t maxcol_int[3] = {
    //     quantize_5b(&maxcol.x),
    //     quantize_5b(&maxcol.y),
    //     quantize_5b(&maxcol.z)
    // };

#if ASTC_SELECT_DIAG == 1
    if ( (mincol_int.x + mincol_int.y + mincol_int.z)
        > (maxcol_int.x + maxcol_int.y + maxcol_int.z) )
    {
        Vec3i tmpi = mincol_int;
        mincol_int = maxcol_int;
        maxcol_int = tmpi;
        Vec3f tmpf = mincol;
        mincol = maxcol;
        maxcol = tmpf;
        // print_minmax("swp", mincol, maxcol);
    }
#endif
    // printf("int mincol:                    %3d %3d %3d\n",
    //     mincol_int.x, mincol_int.y, mincol_int.z
    // );
    // printf("int maxcol:                    %3d %3d %3d\n",
    //     maxcol_int.x, maxcol_int.y, maxcol_int.z
    // );

    decimal downsampled_weights[MAX_PIXEL_COUNT] = { 0.0 };

    // TODO: use void extent to handle constant block
    if (mincol_int != maxcol_int)
    {
        // Move the endpoints line segment such that mincol is at zero
        Vec3f ep_vec = maxcol - mincol;
        // Normalize the endpoint vector
        // decimal norm = 1.0 / std::sqrt(ep_vec.dot(ep_vec));
        // Vec3f ep_vec_norm = ep_vec * norm;

        // It works when the norm is squared, why?
        Vec3f ep_vec_scaled = ep_vec / ep_vec.dot(ep_vec);
        assert(ep_vec.dot(ep_vec) != 0.0);

        // Project all pixels onto the endpoint vector. For each pixel, the result
        // tells how far it goes into the endpoint vector direction. Small values
        // (-> 0.0) mean closer to mincol, large values (-> 1.0) mean closer to
        // maxcol.
        // In other words, the resulting array is the array of ideal weights,
        // assuming there is no quantization.
        decimal ideal_weights[MAX_PIXEL_COUNT];
        for (int i = 0; i < pixel_count; ++i)
        {
            ideal_weights[i] = fclamp(
                (block_flt[i] - mincol).dot(ep_vec_scaled),
                0.0,
                1.0
            );
        }

        // We downsample the weight grid before quantization
        bilinear_downsample(
            ideal_weights,
            block_size_x,
            block_size_y,
            &bilin_weights,
            downsampled_weights,
            wgt_grid_w,
            wgt_grid_h
        );
    }

    // Quantize weights
    uint8_t quantized_weights[MAX_PIXEL_COUNT];
    for (int i = 0; i < wgt_count; ++i)
    {
        quantized_weights[i] = quantize_2b(downsampled_weights[i]);
    }

    // Output buffers for quantized weights and output data
	uint8_t wgt_buf[16];
	uint8_t out_buf[16];
    for (int i = 0; i < 16; ++i)
    {
        wgt_buf[i] = 0;
        out_buf[i] = 0;
    }

    // weights ISE encoding
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

    // Quantized endpoint output data (layout is R0 R1 G0 G1 B0 B1)
    uint8_t endpoints_q[6] = {
        (uint8_t)(mincol_int.x),
        (uint8_t)(maxcol_int.x),
        (uint8_t)(mincol_int.y),
        (uint8_t)(maxcol_int.y),
        (uint8_t)(mincol_int.z),
        (uint8_t)(maxcol_int.z),
    };
    // for (int i = 0; i < 3; ++i)
    // {
    //     endpoints_q[2*i] = mincol_int[i] >> 3;
    //     endpoints_q[2*i+1] = maxcol_int[i] >> 3;
    // }

    // write out endpoints
    off = 17;  // starting bit position for endpoints data
    for (int i = 0; i < 6; ++i)
    {
        write_bits(endpoints_q[i], ep_bits, off, out_buf);
        off += ep_bits;
    }

    memcpy(out, out_buf, 16);

    // DEBUG prints
#if 0
    printf("Pixels R:\n");
    for (int i = 0; i < pixel_count; ++i)
    {
        printf("%6.3f" , (double)block_flt[i].x);
        if ((i % block_size_x) == (block_size_x - 1))
        {
            printf("\n");
        }
    }

    printf("ep_vec: %5.3f %5.3f %5.3f  ep_vec_scaled:  %5.3f %5.3f %5.3f\n",
        (double)ep_vec.x, (double)ep_vec.y, (double)ep_vec.z,
        (double)ep_vec_scaled.x, (double)ep_vec_scaled.y, (double)ep_vec_scaled.z);

    printf("Weights:\n");
    for (int i = 0; i < pixel_count; ++i)
    {
        printf("%6.3f" , (double)ideal_weights[i]);
        if ((i % block_size_x) == (block_size_x - 1))
        {
            printf("\n");
        }
    }

    printf("Decimated weights:\n");
    for (int i = 0; i < wgt_count; ++i)
    {
        printf("%6.3f" , (double)downsampled_weights[i]);
        if ((i % wgt_grid_w) == (wgt_grid_w - 1))
        {
            printf("\n");
        }
    }
// #endif
//     printf("Quantized weights:\n");
//     for (int i = 0; i < wgt_count; ++i)
//     {
//         printf("%4d" , quantized_weights[i]);
//         if ((i % wgt_grid_w) == (wgt_grid_w - 1))
//         {
//             printf("\n");
//         }
//     }
// #if 0
    if ( (mincol_int[0] + mincol_int[1] + mincol_int[2])
        > (maxcol_int[0] + maxcol_int[1] + maxcol_int[2]) )
    {
        printf("ERROR: mincol must be smaller than maxcol\n");
    }

    printf("Endpoints:\n");
    for (int i = 0; i < 6; ++i)
    {
        printf("ep[%d]: %d\n", i, endpoints_q[i]);
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
#endif
}

// void decode_block_astc(
//     const uint32_t enc_block[4],
//     uint8_t out_pixels[NCH_RGB*16]
// ){
// }
