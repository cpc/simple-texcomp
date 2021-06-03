#include <cstring>
#include <cassert>
#include <fstream>

#include "platform.hpp"
#include "simple_mathlib.hpp"
#include "simple_texcomp.hpp"

#include <Tracy.hpp>

namespace simple::astc {

/**
 * TODO: Optimize this for different quantizations
 */
constexpr decimal INSET_MARGIN = (F(8.0) / F(255.0)) / F(16.0);

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

// Quantization midpoints for better rounding
// taken from: https://gist.github.com/castano/c92c7626f288f9e99e158520b14a61cf
constexpr decimal quant_midpoints_2b[4] = {
    F(0.16666667), F(0.50000000), F(0.83333333), F(1.00000000)
};

constexpr decimal quant_midpoints_3b[8] = {
    F(0.07058824), F(0.21372549), F(0.35686275), F(0.50000000), F(0.64313725), F(0.78627451), F(0.92941176), F(1.00000000)
};

constexpr decimal quant_midpoints_5b[32] = {
    F(0.01568627), F(0.04705882), F(0.07843137), F(0.11176471), F(0.14509804), F(0.17647059), F(0.20784314), F(0.24117647),
    F(0.27450980), F(0.30588235), F(0.33725490), F(0.37058824), F(0.40392157), F(0.43529412), F(0.46666667), F(0.50000000),
    F(0.53333333), F(0.56470588), F(0.59607843), F(0.62941176), F(0.66274510), F(0.69411765), F(0.72549020), F(0.75882353),
    F(0.79215686), F(0.82352941), F(0.85490196), F(0.88823529), F(0.92156863), F(0.95294118), F(0.98431373), F(1.00000000)
};

// calculate the tables above
void init_tables() {
    for (int i = 0; i < 3; i++) {
        decimal f0 = decimal(((i+0) << 6) | ((i+0) << 4) | ((i+0) << 2) | (i+0)) / F(255.0);
        decimal f1 = decimal(((i+1) << 6) | ((i+1) << 4) | ((i+1) << 2) | (i+1)) / F(255.0);
        printf("%.8f, ", (double)((f0 + f1) * F(0.5)));
    }
    printf("%.8f\n\n", 1.0);

    for (int i = 0; i < 7; i++) {
        decimal f0 = decimal(((i+0) << 5) | ((i+0) << 2) | ((i+0) >> 1)) / F(255.0);
        decimal f1 = decimal(((i+1) << 5) | ((i+1) << 2) | ((i+1) >> 1)) / F(255.0);
        printf("%.8f, ", (double)((f0 + f1) * F(0.5)));
    }
    printf("%.8f\n\n", 1.0);

    for (int i = 0; i < 31; i++) {
        decimal f0 = decimal(((i+0) << 3) | ((i+0) >> 2)) / F(255.0);
        decimal f1 = decimal(((i+1) << 3) | ((i+1) >> 2)) / F(255.0);
        printf("%.8f, ", (double)((f0 + f1) * F(0.5)));
        if ((i % 8) == 7) {
            printf("\n");
        }
    }
    printf("%.8f\n", 1.0);
}

/* Pre-computed bilinear filter weights */
bilin::bilinear_weights bilin_weights;

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

void print_bin(unsigned int num, unsigned int nb)
{
	for (uint b = 0; b < nb; ++b)
	{
		unsigned int x = (num >> (nb-1-b)) & 1;
		printf("%d", x);
	}
}

void print_minmax(const char *pre, Vec3f mincol, Vec3f maxcol)
{
    Vec3f mincol2 = mincol * F(255.0);
    Vec3f maxcol2 = maxcol * F(255.0);

    printf("%s mincol: %5.3f %5.3f %5.3f  %3.0f %3.0f %3.0f\n",
        pre,
        (double)mincol.x, (double)mincol.y, (double)mincol.z,
        (double)mincol2.x, (double)mincol2.y, (double)mincol2.z
    );
    printf("    maxcol: %5.3f %5.3f %5.3f  %3.0f %3.0f %3.0f\n",
        (double)maxcol.x, (double)maxcol.y, (double)maxcol.z,
        (double)maxcol2.x, (double)maxcol2.y, (double)maxcol2.z
    );
}

/* Find min/max color as a corners of a bounding box of the block */
// static void find_minmaxcolor_bbox_astc(
//     const Vec3f* block,
//     int pixel_count,
//     Vec3f *mincol,
//     Vec3f *maxcol
// ){
//     *mincol = { F(1.0), F(1.0), F(1.0) };
//     *maxcol = { F(0.0), F(0.0), F(0.0) };

//     for (int i = 0; i < pixel_count; ++i)
//     {
//         assert(!std::isnan(block[i].x));
//         assert(!std::isnan(block[i].y));
//         assert(!std::isnan(block[i].z));
//         assert((block[i].x >= F(0.0)) && (block[i].x <= F(1.0)));
//         assert((block[i].y >= F(0.0)) && (block[i].y <= F(1.0)));
//         assert((block[i].z >= F(0.0)) && (block[i].z <= F(1.0)));

//         *mincol = min3f(*mincol, block[i]);
//         *maxcol = max3f(*maxcol, block[i]);
//     }
// }

/* Shrink the bounding box */
void inset_bbox(Vec3f *mincol, Vec3f *maxcol)
{
    Vec3f inset = (*maxcol - *mincol) * (F(1.0) / F(16.0)) - INSET_MARGIN;
    *mincol = clamp3f(*mincol + inset, F(0.0), F(1.0));
    *maxcol = clamp3f(*maxcol - inset, F(0.0), F(1.0));
}

/** Find min/max values but consider only values within 2.0*stddev
 */
void find_minmax_trimmed(
    const Vec3f *block,
    uint8_t pixel_count,
    const Vec3f &avg,
    const Vec3f &std,
    Vec3f *mincol,
    Vec3f *maxcol
){
    ZoneScopedN("minmax_tr");
    // Vec3f avg2 = avg * F(255.0);
    // Vec3f std2 = std * F(255.0);
    // printf("       avg: %5.3f %5.3f %5.3f  %3.0f %3.0f %3.0f\n",
    //     (double)avg.x, (double)avg.y, (double)avg.z,
    //     (double)avg2.x, (double)avg2.y, (double)avg2.z);
    // printf("       std: %5.3f %5.3f %5.3f  %3.0f %3.0f %3.0f\n",
    //     (double)std.x, (double)std.y, (double)std.z,
    //     (double)std2.x, (double)std2.y, (double)std2.z);

    // TODO rlcamp, lclamp
    Vec3f stdmax = clamp3f(avg + std*F(2.0), F(0.0), F(1.0));
    Vec3f stdmin = clamp3f(avg - std*F(2.0), F(0.0), F(1.0));
    *mincol = avg;
    *maxcol = avg;
    // print_minmax("s2 ", stdmin, stdmax);
    for (int i = 0; i < pixel_count; ++i)
    {
        if ((block[i].x < stdmax.x) && (block[i].x > stdmin.x))
        {
            maxcol->x = fmax(maxcol->x, block[i].x);
            mincol->x = fmin(mincol->x, block[i].x);
        }

        if ((block[i].y < stdmax.y) && (block[i].y > stdmin.y))
        {
            maxcol->y = fmax(maxcol->y, block[i].y);
            mincol->y = fmin(mincol->y, block[i].y);
        }

        if ((block[i].z < stdmax.z) && (block[i].z > stdmin.z))
        {
            maxcol->z = fmax(maxcol->z, block[i].z);
            mincol->z = fmin(mincol->z, block[i].z);
        }
    }
}

#if ASTC_SELECT_DIAG == 1
/* Optional selection of either current or oposite diagonal - small potential
 * quality improvement at a small runtime cost
 */
bool select_diagonal(
    const Vec3f *block,
    uint8_t pixel_count,
    Vec3f *mincol,
    Vec3f *maxcol
){
    ZoneScopedN("sel_diag");

    bool swapped = false;
    Vec3f center = (*mincol + *maxcol) * F(0.5);

    Vec2f cov = { F(0.0), F(0.0) };
    for (int i = 0; i < pixel_count; ++i) {
        Vec3f t = block[i] - center;
        cov.x += t.x * t.z;
        cov.y += t.y * t.z;
    }

    // printf("     cov_x: %8.5f\n", (double)cov.x);
    // printf("     cov_y: %8.5f\n", (double)cov.y);

    if (cov.x < F(0.0)) {
        decimal tmp = maxcol->x;
        maxcol->x = mincol->x;
        mincol->x = tmp;
        swapped = true;
    }

    if (cov.y < F(0.0)) {
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
uint8_t quantize_2b(decimal x)
{
    assert((x >= F(0.0)) && (x <= F(1.0)));
    assert(!std::isnan(x));

    // There is a couple of ways to do this. The correct one would be:
    //
    // uint8_t quant = (uint8_t)(x * F(3.0));
    // quant += (x > quant_midpoints_2b[quant]);
    //
    // or without ideal rounding (about 0.0016 dB loss):
    //
    // uint8_t quant = u8clamp((uint8_t)(x * F(3.0) + F(0.5)), 0, 3);
    //
    // However, I found that for some reason, quantizing first to 3 bits, then
    // bit-shifting to 2 bits yealds better quality.
    // The ideal rounding yields about .001 dB improvement, and only when it is
    // performed after the first rounding (see the +0.5), otherwise has no
    // effect.

    uint8_t quant = u8clamp((uint8_t)(x * F(7.0) + F(0.5)), 0, 7);
    // uint8_t quant = u8clamp((uint8_t)(x * F(7.0)), 0, 7);
    // quant += (x > quant_midpoints_3b[quant]);
    quant >>= 1;

    assert((quant >= 0) && (quant <= 3));
    return quant;
}

/** Quantize a decimal value into 5 bits (32 values)
 *
 * Returns the quantized input decimal as 8-bit integer and also modifies the
 * input decimal to the new value
 */
Vec3i quantize_5b(Vec3f *vec)
{
    assert(!std::isnan(vec->x));
    assert(!std::isnan(vec->y));
    assert(!std::isnan(vec->z));
    assert((vec->x >= F(0.0)) && (vec->x <= F(1.0)));
    assert((vec->y >= F(0.0)) && (vec->y <= F(1.0)));
    assert((vec->z >= F(0.0)) && (vec->z <= F(1.0)));

    // Ideal rounding yields about 0.002 dB improvement. Otherwise, non-ideal
    // rounding:
    int quant_x = iclamp((int)(vec->x * F(31.0) + F(0.5)), 0, 31);
    int quant_y = iclamp((int)(vec->y * F(31.0) + F(0.5)), 0, 31);
    int quant_z = iclamp((int)(vec->z * F(31.0) + F(0.5)), 0, 31);

    // int quant_x = (int)(vec->x * F(31.0));
    // int quant_y = (int)(vec->y * F(31.0));
    // int quant_z = (int)(vec->z * F(31.0));

    // quant_x += (vec->x > quant_midpoints_5b[quant_x]);
    // quant_y += (vec->y > quant_midpoints_5b[quant_y]);
    // quant_z += (vec->z > quant_midpoints_5b[quant_z]);

    int dequant_x = (quant_x << 3) | (quant_x >> 2);
    int dequant_y = (quant_y << 3) | (quant_y >> 2);
    int dequant_z = (quant_z << 3) | (quant_z >> 2);

    vec->x = (decimal)(dequant_x) * (F(1.0) / F(255.0));
    vec->y = (decimal)(dequant_y) * (F(1.0) / F(255.0));
    vec->z = (decimal)(dequant_z) * (F(1.0) / F(255.0));

    assert((vec->x >= F(0.0)) && (vec->x <= F(1.0)));
    assert((vec->y >= F(0.0)) && (vec->y <= F(1.0)));
    assert((vec->z >= F(0.0)) && (vec->z <= F(1.0)));
    return Vec3i { quant_x, quant_y, quant_z };
}

void encode_block(
    const uint8_t block_pixels[NCH_RGB*BLOCK_X*BLOCK_Y],
    uint32_t out[4]
){
    ZoneScopedN("enc_blk_astc");

    const uint16_t block_mode = 102;
    const uint8_t wgt_grid_w = 8;
    const uint8_t wgt_grid_h = 5;
    const uint8_t wgt_count = wgt_grid_w * wgt_grid_h;

    // const uint8_t color_quant_level = 11;  // range 32
    const uint8_t ep_bits = 5;
    // const uint8_t wgt_quant_mode = 2;      // range 4
    const uint8_t wgt_bits = 2;

    constexpr uint8_t block_size_x = 12;
    constexpr uint8_t block_size_y = 12;
    constexpr uint8_t pixel_count = block_size_x * block_size_y;

    constexpr uint8_t MAX_PIXEL_COUNT = MAX_BLOCK_DIM*MAX_BLOCK_DIM;

    // printf("=== BLOCK ===\n");

    // Convert the block into floating point and determine the line through
    // color space
    Vec3f block_flt[MAX_PIXEL_COUNT];
#if ASTC_TRIM_ENDPOINTS == 0
    // Default method, select min/max directly
    // using tmp values to allow vectorization
    decimal tmp_min[3] = { F(1.0), F(1.0), F(1.0) };
    decimal tmp_max[3] = { F(0.0), F(0.0), F(0.0) };
    {
        ZoneScopedN("minmax");
        for (int i = 0; i < pixel_count; ++i) {
            block_flt[i].x = (decimal) block_pixels[NCH_RGB * i] / F(255.0);
            block_flt[i].y = (decimal) block_pixels[NCH_RGB * i + 1] / F(255.0);
            block_flt[i].z = (decimal) block_pixels[NCH_RGB * i + 2] / F(255.0);

            tmp_min[0] = fmin(tmp_min[0], block_flt[i].x);
            tmp_min[1] = fmin(tmp_min[1], block_flt[i].y);
            tmp_min[2] = fmin(tmp_min[2], block_flt[i].z);
            tmp_max[0] = fmax(tmp_max[0], block_flt[i].x);
            tmp_max[1] = fmax(tmp_max[1], block_flt[i].y);
            tmp_max[2] = fmax(tmp_max[2], block_flt[i].z);
        }
    }
    Vec3f mincol = { tmp_min[0], tmp_min[1], tmp_min[2] };
    Vec3f maxcol = { tmp_max[0], tmp_max[1], tmp_max[2] };
    // print_minmax("   ", mincol, maxcol);
#else
    // Trimmed method, min/max are selected so they are within avg+=2.0*stddev
    Vec3f sum = { F(0.0), F(0.0), F(0.0) };
    Vec3f sq_sum = { F(0.0), F(0.0), F(0.0) };
    for (int i = 0; i < pixel_count; ++i)
    {
        block_flt[i].x = (decimal)block_pixels[NCH_RGB*i] / F(255.0);
        block_flt[i].y = (decimal)block_pixels[NCH_RGB*i+1] / F(255.0);
        block_flt[i].z = (decimal)block_pixels[NCH_RGB*i+2] / F(255.0);

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

    Vec3f mincol = avg;
    Vec3f maxcol = avg;
    find_minmax_trimmed(block_flt, pixel_count, avg, std, &mincol, &maxcol);
    // print_minmax("m2 ", mincol, maxcol);
#endif

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

    uint8_t quantized_weights[MAX_PIXEL_COUNT];

    if (mincol_int == maxcol_int)
    {
        // When both endpoints are the same, encode all weights and zeros and
        // skip the calculation
        for (int i = 0; i < wgt_count; ++i)
        {
            quantized_weights[i] = 0;
        }
    } else {
        // Downsample and quantize the weights
        decimal downsampled_weights[MAX_PIXEL_COUNT];

        // Move the endpoints line segment such that mincol is at zero
        Vec3f ep_vec = maxcol - mincol;
        // Normalize the endpoint vector
        // decimal norm = F(1.0) / std::sqrt(ep_vec.dot(ep_vec));
        // Vec3f ep_vec_norm = ep_vec * norm;

        // It works when the norm is squared, why?
        Vec3f ep_vec_scaled = ep_vec / ep_vec.dot(ep_vec);
        assert(ep_vec.dot(ep_vec) != F(0.0));

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
                F(0.0),
                F(1.0)
            );
        }

        // We downsample the weight grid before quantization
        // TODO: Investigate possible shift in values
        bilin::downsample(
            ideal_weights,
            block_size_x,
            block_size_y,
            &bilin_weights,
            downsampled_weights,
            wgt_grid_w,
            wgt_grid_h
        );

        // Quantize weights
        for (int i = 0; i < wgt_count; ++i)
        {
            quantized_weights[i] = quantize_2b(downsampled_weights[i]);
        }
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

    printf("Quantized weights:\n");
    for (int i = 0; i < wgt_count; ++i)
    {
        printf("%4d" , quantized_weights[i]);
        if ((i % wgt_grid_w) == (wgt_grid_w - 1))
        {
            printf("\n");
        }
    }

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

} // namespace
