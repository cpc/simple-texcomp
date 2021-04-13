#ifndef SIMPLE_TEXCOMP_HPP
#define SIMPLE_TEXCOMP_HPP

#include <cstdint>

/* Number of channels if a RGB pixel. Used for calculating data sizes */
#define NCH_RGB  3

/* Compile-time switch between float and double for the calculations
 *
 * EPSILON is used for comparing floating point numbers
 */
#if USE_DOUBLE == 1
#define F(x) (x)     // do not append f to decimal literals
typedef double decimal;
const decimal EPSILON = 1e-12;
#else
#define F(x) (x##f)  // append f to decimal literals
typedef float decimal;
const decimal EPSILON = 1e-6f;
#endif

/******************************************************************************
 * ASTC
 *
 * Decoding function is purposefully missing. Instead, encoded images are saved
 * as .astc files and can be decoded with the reference ARM decoder.
 *****************************************************************************/

/* Input ASTC block size */
#define ASTC_BLOCK_X  12
#define ASTC_BLOCK_Y  12

/* Maximum input block size */
#define ASTC_MAX_BLOCK_DIM 12
/* Maximum dimension of the weight grid stored in an ASTC block
 *
 * 12x12 weights do not fit into a block but there could be e.g. 12x4
 */
#define ASTC_MAX_GRID_DIM  12

/* Optional refinement (can improve quality at small runtime cost)
 * 1 - enable, 0 - disable
 */
#ifndef ASTC_SELECT_DIAG
#define ASTC_SELECT_DIAG  1
#endif

/* Stuff for seamless switching between float and double
 */
#if USE_DOUBLE == 1
typedef double decimal;
#else
typedef float decimal;
#endif

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

/* Initialize ASTC-specific stuff */
int init_astc(
    int block_size_x,
    int block_size_y,
    int weight_grid_x,
    int weight_grid_y
);

/* Encode a block of pixels into the ASTC format */
void encode_block_astc(
    const uint8_t block_pixels[NCH_RGB*ASTC_BLOCK_X*ASTC_BLOCK_Y],
    uint32_t out[4]
);


/******************************************************************************
 * BC1
 *****************************************************************************/

/* Optional refinement (can improve quality at small runtime cost)
 * 1 - enable, 0 - disable
 */
#ifndef BC1_SELECT_DIAG
#define BC1_SELECT_DIAG  1
#endif

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


/******************************************************************************
 * Bilinear interpolation
 *****************************************************************************/

/* Data structure for keeping pre-computed bilinear weights.
 *
 * They are stored separately for X and Y dimension (requires M+N instead of
 * M*N memory).
 */
struct bilinear_weights {
    // How many pixels go to a calculation of each weight
    uint8_t bilin_pixel_count_x[ASTC_MAX_GRID_DIM];
    uint8_t bilin_pixel_count_y[ASTC_MAX_GRID_DIM];
    // Each bilin_idx_x[i] stores an array of indices that point to locations
    // of pixels used to calculate that weight
    uint8_t bilin_idx_x[ASTC_MAX_GRID_DIM][ASTC_MAX_BLOCK_DIM];
    uint8_t bilin_idx_y[ASTC_MAX_GRID_DIM][ASTC_MAX_BLOCK_DIM];
    // Each bilin_weights_x[i] stores the actual weight values the previous
    // array points at
    decimal bilin_weights_x[ASTC_MAX_GRID_DIM][ASTC_MAX_BLOCK_DIM];
    decimal bilin_weights_y[ASTC_MAX_GRID_DIM][ASTC_MAX_BLOCK_DIM];
};

/** Pre-compute bilinear interpolation weights
 */
int populate_bilinear_weights(
    int w_inp,             // input block width
    int h_inp,             // input block height
    bilinear_weights* bw,  // output data structure
    int w_out,             // interpolated block width
    int h_out              // interpolated block height
);

/** Downsample input block XxY into the size of MxN
 *
 * Assumes only 1 channel
 * Returns 1 in case of error, 0 on success
 */
void bilinear_downsample(
    const decimal* inp,  // input values
    int w_inp,           // width of the input block (X)
    int h_inp,           // height of the input block (Y)
    const bilinear_weights* bw, // table with pre-computed filter weights
    decimal* out,        // output values
    int w_out,           // width of the output block (M; M <= X)
    int h_out            // height of the output block (N; N <= Y)
);

#endif // SIMPLE_TEXCOMP_HPP
