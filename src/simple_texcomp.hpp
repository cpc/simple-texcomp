#ifndef SIMPLE_TEXCOMP_HPP
#define SIMPLE_TEXCOMP_HPP

#include <cstdint>

namespace simple {

/* Compile-time switch between float and double for the calculations
 *
 * EPSILON is used for comparing floating point numbers
 */
#if FLOAT_PRECISION == 64
#define F(x) (x)
typedef double decimal;
constexpr decimal EPSILON = 1e-12;
#elif FLOAT_PRECISION == 16
#define F(x) (x##f16)
typedef _Float16 decimal;
constexpr decimal EPSILON = 1e-6f16;
#else
#define F(x) (x##f)
typedef float decimal;
constexpr decimal EPSILON = 1e-6f;
#endif

/* Number of channels if a RGB pixel. Used for calculating data sizes */
constexpr int NCH_RGB  = 4;

/******************************************************************************
 * ASTC
 *
 * Decoding function is purposefully missing. Instead, encoded images are saved
 * as .astc files and can be decoded with the reference ARM decoder.
 *****************************************************************************/

/* Optional refinement (can improve quality at small runtime cost)
 * 1 - enable, 0 - disable
 */
#ifndef ASTC_SELECT_DIAG
#define ASTC_SELECT_DIAG  1
#endif

/* Advanced endpoint selection method that minimizes the effect of outliers
 * 1 - enable, 0 - disable
 */
#ifndef ASTC_TRIM_ENDPOINTS
#define ASTC_TRIM_ENDPOINTS  0
#endif

namespace astc {

    /* Input ASTC block size */
    constexpr int BLOCK_X = 12;
    constexpr int BLOCK_Y = 12;

    /* Maximum input block size */
    constexpr int MAX_BLOCK_DIM = 12;
    /* Maximum dimension of the weight grid stored in an ASTC block
     *
     * 12x12 weights do not fit into a block but there could be e.g. 12x4
     */
    constexpr int MAX_GRID_DIM = 12;

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
    // int init_astc(
    //     int block_size_x,
    //     int block_size_y,
    //     int weight_grid_x,
    //     int weight_grid_y
    // );

    /* Encode a block of pixels into the ASTC format */
    void encode_block(
        const uint8_t block_pixels[NCH_RGB*BLOCK_X*BLOCK_Y],
        uint32_t out[4]
    );

    /* Encode a block of pixels into the ASTC format using integer arithmetic */
    void encode_block_int(
        const uint8_t block_pixels[NCH_RGB*BLOCK_X*BLOCK_Y],
        uint32_t out[4]
    );

}

/******************************************************************************
 * BC1
 *****************************************************************************/

/* Optional refinement (can improve quality at small runtime cost)
 * 1 - enable, 0 - disable
 */
#ifndef BC1_SELECT_DIAG
#define BC1_SELECT_DIAG  1
#endif

namespace bc1 {

    /* Encode a block of 4x4 pixels into the BC1 format */
    void encode_block(
        const uint8_t block_pixels[NCH_RGB*16],
        uint32_t out[2]
    );

    /* Decode an encoded block into an array of 16 pixels */
    void decode_block(
        const uint32_t enc_block[2],
        uint8_t out_pixels[NCH_RGB*16]
    );

}

/******************************************************************************
 * YCoCg-BC3
 *****************************************************************************/

namespace ycocg_bc3 {

    /* Encode a block of 4x4 pixels into the YCoCg-BC3 format */
    void encode_block(
        const uint8_t block_pixels[NCH_RGB*16],
        uint32_t out[4]
    );

    /* Decode an encoded block into an array of 16 pixels */
    void decode_block(
        const uint32_t enc_block[4],
        uint8_t out_pixels[NCH_RGB*16]
    );

}

/******************************************************************************
 * Bilinear interpolation
 *****************************************************************************/

namespace bilin {

    /* Data structure for keeping pre-computed bilinear weights.
     *
     * They are stored separately for X and Y dimension (requires M+N instead of
     * M*N memory).
     */
    struct bilinear_weights {
        // How many pixels go to a calculation of each weight
        uint8_t bilin_pixel_count_x[astc::MAX_GRID_DIM];
        uint8_t bilin_pixel_count_y[astc::MAX_GRID_DIM];
        // Each bilin_idx_x[i] stores an array of indices that point to locations
        // of pixels used to calculate that weight
        uint8_t bilin_idx_x[astc::MAX_GRID_DIM][astc::MAX_BLOCK_DIM];
        uint8_t bilin_idx_y[astc::MAX_GRID_DIM][astc::MAX_BLOCK_DIM];
        // Each bilin_weights_x[i] stores the actual weight values the previous
        // array points at
        decimal bilin_weights_x[astc::MAX_GRID_DIM][astc::MAX_BLOCK_DIM];
        decimal bilin_weights_y[astc::MAX_GRID_DIM][astc::MAX_BLOCK_DIM];
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
    void downsample(
        const decimal *__restrict__ inp, // input values
        int w_inp,                       // width of the input block (X)
        int h_inp,                       // height of the input block (Y)
        const bilinear_weights *__restrict__ bw,  // pre-computed values
        decimal *__restrict__ out,       // output values
        int w_out,                       // width of the output block (M)
        int h_out                        // height of the output block (N)
    );

    /** Same as downsample but with fixed inp/out size */
    void downsample_12x12_to_8x5(
        const decimal *__restrict__ inp,  // input values
        decimal *__restrict__ out         // output values
    );

    /** Same as downsample but with fixed inp/out size */
    void downsample_12x12_to_8x5_u8(
        const uint8_t *__restrict__ inp,  // input values
        uint8_t *__restrict__ out         // output values
    );

}

} // namespace simple

#endif // SIMPLE_TEXCOMP_HPP
