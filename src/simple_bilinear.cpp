#include <cstdio>
#include <cmath>

#include "simple_texcomp.hpp"
#include "simple_mathlib.hpp"

#include <Tracy.hpp>

namespace simple::bilin {

/* Used only in downsample_ref
// Minimum representable positive floating point number
#if FLOAT_PRECISION < 32
// _Float16 is not in std
static constexpr decimal FLT_MINVAL = 0.000000059604645f16;
#else
static constexpr decimal FLT_MINVAL = std::numeric_limits<decimal>::min();
#endif // FLAOT_PRECISION < 32
*/

// Full table for 8x5, used just for reference
/*
static constexpr bilinear_weights BILIN_WEIGHTS_8x5 = {
    { 2,  3,  3,  3,  3,  3,  3,  2,  0,  0,  0,  0, }, // bilin_pixel_count_x
    { 3,  5,  6,  5,  3,  0,  0,  0,  0,  0,  0,  0, }, // bilin_pixel_count_y
    {
        {  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  1,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  4,  5,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  5,  6,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  7,  8,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  8,  9, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        { 10, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
    }, // bilin_idx_x
    {
        {  0,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  1,  2,  3,  4,  5,  0,  0,  0,  0,  0,  0,  0, },
        {  3,  4,  5,  6,  7,  8,  0,  0,  0,  0,  0,  0, },
        {  6,  7,  8,  9, 10,  0,  0,  0,  0,  0,  0,  0, },
        {  9, 10, 11,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
        {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, },
    }, // bilin_idx_y
    {
        {
            1.00000, 0.36364, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.63636, 0.72727, 0.09091, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.27273, 0.90909, 0.45455, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.54545, 0.81818, 0.18182, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.18182, 0.81818, 0.54545, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.45455, 0.90909, 0.27273, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.09091, 0.72727, 0.63636, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.36364, 1.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
    }, // bilin_weights_x
    {
        {
            1.00000, 0.63636, 0.27273, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.36364, 0.72727, 0.90909, 0.54545, 0.18182, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.09091, 0.45455, 0.81818, 0.81818, 0.45455, 0.09091,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.18182, 0.54545, 0.90909, 0.72727, 0.36364, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.27273, 0.63636, 1.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
        {
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        },
    }, // bilin_weights_y
};
*/

static constexpr uint8_t BILIN_IDX_X_8[8][3] = {
    {  0,  1,  0, },
    {  1,  2,  3, },
    {  2,  3,  4, },
    {  4,  5,  6, },
    {  5,  6,  7, },
    {  7,  8,  9, },
    {  8,  9, 10, },
    { 10, 11,  0, },
};

static constexpr uint8_t BILIN_IDX_Y_5[5][6] = {
    {  0,  1,  2,  0,  0,  0, },
    {  1,  2,  3,  4,  5,  0, },
    {  3,  4,  5,  6,  7,  8, },
    {  6,  7,  8,  9, 10,  0, },
    {  9, 10, 11,  0,  0,  0, },
};

static constexpr decimal BILIN_WEIGHTS_X_8[3*8] = {
    // not normalized
    // 1.00000, 0.36364, 0.00000,
    // 0.63636, 0.72727, 0.09091,
    // 0.27273, 0.90909, 0.45455,
    // 0.54545, 0.81818, 0.18182,
    // 0.18182, 0.81818, 0.54545,
    // 0.45455, 0.90909, 0.27273,
    // 0.09091, 0.72727, 0.63636,
    // 0.36364, 1.00000, 0.00000,
    // normalized
    F(0.73333), F(0.26667), F(0.00000),
    F(0.43750), F(0.50000), F(0.06250),
    F(0.16667), F(0.55556), F(0.27778),
    F(0.35294), F(0.52941), F(0.11765),
    F(0.11765), F(0.52941), F(0.35294),
    F(0.27778), F(0.55556), F(0.16667),
    F(0.06250), F(0.50000), F(0.43750),
    F(0.26667), F(0.73333), F(0.00000),
};

static constexpr decimal BILIN_WEIGHTS_Y_5[6*5] = {
    // not normalized
    // 1.00000, 0.63636, 0.27273, 0.00000, 0.00000, 0.00000,
    // 0.36364, 0.72727, 0.90909, 0.54545, 0.18182, 0.00000,
    // 0.09091, 0.45455, 0.81818, 0.81818, 0.45455, 0.09091,
    // 0.18182, 0.54545, 0.90909, 0.72727, 0.36364, 0.00000,
    // 0.27273, 0.63636, 1.00000, 0.00000, 0.00000, 0.00000,
    // normalized
    F(0.52381), F(0.33333), F(0.14286), F(0.00000), F(0.00000), F(0.00000),
    F(0.13333), F(0.26667), F(0.33333), F(0.20000), F(0.06667), F(0.00000),
    F(0.03333), F(0.16667), F(0.30000), F(0.30000), F(0.16667), F(0.03333),
    F(0.06667), F(0.20000), F(0.33333), F(0.26667), F(0.13333), F(0.00000),
    F(0.14286), F(0.33333), F(0.52381), F(0.00000), F(0.00000), F(0.00000),
};

/** Test if two floating point numbers are almost equal */
static inline bool is_close_enough(decimal a, decimal b)
{
    return ( fabs(a - b) <= EPSILON );
}

/** Calculate "tent" function, i.e., piece-wise linear interpolation
 *
 * The function is 1.0 at the center and goes linearly to 0.0 towards
 * (center - half_span) and (center + half_span).
 */
static decimal tent(decimal x, decimal center, decimal half_span)
{
    if ( is_close_enough(x, center) )
    {
        return F(1.0);
    }

    if ( x > center )
    {
        // mirror x across center and follow with left-side calculation
        x = center - (x - center);
    }

    decimal left = center - half_span;
    if ( is_close_enough(x, left) || (x < left) )
    {
        return F(0.0);
    }

    return (x - left) / (center - left);
}

/** See top header for description */
int populate_bilinear_weights(
    int w_inp,
    int h_inp,
    bilinear_weights* bw,
    int w_out,
    int h_out
){
    if ( (w_out > w_inp)
        || (h_out > h_inp)
        || (w_inp > astc::MAX_BLOCK_DIM)
        || (h_inp > astc::MAX_BLOCK_DIM)
        || (w_out > astc::MAX_GRID_DIM)
        || (h_out > astc::MAX_GRID_DIM) )
    {
        return 1;
    }

    for (int i = 0; i < astc::MAX_GRID_DIM; ++i)
    {
        bw->bilin_pixel_count_x[i] = 0;
        bw->bilin_pixel_count_y[i] = 0;
        for (int j = 0; j < astc::MAX_BLOCK_DIM; ++j)
        {
            bw->bilin_idx_x[i][j] = 0;
            bw->bilin_idx_y[i][j] = 0;
            bw->bilin_weights_x[i][j] = 0;
            bw->bilin_weights_y[i][j] = 0;
        }
    }

    decimal step_x = F(1.0) / ((decimal)(w_inp - 1));
    decimal step_y = F(1.0) / ((decimal)(h_inp - 1));
    decimal step_m = F(1.0) / ((decimal)(w_out - 1));
    decimal step_n = F(1.0) / ((decimal)(h_out - 1));

    // horizontal pass
    decimal m = F(0.0);
    for (int i = 0; i < w_out; ++i)
    {
        decimal x = F(0.0);
        uint8_t count = 0;
        decimal weight_sum = F(0.0);
        for (int j = 0; j < w_inp; ++j)
        {
            decimal weight = tent(x, m, step_m);
            if (weight > F(0.0))
            {
                bw->bilin_idx_x[i][count] = j;
                bw->bilin_weights_x[i][count] = weight;
                weight_sum += weight;
                ++count;
            }
            x += step_x;
        }
        bw->bilin_pixel_count_x[i] = count;

        // normalize weights to sum up to 1 (not doing this at runtime)
        if (weight_sum > F(0.0))
        {
            for (int j = 0; j < count; ++j)
            {
                bw->bilin_weights_x[i][j] /= weight_sum;
            }
        }

        m += step_m;
    }

    // vertical pass
    decimal n = F(0.0);
    for (int i = 0; i < h_out; ++i)
    {
        decimal y = F(0.0);
        uint8_t count = 0;
        decimal weight_sum = F(0.0);
        for (int j = 0; j < h_inp; ++j)
        {
            decimal weight = tent(y, n, step_n);
            if (weight > F(0.0))
            {
                bw->bilin_idx_y[i][count] = j;
                bw->bilin_weights_y[i][count] = weight;
                weight_sum += weight;
                ++count;
            }
            y += step_y;
        }
        bw->bilin_pixel_count_y[i] = count;

        // normalize weights to sum up to 1 (not doing this at runtime)
        if (weight_sum > F(0.0))
        {
            for (int j = 0; j < count; ++j)
            {
                bw->bilin_weights_y[i][j] /= weight_sum;
            }
        }

        n += step_n;
    }
    printf("\n");

    return 0;
}

/** Old code for reference */
/*
static void downsample_ref(
    const decimal *__restrict__ inp,
    int w_inp,
    int h_inp,
    const bilinear_weights *__restrict__ bw,
    decimal *__restrict__ out,
    int w_out,
    int h_out
){
    ZoneScopedN("bilin");

    // Buffer for holding intermediate results
    decimal tmp[astc::MAX_BLOCK_DIM*astc::MAX_GRID_DIM];

    // First, interpolate rows.
    { ZoneScopedN("rows");
    for (int y = 0; y < h_inp; ++y)
    {
        for (int m = 0; m < w_out; ++m)
        {
            // uint8_t pixel_count = bw->bilin_pixel_count_x[m];
            constexpr uint8_t pixel_count = 3;
            decimal weight_sum = FLT_MINVAL;  // prevent division by 0
            decimal out_pixel = F(0.0);
            for (uint8_t x = 0; x < pixel_count; ++x)
            {
                uint8_t idx = bw->bilin_idx_x[m][x];
                decimal inp_pixel = inp[y*w_inp+idx];
                decimal weight = bw->bilin_weights_x[m][x];
                out_pixel += inp_pixel * weight;
                weight_sum += weight;
            }
            // the weights do not sum up to 1 => we need to normalize
            out_pixel /= weight_sum;
            tmp[y*w_out+m] = out_pixel;
        }
    }
    }


    // Next, columns
    { ZoneScopedN("cols");
    for (int m = 0; m < w_out; ++m)
    {
        for (int n = 0; n < h_out; ++n)
        {
            // uint8_t pixel_count = bw->bilin_pixel_count_y[n];
            constexpr uint8_t pixel_count = 6;
            decimal weight_sum = FLT_MINVAL;  // prevent division by 0
            decimal out_pixel = F(0.0);
            for (uint8_t y = 0; y < pixel_count; ++y)
            {
                uint8_t idx = bw->bilin_idx_y[n][y];
                decimal inp_pixel = tmp[idx*w_out+m];
                decimal weight = bw->bilin_weights_y[n][y];
                out_pixel += inp_pixel * weight;
                weight_sum += weight;
            }
            // the weights do not sum up to 1 => we need to normalize
            out_pixel /= weight_sum;
            out[n*w_out+m] = out_pixel;
        }
    }
    }
}
*/

/** See top header for description */
void downsample(
    const decimal *__restrict__ inp,
    int w_inp,
    int h_inp,
    const bilinear_weights *__restrict__ bw,
    decimal *__restrict__ out,
    int w_out,
    int h_out
){
    ZoneScopedN("bilin");

    // Buffer for holding intermediate results
    decimal tmp[astc::MAX_BLOCK_DIM*astc::MAX_GRID_DIM];

    // First, interpolate rows.
    { ZoneScopedN("rows");
    for (int y = 0; y < h_inp; ++y)
    {
        // this loop vectorizes on x86_64 but not on armv8-a
        for (int m = 0; m < w_out; ++m)
        {
            uint8_t pixel_count = bw->bilin_pixel_count_x[m];
            decimal out_pixel = F(0.0);
            for (uint8_t x = 0; x < pixel_count; ++x)
            {
                const uint8_t idx = bw->bilin_idx_x[m][x];
                const decimal inp_pixel = inp[y*w_inp+idx];
                const decimal weight = bw->bilin_weights_x[m][x];
                out_pixel += inp_pixel * weight;
            }
            // the weights do not sum up to 1 => we need to normalize
            tmp[y*w_out+m] = out_pixel;
        }
    }
    }

    // Next, columns
    { ZoneScopedN("cols");
    for (int n = 0; n < h_out; ++n)
    {
        for (int m = 0; m < w_out; ++m)
        {
            uint8_t pixel_count = bw->bilin_pixel_count_y[n];
            decimal out_pixel = F(0.0);
            for (uint8_t y = 0; y < pixel_count; ++y)
            {
                const uint8_t idx = bw->bilin_idx_y[n][y];
                const decimal inp_pixel = tmp[idx*w_out+m];
                const decimal weight = bw->bilin_weights_y[n][y];
                out_pixel += inp_pixel * weight;
            }
            // the weights do not sum up to 1 => we need to normalize
            out[n*w_out+m] = out_pixel;
        }
    }
    }
}

/** See top header for description */
void downsample_12x12_to_8x5(
    const decimal *__restrict__ inp,
    decimal *__restrict__ out
){
    ZoneScopedN("bilin");

    constexpr int w_inp = 12;
    constexpr int h_inp = 12;
    constexpr int w_out = 8;
    constexpr int h_out = 5;

    // size of the bilin. kernel
    constexpr uint8_t pixel_count_x = 3;
    constexpr uint8_t pixel_count_y = 6;

    // Buffer for holding intermediate results
    decimal tmp[astc::MAX_BLOCK_DIM*astc::MAX_GRID_DIM];

    // First, interpolate rows.
    { ZoneScopedN("rows");
    for (int y = 0; y < h_inp; ++y)
    {
        for (int m = 0; m < w_out; ++m)
        {
            decimal out_pixel = F(0.0);
            for (uint8_t x = 0; x < pixel_count_x; ++x)
            {
                const uint8_t idx = BILIN_IDX_X_8[m][x];
                const decimal inp_pixel = inp[y*w_inp+idx];
                const decimal weight = BILIN_WEIGHTS_X_8[m*pixel_count_x+x];
                out_pixel += inp_pixel * weight;
            }
            tmp[y*w_out+m] = out_pixel;
        }
    }
    }

    // Next, columns
    { ZoneScopedN("cols");
    for (int n = 0; n < h_out; ++n)
    {
        for (int m = 0; m < w_out; ++m)
        {
            decimal out_pixel = F(0.0);
            for (uint8_t y = 0; y < pixel_count_y; ++y)
            {
                const uint8_t idx = BILIN_IDX_Y_5[n][y];
                const decimal inp_pixel = tmp[idx*w_out+m];
                const decimal weight = BILIN_WEIGHTS_Y_5[n*pixel_count_y+y];
                out_pixel += inp_pixel * weight;
            }
            out[n*w_out+m] = out_pixel;
        }
    }
    }
}

} // namespace simple::bilin
