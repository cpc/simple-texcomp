#include <cassert>
#include <cstdio>
#include <cmath>

#include "simple_texcomp.hpp"

namespace simple::bilin {

/** Test if two floating point numbers are almost equal */
static inline bool is_close_enough(decimal a, decimal b)
{
    return ( std::fabs(a - b) <= EPSILON );
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
        for (int j = 0; j < w_inp; ++j)
        {
            decimal weight = tent(x, m, step_m);
            if (weight > F(0.0))
            {
                bw->bilin_idx_x[i][count] = j;
                bw->bilin_weights_x[i][count] = weight;
                ++count;
            }
            x += step_x;
        }
        bw->bilin_pixel_count_x[i] = count;
        m += step_m;
    }

    // vertical pass
    decimal n = F(0.0);
    for (int i = 0; i < h_out; ++i)
    {
        decimal y = F(0.0);
        uint8_t count = 0;
        for (int j = 0; j < h_inp; ++j)
        {
            decimal weight = tent(y, n, step_n);
            if (weight > F(0.0))
            {
                bw->bilin_idx_y[i][count] = j;
                bw->bilin_weights_y[i][count] = weight;
                ++count;
            }
            y += step_y;
        }
        bw->bilin_pixel_count_y[i] = count;
        n += step_n;
    }
    printf("\n");

    return 0;
}

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
    // Buffer for holding intermediate results
    static decimal tmp[astc::MAX_BLOCK_DIM*astc::MAX_GRID_DIM];

    // First, interpolate rows.
    for (int y = 0; y < h_inp; ++y)
    {
        for (int m = 0; m < w_out; ++m)
        {
            uint8_t pixel_count = bw->bilin_pixel_count_x[m];
            decimal weight_sum = std::numeric_limits<decimal>::min();  // prevent division by 0
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
            assert(!std::isnan(out_pixel));
            tmp[y*w_out+m] = out_pixel;
        }
    }

    // Next, columns
    for (int m = 0; m < w_out; ++m)
    {
        for (int n = 0; n < h_out; ++n)
        {
            uint8_t pixel_count = bw->bilin_pixel_count_y[n];
            decimal weight_sum = std::numeric_limits<decimal>::min();  // prevent division by 0
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
            assert(!std::isnan(out_pixel));
            out[n*w_out+m] = out_pixel;
        }
    }
}

} // namespace
