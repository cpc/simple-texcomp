/** Small utility to downsample an image using bilinear filtering
 *
 * Designed to work on small images.
 */

#include <cmath>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

// maximum dimension of weight grid stored in ASTC block
#define MAX_WEIGHT_GRID_SIZE  10
// maximum input block size
#define MAX_BLOCK_SIZE        12

/* Exit with error, optionally printing usage */
void err_exit(const char* err_msg, bool print_usage=false)
{
    fprintf(stderr, "ERROR: %s\n", err_msg);
    if (print_usage)
    {
        fprintf(stderr,
            "Usage: downsampler <inp_img> <out_img> M N\n"
            "\nM, N ... output size\n"
        );
    }
    exit(1);
}

/** Test if two floating point numbers are almost equal */
static inline bool is_close_enough(float a, float b)
{
    static float epsilon = 1e-6;
    return ( std::fabs(a - b) <= epsilon );
}

/** Calculate "tent" function, i.e. piece-wise linear interpolation
 *
 * The function is 1.0 at the center and goes linearly to 0.0 towards
 * (center - half_span) and (center + half_span).
 */
static float tent(float x, float center, float half_span)
{
    if ( is_close_enough(x, center) )
    {
        return 1.0f;
    }

    if ( x > center )
    {
        // mirror x across center and follow with left-side calculation
        x = center - (x - center);
    }

    float left = center - half_span;
    if ( is_close_enough(x, left) || (x < left) )
    {
        return 0.0f;
    }

    return (x - left) / (center - left);
}

/** Downsample input block XxY into the size of MxN
 *
 * Returns 1 in case of error, 0 on success
 */
int bilinear_downsample(
    const float* inp,  // input values
    int nch,           // number of channels (e.g. 3 for RGB or 1 for grayscale)
    int w_inp,         // width of the input block (X)
    int h_inp,         // height of the input block (Y)
    float* out,        // output values
    int w_out,         // width of the output block (M; M <= X)
    int h_out          // height of the output block (N; N <= Y)
){
    if ( (w_out > w_inp)
        || (h_out > h_inp)
        || (w_inp > MAX_BLOCK_SIZE)
        || (h_inp > MAX_BLOCK_SIZE)
        || (w_out > MAX_WEIGHT_GRID_SIZE)
        || (h_out > MAX_WEIGHT_GRID_SIZE) )
    {
        return 1;
    }

    if (nch != 1)
    {
        return 1;
    }

    // Data structures for keeping pre-computed bilinear weights. Stored
    // separately for X and Y dimension (to use M+N instead of M*N memory)
    //
    // How many pixels go to a calculation of each weight
    uint8_t bilin_pixel_count_x[MAX_WEIGHT_GRID_SIZE];
    uint8_t bilin_pixel_count_y[MAX_WEIGHT_GRID_SIZE];
    // Each bilin_idx_x[i] stores an array of indices that point to locations
    // of pixels used to calculate that weight
    uint8_t bilin_idx_x[MAX_WEIGHT_GRID_SIZE][MAX_BLOCK_SIZE];
    uint8_t bilin_idx_y[MAX_WEIGHT_GRID_SIZE][MAX_BLOCK_SIZE];
    // Each bilin_weights_x[i] stores the actual weight values the previous
    // array points at
    float bilin_weights_x[MAX_WEIGHT_GRID_SIZE][MAX_BLOCK_SIZE];
    float bilin_weights_y[MAX_WEIGHT_GRID_SIZE][MAX_BLOCK_SIZE];

    for (int i = 0; i < MAX_WEIGHT_GRID_SIZE; ++i)
    {
        bilin_pixel_count_x[i] = 0;
        bilin_pixel_count_y[i] = 0;
        for (int j = 0; j < MAX_BLOCK_SIZE; ++j)
        {
            bilin_idx_x[i][j] = 0;
            bilin_idx_y[i][j] = 0;
            bilin_weights_x[i][j] = 0;
            bilin_weights_y[i][j] = 0;
        }
    }

    // Popuate the data structures above -- can be done at init phase or even
    // set as a const array at compile time if the block sizes are fixed.
    // Also, the weights can be normalized  at this step instead of runtime.
    float step_x = 1.0f / ((float)(w_inp - 1));
    float step_y = 1.0f / ((float)(h_inp - 1));
    float step_m = 1.0f / ((float)(w_out - 1));
    float step_n = 1.0f / ((float)(h_out - 1));

    printf("\ninp size: %dx%d, out size: %dx%d\n", w_inp, h_inp, w_out, h_out);

    float m = 0.0f;
    printf("\nstep_x: %6.3f, step_m: %6.3f\n", step_x, step_m);
    for (int i = 0; i < w_out; ++i)
    {
        printf("  m[%2d]: %6.3f\n", i, m);
        float x = 0.0f;
        uint8_t count = 0;
        for (int j = 0; j < w_inp; ++j)
        {
            float weight = tent(x, m, step_m);
            if (weight > 0.0f)
            {
                printf("- x[%2d]: %13.10f  tent: %13.10f\n", j, x, weight);
                bilin_idx_x[i][count] = j;
                bilin_weights_x[i][count] = weight;
                ++count;
            }
            x += step_x;
        }
        bilin_pixel_count_x[i] = count;
        m += step_m;
    }

    float n = 0.0f;
    printf("\nstep_y: %6.3f, step_n: %6.3f\n", step_y, step_n);
    for (int i = 0; i < h_out; ++i)
    {
        printf("  n[%2d]: %6.3f\n", i, n);
        float y = 0.0f;
        uint8_t count = 0;
        for (int j = 0; j < h_inp; ++j)
        {
            float weight = tent(y, n, step_n);
            if (weight > 0.0f)
            {
                printf("- y[%2d]: %13.10f  tent: %13.10f\n", j, y, weight);
                bilin_idx_y[i][count] = j;
                bilin_weights_y[i][count] = weight;
                ++count;
            }
            y += step_y;
        }
        bilin_pixel_count_y[i] = count;
        n += step_n;
    }
    printf("\n");

    // Print out the data structures to double check or use for static init
    bool printout = false;
    if (printout)
    {
        printf("bilin_pixel_count_x = {\n");
        for (int i = 0; i < MAX_WEIGHT_GRID_SIZE; ++i)
        {
            printf("%3d,", bilin_pixel_count_x[i]);
        }
        printf("\n};\n");

        printf("\nbilin_pixel_count_y = {\n");
        for (int i = 0; i < MAX_WEIGHT_GRID_SIZE; ++i)
        {
            printf("%3d,", bilin_pixel_count_y[i]);
        }
        printf("\n};\n");

        printf("\nbilin_idx_x = {\n");
        for (int i = 0; i < MAX_WEIGHT_GRID_SIZE; ++i)
        {
            printf("    {");
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j)
            {
                printf("%3d,", bilin_idx_x[i][j]);
            }
            printf(" },\n");
        }
        printf("};\n");

        printf("\nbilin_idx_y = {\n");
        for (int i = 0; i < MAX_WEIGHT_GRID_SIZE; ++i)
        {
            printf("    {");
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j)
            {
                printf("%3d,", bilin_idx_y[i][j]);
            }
            printf(" },\n");
        }
        printf("};\n");

        printf("\nbilin_weights_x = {\n");
        for (int i = 0; i < MAX_WEIGHT_GRID_SIZE; ++i)
        {
            printf("    {\n");
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j)
            {
                if (j % 6 == 0)
                {
                    printf("        ");
                }
                printf("%7.5f, ", bilin_weights_x[i][j]);
                if (j % 6 == 5)
                {
                    printf("\n");
                }
            }
            printf("    },\n");
        }
        printf("};\n");

        printf("\nbilin_weights_y = {\n");
        for (int i = 0; i < MAX_WEIGHT_GRID_SIZE; ++i)
        {
            printf("    {\n");
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j)
            {
                if (j % 6 == 0)
                {
                    printf("        ");
                }
                printf("%7.5f, ", bilin_weights_y[i][j]);
                if (j % 6 == 5)
                {
                    printf("\n");
                }
            }
            printf("    },\n");
        }
        printf("};\n");
    }

    // Now go ahead and finally filter the values
    // First, rows.
    float tmp[h_inp*w_out];  // rows are filtered but columns are not

    printf("Input image:\n");
    for (int y = 0; y < h_inp; ++y)
    {
        for (int x = 0; x < w_inp; ++x)
        {
            printf("%4d", (int)(inp[y*w_inp+x] * 255.0f));
        }
        printf("\n");
    }

    for (int y = 0; y < h_inp; ++y)
    {
        // printf("row %d\n", y);
        for (int m = 0; m < w_out; ++m)
        {
            uint8_t pixel_count = bilin_pixel_count_x[m];
            // printf("weight %d, npixels: %d\n", m, pixel_count);
            float weight_sum = 1e-11f;  // prevent division by 0
            float out_pixel = 0.0f;
            for (uint8_t x = 0; x < pixel_count; ++x)
            {
                uint8_t idx = bilin_idx_x[m][x];
                float inp_pixel = inp[y*w_inp+idx];
                float weight = bilin_weights_x[m][x];
                out_pixel += inp_pixel * weight;
                weight_sum += weight;
                // printf("%3d %3d %3d %6.3f\n", x, idx, y*w_inp+idx, inp_pixel);
            }
            // the weights do not sum up to 1 => we need to normalize
            out_pixel /= weight_sum;
            tmp[y*w_out+m] = out_pixel;
            // printf("%6.3f\n", out_pixel);
        }
        // printf("\n");
    }

    printf("Tmp:\n");
    for (int y = 0; y < h_inp; ++y)
    {
        for (int x = 0; x < w_out; ++x)
        {
            printf("%4d", (int)(tmp[y*w_out+x] * 255.0f));
        }
        printf("\n");
    }

    // Next, columns
    for (int m = 0; m < w_out; ++m)
    {
        // printf("column %d\n", m);
        for (int n = 0; n < h_out; ++n)
        {
            uint8_t pixel_count = bilin_pixel_count_y[n];
            // printf("weight %d, npixels: %d\n", n, pixel_count);
            float weight_sum = 1e-11f;  // prevent division by 0
            float out_pixel = 0.0f;
            for (uint8_t y = 0; y < pixel_count; ++y)
            {
                uint8_t idx = bilin_idx_y[n][y];
                float inp_pixel = tmp[idx*w_out+m];
                float weight = bilin_weights_y[n][y];
                out_pixel += inp_pixel * weight;
                weight_sum += weight;
                // printf("%3d %3d %3d %6.3f %3d\n", n, idx, idx*w_inp+m,
                //     inp_pixel, (int)(inp_pixel * 255.0f));
            }
            // the weights do not sum up to 1 => we need to normalize
            out_pixel /= weight_sum;
            out[n*w_out+m] = out_pixel;
            // printf("%6.3f (%3d)\n", out_pixel, (int)(out_pixel * 255.0f));
        }
        // printf("\n");
    }

    printf("Out:\n");
    for (int y = 0; y < h_out; ++y)
    {
        for (int x = 0; x < w_out; ++x)
        {
            printf("%4d", (int)(out[y*w_out+x] * 255.0f));
        }
        printf("\n");
    }

    return 0;
}

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        err_exit("Provide 4 arguments", true);
    }

    const std::string inp_file = argv[1];
    const std::string out_file = argv[2];

    int M = atoi(argv[3]);
    int N = atoi(argv[4]);

    if ( (M < 2) || (N < 2) )
    {
        err_exit("Output size must be >= 2");
    }

    // Read input image
    int inp_w, inp_h, nch;
    uint8_t *inp_pixels = stbi_load(
        inp_file.data(),
        &inp_w,
        &inp_h,
        &nch,
        0
    );
    if (inp_pixels == NULL)
    {
        err_exit("Can't open input file");
    }

    printf("Input image: %s\n", inp_file.c_str());
    printf("%dx%d, %d channels\n", inp_w, inp_h, nch);

    // Downsample
    std::vector<float> inp_pixels_flt(inp_w*inp_h*nch);
    for (int i = 0; i < inp_h*inp_w*nch; ++i)
    {
        inp_pixels_flt.at(i) = (float)(inp_pixels[i]) / 255.0f;
    }

    std::vector<float> out_pixels_flt(M*N*nch);

    int res = bilinear_downsample(
        inp_pixels_flt.data(),
        nch, inp_w, inp_h,
        out_pixels_flt.data(),
        M, N
    );
    if (res != 0)
    {
        err_exit("Error downsampling input image");
    }

    // Save output image
    printf("Saving result to: %s\n", out_file.c_str());
    std::vector<uint8_t> out_pixels(M*N*nch);
    for (size_t i = 0; i < out_pixels.size(); ++i)
    {
        float out_pixel_flt = out_pixels_flt.at(i);
        if (out_pixel_flt > 1.0f)
        {
            printf("ERROR: Result pixel > 1.0: %ld, %.6f\n", i, out_pixel_flt);
        }
        out_pixels.at(i) = (uint8_t)(out_pixel_flt * 255.0f);
    }
    res = stbi_write_png(
        out_file.data(),
        M,
        N,
        nch,
        out_pixels.data(),
        M*nch
    );
    if (res == 0)
    {
        err_exit("Can't save output image");
    }

    stbi_image_free(inp_pixels);
}
