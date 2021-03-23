/** Small utility to downsample an image using bilinear filtering
 *
 * Designed to work on small images.
 */

#include <cmath>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

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
    if ((w_out > w_inp) || (h_out > h_inp))
    {
        return 1;
    }

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
        for (int j = 0; j < w_inp; ++j)
        {
            float weight = tent(x, m, step_m);
            if (weight > 0.0f)
            {
                printf("- x[%2d]: %13.10f  tent: %13.10f\n", j, x, weight);
            }
            x += step_x;
        }
        m += step_m;
    }

    float n = 0.0f;
    printf("\nstep_y: %6.3f, step_n: %6.3f\n", step_y, step_n);
    for (int i = 0; i < h_out; ++i)
    {
        printf("  n[%2d]: %6.3f\n", i, n);
        float y = 0.0f;
        for (int j = 0; j < h_inp; ++j)
        {
            float weight = tent(y, n, step_n);
            if (weight > 0.0f)
            {
                printf("- y[%2d]: %13.10f  tent: %13.10f\n", j, y, weight);
            }
            y += step_y;
        }
        n += step_n;
    }
    printf("\n");

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
    int res = bilinear_downsample(
        nullptr,
        nch, inp_w, inp_h,
        nullptr,
        M, N
    );
    if (res != 0)
    {
        err_exit("Error downsampling input image");
    }

    // Save output image
    printf("Saving result to: %s\n", out_file.c_str());

    stbi_image_free(inp_pixels);
}
