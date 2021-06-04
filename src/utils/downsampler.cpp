/** Small utility to downsample an image using bilinear filtering
 *
 * Designed to work on small images.
 */

#include <array>
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

#include "simple_texcomp.hpp"

using namespace simple;

/* Exit with error, optionally printing usage */
static void err_exit(const char* err_msg, bool print_usage=false)
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

static void print_bilinear_weights(const bilin::bilinear_weights* bw)
{
    printf("bilin_pixel_count_x = {\n");
    for (int i = 0; i < astc::MAX_GRID_DIM; ++i)
    {
        printf("%3d,", bw->bilin_pixel_count_x[i]);
    }
    printf("\n};\n");

    printf("\nbilin_idx_x = {\n");
    for (int i = 0; i < astc::MAX_GRID_DIM; ++i)
    {
        printf("    {");
        for (int j = 0; j < astc::MAX_BLOCK_DIM; ++j)
        {
            printf("%3d,", bw->bilin_idx_x[i][j]);
        }
        printf(" },\n");
    }
    printf("};\n");

    printf("\nbilin_weights_x = {\n");
    for (int i = 0; i < astc::MAX_GRID_DIM; ++i)
    {
        printf("    {\n");
        for (int j = 0; j < astc::MAX_BLOCK_DIM; ++j)
        {
            if (j % 6 == 0)
            {
                printf("        ");
            }
            printf("%7.5f, ", (double)bw->bilin_weights_x[i][j]);
            if (j % 6 == 5)
            {
                printf("\n");
            }
        }
        printf("    },\n");
    }
    printf("};\n");

    printf("\nbilin_pixel_count_y = {\n");
    for (int i = 0; i < astc::MAX_GRID_DIM; ++i)
    {
        printf("%3d,", bw->bilin_pixel_count_y[i]);
    }
    printf("\n};\n");

    printf("\nbilin_idx_y = {\n");
    for (int i = 0; i < astc::MAX_GRID_DIM; ++i)
    {
        printf("    {");
        for (int j = 0; j < astc::MAX_BLOCK_DIM; ++j)
        {
            printf("%3d,", bw->bilin_idx_y[i][j]);
        }
        printf(" },\n");
    }
    printf("};\n");

    printf("\nbilin_weights_y = {\n");
    for (int i = 0; i < astc::MAX_GRID_DIM; ++i)
    {
        printf("    {\n");
        for (int j = 0; j < astc::MAX_BLOCK_DIM; ++j)
        {
            if (j % 6 == 0)
            {
                printf("        ");
            }
            printf("%7.5f, ", (double)bw->bilin_weights_y[i][j]);
            if (j % 6 == 5)
            {
                printf("\n");
            }
        }
        printf("    },\n");
    }
    printf("};\n");
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

    // Pre-compute bilinear weights
    bilin::bilinear_weights bw;
    int ret = populate_bilinear_weights(inp_w, inp_h, &bw, M, N);
    if (ret != 0)
    {
        err_exit("Error pre-computing weights. Check dimensions.");
    }

    // Print out the data structures
    print_bilinear_weights(&bw);

    // Downsample
    std::vector<decimal> inp_pixels_flt(inp_w*inp_h*nch);
    for (int i = 0; i < inp_h*inp_w*nch; ++i)
    {
        inp_pixels_flt.at(i) = (decimal)(inp_pixels[i]) / F(255.0);
    }

    std::vector<decimal> out_pixels_flt(M*N*nch);

    bilin::downsample(
        inp_pixels_flt.data(),
        inp_w, inp_h,
        // &bw,
        out_pixels_flt.data(),
        M, N
    );

    // Save output image
    printf("Saving result to: %s\n", out_file.c_str());
    std::vector<uint8_t> out_pixels(M*N*nch);
    for (size_t i = 0; i < out_pixels.size(); ++i)
    {
        decimal out_pixel_flt = out_pixels_flt.at(i);
        if (out_pixel_flt > F(1.0))
        {
            printf("ERROR: Result pixel > 1.0: %ld, %.6f\n",
                i, (double)out_pixel_flt);
        }
        out_pixels.at(i) = (uint8_t)(out_pixel_flt * F(255.0));
    }

    ret = stbi_write_png(
        out_file.data(),
        M,
        N,
        nch,
        out_pixels.data(),
        M*nch
    );
    if (ret == 0)
    {
        err_exit("Can't save output image");
    }

    stbi_image_free(inp_pixels);
}
