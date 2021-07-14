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

#include <filesystem.hpp>
namespace fs = ghc::filesystem;

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

    printf("\nbilin_weights_x_u8 = {\n");
    for (int i = 0; i < astc::MAX_GRID_DIM; ++i)
    {
        printf("    {\n");
        for (int j = 0; j < astc::MAX_BLOCK_DIM; ++j)
        {
            if (j % 6 == 0)
            {
                printf("        ");
            }
            printf("%3d, ", (int)(bw->bilin_weights_x[i][j] * F(256.0)));
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

    printf("\nbilin_weights_y_u8 = {\n");
    for (int i = 0; i < astc::MAX_GRID_DIM; ++i)
    {
        printf("    {\n");
        for (int j = 0; j < astc::MAX_BLOCK_DIM; ++j)
        {
            if (j % 6 == 0)
            {
                printf("        ");
            }
            printf("%3d, ", (int)(bw->bilin_weights_y[i][j] * F(256.0)));
            if (j % 6 == 5)
            {
                printf("\n");
            }
        }
        printf("    },\n");
    }
    printf("};\n");
}

static int save_image(
    std::vector<decimal> img_pixels,
    int out_w,
    int out_h,
    int nch,
    fs::path out_file
){
    std::vector<uint8_t> out_pixels(img_pixels.size());
    for (size_t i = 0; i < out_pixels.size(); ++i)
    {
        decimal out_pixel_flt = img_pixels.at(i);
        if (out_pixel_flt > F(1.0))
        {
            printf("ERROR: Result pixel > 1.0: %ld, %.6f\n",
                i, (double)out_pixel_flt);
        }
        out_pixels.at(i) = (uint8_t)(out_pixel_flt * F(255.0));
    }

    int ret = stbi_write_png(
        out_file.c_str(),
        out_w,
        out_h,
        nch,
        out_pixels.data(),
        out_w*nch
    );

    return ret;
}

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        err_exit("Provide 4 arguments", true);
    }

    const fs::path inp_file = argv[1];
    const fs::path out_file = argv[2];
    const fs::path out_file_12x12_to_8x5 = fs::path(argv[2]).replace_filename(
        out_file.stem().string() + "_12x12_to_8x5" + out_file.extension().string()
    );

    int M = atoi(argv[3]);
    int N = atoi(argv[4]);

    if ( (M < 2) || (N < 2) )
    {
        err_exit("Output size must be >= 2");
    }

    // Read input image
    int inp_w, inp_h, nch;
    uint8_t *inp_pixels = stbi_load(
        inp_file.c_str(),
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

    // Downsample (assume grayscale input)
    if ((nch != 1) && (nch != 3))
    {
        err_exit("Unsupported number of channels");
    }
    if (nch == 3)
    {
        printf("WARNING! Image seems not to be grayscale, sampling ITU BT.709 luminance\n");
    }

    std::vector<decimal> inp_pixels_flt(inp_w*inp_h);
    std::vector<uint8_t> inp_pixels_u8(inp_w*inp_h);

    for (int i = 0; i < inp_h*inp_w; ++i)
    {
        decimal y = F(0.0);

        if (nch == 3)
        {
            uint8_t r = inp_pixels[nch*i+0];
            uint8_t g = inp_pixels[nch*i+1];
            uint8_t b = inp_pixels[nch*i+2];

            y = (decimal)(r) / F(255.0) * F(0.2126) +
                (decimal)(g) / F(255.0) * F(0.7152) +
                (decimal)(b) / F(255.0) * F(0.0722);
        }
        else
        {
            y = (decimal)(inp_pixels[nch*i]) / F(255.0);
        }

        inp_pixels_flt.at(i) = y;
        inp_pixels_u8.at(i) = (uint8_t)(y * F(255.0));
    }

    std::vector<decimal> out_pixels_flt(M*N);

    bilin::downsample(
        inp_pixels_flt.data(),
        inp_w, inp_h,
        &bw,
        out_pixels_flt.data(),
        M, N
    );

    // Save output image
    printf("\nSaving result to: %s\n", out_file.c_str());
    ret = save_image(
        out_pixels_flt,
        M, N, 1,
        out_file
    );
    if (ret == 0)
    {
        err_exit("Can't save output image");
    }

    // Downsample & save fixed-sized images
    if ((inp_w == 12) && (inp_h == 12))
    {
        std::vector<decimal> out_pixels_flt_12x12_to_8x5(M*N);

        bilin::downsample_12x12_to_8x5(
            inp_pixels_flt.data(),
            out_pixels_flt_12x12_to_8x5.data()
        );

        printf("Saving result to: %s\n", out_file_12x12_to_8x5.c_str());
        ret = save_image(
            out_pixels_flt_12x12_to_8x5,
            M, N, 1,
            out_file_12x12_to_8x5
        );
        if (ret == 0)
        {
            err_exit("Can't save output image");
        }

        // Compute the integer variant
        printf("\n8-bit unsigned version comparison:\n");
        std::vector<uint8_t> out_pixels_u8_12x12_to_8x5(M*N);

        bilin::downsample_12x12_to_8x5_u8(
            inp_pixels_u8.data(),
            out_pixels_u8_12x12_to_8x5.data()
        );

        decimal sqerr = F(0.0);
        for (int i = 0; i < M*N; ++i)
        {
            decimal ref = out_pixels_flt_12x12_to_8x5.at(i);
            decimal tgt = (decimal)(out_pixels_u8_12x12_to_8x5.at(i)) / F(255.0);

            decimal diff = tgt - ref;
            printf("%3d:  ref: %10.8f  u8: %10.8f  diff: %+11.8f\n",
                i, (double)ref, (double)tgt, (double)diff);

            sqerr += diff * diff;
        }
        printf("sqerr: %.8f\n", (double)sqerr);
    }
    else
    {
        printf("WARNING: Input size is not 12x12, not downsampling with "
            "12x12_to_8x5\n");
    }

    stbi_image_free(inp_pixels);
}
