#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "simple_bcn_common.h"
#include "simple_bc1.h"
#include "simple_ycocg_bc3.h"

namespace fs = std::filesystem;

/* Possible encoding formats */
typedef enum enc_format_t {
    BC1,
    YCOCG_BC3,
} enc_format_t;

/* Define encoding format here */
const enc_format_t ENC_FORMAT = YCOCG_BC3;

/* Print usage */
void err_exit(const std::string& err_msg, bool print_usage=false)
{
    fprintf(stderr, "ERROR: %s\n", err_msg.data());
    if (print_usage)
    {
        fprintf(stderr,
            "Usage: transcoder input_image [<input_image> ...] output_folder\n"
        );
    }
    exit(1);
}

/* Encode the whole image according to ENC_FORMAT
 *
 * Returns zero on success, non-zero on failure.
 */
int encode_image(
    const uint8_t *inp_pixels,
    int img_w,
    int img_h,
    std::vector<uint32_t>& enc_data
){
    // Number of blocks in horizontal and vertical dimension
    const int nblocks_x = img_w / 4;
    const int nblocks_y = img_h / 4;

    void (*encode_block)( const uint8_t[NCH_RGB*16], uint32_t[2] );

    // Based on encoding format:
    //   1. Set encoded block size as a number of 32-bit integers
    //   2. Assign decoding function
    int block_nints = 0;
    switch (ENC_FORMAT) {
    case BC1:
        block_nints = 2;  // 64 bits per 4x4 block
        encode_block = encode_block_bc1;
        break;
    case YCOCG_BC3:
        block_nints = 4; // 128 bits per 4x4 block
        encode_block = encode_block_ycocg_bc3;
        break;
    default:
        fprintf(stderr, "ERROR: Unsupported encoding format\n");
        return 1;
    };

	// Resize the encoded data buffer
	enc_data.resize(nblocks_x * nblocks_y * block_nints);

    // Iterate through blocks and encode them one-by-one (this can be easily
    // parallelized, all blocks are independent to each other)
    for (int block_y = 0; block_y < nblocks_y; ++block_y)
    {
        for (int block_x = 0; block_x < nblocks_x; ++block_x)
        {
            // Read 4x4 block of pixels into an array
            uint8_t block_pixels[NCH_RGB*16];
            const int x = block_x * 4;
            const int y = block_y * 4;
            for (int i = 0; i < 4; ++i)
            {
                memcpy(
                    block_pixels + (NCH_RGB * i * 4),
                    inp_pixels + (NCH_RGB * (img_w * (y+i) + x)),
                    NCH_RGB * 4
                );
            }

            // Encode the block
            const int offset = (nblocks_x * block_y + block_x) * block_nints;
            encode_block(block_pixels, enc_data.data() + offset);
        }
    }

    return 0;
}

/* Decode the whole image according to ENC_FORMAT
 *
 * Returns zero on success, non-zero on failure.
 */
int decode_image(
    const std::vector<uint32_t>& enc_data,
    int img_w,
    int img_h,
    std::vector<uint8_t>& out_pixels
){
    // Number of blocks in horizontal and vertical dimension
    const int nblocks_x = img_w / 4;
    const int nblocks_y = img_h / 4;

    void (*decode_block)( const uint32_t[NCH_RGB*16], uint8_t[2] );

    // Based on encoding format:
    //   1. Set encoded block size as a number of 32-bit integers
    //   2. Assign decoding function
    int block_nints = 0;
    switch (ENC_FORMAT) {
    case BC1:
        block_nints = 2;  // 64 bits per 4x4 block
        decode_block = decode_block_bc1;
        break;
    case YCOCG_BC3:
        block_nints = 4; // 128 bits per 4x4 block
        decode_block = decode_block_ycocg_bc3;
        break;
    default:
        fprintf(stderr, "ERROR: Unsupported encoding format\n");
        return 1;
    };

    // Iterate through blocks and encode them one-by-one (this could be easily
    // parallelized, all blocks are independent to each other)
    for (int block_y = 0; block_y < nblocks_y; ++block_y)
    {
        for (int block_x = 0; block_x < nblocks_x; ++block_x)
        {
            const int x = block_x * 4;
            const int y = block_y * 4;
            const int offset = (nblocks_x * block_y + block_x) * block_nints;

            // Decode
            uint8_t decoded_pixels[NCH_RGB*16];
            decode_block(enc_data.data() + offset, decoded_pixels);

            // Write decoded pixels into the output image
            for (int i = 0; i < 4; ++i)
            {
                memcpy(
                    out_pixels.data() + (NCH_RGB * (img_w * (y+i) + x)),
                    decoded_pixels + (NCH_RGB * i * 4),
                    NCH_RGB * 4
                );
            }
        }
    }

    return 0;
}

/* Dump raw encoded data into a file (for debugging) */
void dump_enc_data(
    const std::vector<uint32_t>& enc_data,
    const std::string& filename
){
    printf("-- Saving encoded data to '%s'\n", filename.data());
    std::ofstream wf(filename);
    for (const uint32_t& x : enc_data)
    {
        wf << (uint8_t)( x & 0xff );
        wf << (uint8_t)( (x >> 8) & 0xff );
        wf << (uint8_t)( (x >> 16) & 0xff );
        wf << (uint8_t)( (x >> 24) & 0xff );
    }
    wf.close();
}

/* Read the current time and return seconds (taken from astcenc) */
double get_time()
{
    timeval tv;
    gettimeofday(&tv, 0);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1.0e-6;
}

int main(int argc, char **argv)
{
    // Check command line arguments
    if (argc < 3)
    {
        err_exit("Provide at least one image and an output folder", true);
    }

    // Last argument is the output directory
    std::string out_dir = argv[argc-1];
    printf("Output directory: '%s'\n", out_dir.data());
    if ( !(fs::is_directory(out_dir) && fs::exists(out_dir)) )
    {
        err_exit("Can't open output directory", true);
    }

    // Error code
    int err = 0;

    double total_enc_duration = 0.0;
    uint64_t total_num_enc_pixels = 0;
    int num_enc_images = 0;

    // Loop through images one by one
    for (int i = 1; i < argc-1; ++i)
    {
        std::string inp_name = argv[i];
        printf("Image %d/%d: %s\n", i, argc-2, inp_name.data());

        // Read input image
        int inp_w, inp_h, nch;
        uint8_t *inp_pixels = stbi_load(
            inp_name.data(),
            &inp_w,
            &inp_h,
            &nch,
            NCH_RGB
        );
        if (inp_pixels == NULL)
        {
            printf("-- Error opening file, skipping\n");
            err = 1;
            continue;
        }

        // Image width and height must be multiples of 4 (would have to pad...)
        if ( (inp_w % 4 != 0) || (inp_h % 4 != 0) )
        {
            printf("-- Image size must be divisible by 4, skipping\n");
            err = 1;
            continue;
        }

        double enc_start_time = get_time();

        // Encode it into enc_data
        std::vector<uint32_t> enc_data;
        if (encode_image(inp_pixels, inp_w, inp_h, enc_data))
        {
            continue;
        }

        total_enc_duration += ( get_time() - enc_start_time );
        total_num_enc_pixels += ( inp_w * inp_h );
        num_enc_images += 1;

        // Dump encoded data to file as raw bits
        // std::string dump_name = out_dir
        //     / (fs::path(inp_name).filename().replace_extension(".bin"));
        // dump_enc_data(enc_data, dump_name);

        // Decode enc_data into out_pixels
        std::vector<unsigned char> out_pixels(inp_w*inp_h*NCH_RGB);
        if (decode_image(enc_data, inp_w, inp_h, out_pixels))
        {
            continue;
        }

        // Save decoded image into PNG file inside out_dir
        std::string out_name = out_dir / fs::path(inp_name).filename();
        printf("-- Saving decoded image to '%s'\n", out_name.data());
        int ret = stbi_write_png(
            out_name.data(),
            inp_w,
            inp_h,
            NCH_RGB,
            out_pixels.data(),
            inp_w*NCH_RGB
        );
        if (ret == 0)
        {
            err_exit("Can't save output image");
        }

        stbi_image_free(inp_pixels);
    }

    printf("Average encoding time (sec)   : %.6f\n", total_enc_duration / num_enc_images);
    printf("Average encoding rate (Mpx/s) : %.3f\n", total_num_enc_pixels / total_enc_duration / 1e6);

    return err;
}
