#include <cassert>
#include <cstdio>
#include <cstdint>
#include <cstdlib>

#include <sys/time.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "simple_texcomp.hpp"

namespace fs = std::filesystem;
using namespace simple;

/* Possible encoding formats */
typedef enum enc_format_t {
    BC1,
    YCOCG_BC3,
    ASTC,
} enc_format_t;

/* Define encoding format here */
const enc_format_t ENC_FORMAT = ASTC;

/* Exit with error, optionally printing usage */
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

/* Pad image so that width and height are divisible by 4
 *
 * Border pixels are repeated
 */
void pad_image(
    const uint8_t *inp_pixels,
    int img_w,
    int img_h,
    std::vector<uint8_t>& padded_img,
    int& pad_w,
    int& pad_h
){
    // Input block size
    const int block_w = (ENC_FORMAT == ASTC) ? astc::BLOCK_X : 4;
    const int block_h = (ENC_FORMAT == ASTC) ? astc::BLOCK_Y : 4;

    // Number of blocks in horizontal and vertical dimension
    const int nblocks_x = ( img_w + (block_w - 1) ) / block_w;
    const int nblocks_y = ( img_h + (block_h - 1) ) / block_h;

    // Image size after padding
    pad_w = nblocks_x * block_w;
    pad_h = nblocks_y * block_h;

	padded_img.resize(pad_w * pad_h * NCH_RGB);
    for (int i = 0; i < pad_w*pad_h*NCH_RGB; ++i)
    {
        padded_img[i] = 0;
    }

    for (int y = 0; y < pad_h; ++y)
    {
        // Repeat the last row that goes over the image height
        int y_in = std::min(y, img_h-1);

        // Copy the image row
        memcpy(
            padded_img.data() + (y * pad_w * NCH_RGB), // dest
            inp_pixels + (y_in * img_w * NCH_RGB),     // src
            img_w * NCH_RGB                            // num_bytes
        );

        // Repeat the last pixel of the row
        for (int x = img_w; x < pad_w; ++x)
        {
            int i1 = NCH_RGB * (y * pad_w + x);
            int i2 = NCH_RGB * (y * pad_w + img_w - 1);
            for (int n = 0; n < NCH_RGB; ++n)
            {
                padded_img[i1+n] = padded_img[i2+n];
            }
        }
    }
}

/* Trim image of size src_w x src_h to the size of tgt_w x tgt_h
 *
 * Assumes src_w/h >= tgt_w/h, otherwise bad stuff will happen
 */
void trim_image(
    const std::vector<uint8_t>& src_img,
    int src_w,
    std::vector<uint8_t>& tgt_img,
    int tgt_w,
    int tgt_h
){
    // Just copy rows of tgt_w pixels
    for (int y = 0; y < tgt_h; ++y)
    {
        memcpy(
            tgt_img.data() + (y * tgt_w * NCH_RGB),  // dest
            src_img.data() + (y * src_w * NCH_RGB),  // src
            tgt_w * NCH_RGB                          // num_bytes
        );
    }
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
    // Input block size
    const int block_w = (ENC_FORMAT == ASTC) ? astc::BLOCK_X : 4;
    const int block_h = (ENC_FORMAT == ASTC) ? astc::BLOCK_Y : 4;

    // Number of blocks in horizontal and vertical dimension
    const int nblocks_x = img_w / block_w;
    const int nblocks_y = img_h / block_h;

    void (*encode_block)( const uint8_t*, uint32_t* );

    // Based on encoding format:
    //   1. Set encoded block size as a number of 32-bit integers
    //   2. Assign decoding function
    int block_nints = 0;
    switch (ENC_FORMAT) {
    case BC1:
        block_nints = 2;  // 64 bits per 4x4 block
        encode_block = bc1::encode_block;
        break;
    case YCOCG_BC3:
        block_nints = 4; // 128 bits per 4x4 block
        encode_block = ycocg_bc3::encode_block;
        break;
    case ASTC:
        block_nints = 4; // 128 bits per block
        encode_block = astc::encode_block;
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
            uint8_t block_pixels[NCH_RGB*block_w*block_h];
            const int x = block_x * block_w;
            const int y = block_y * block_h;
            for (int i = 0; i < block_h; ++i)
            {
                memcpy(
                    block_pixels + (NCH_RGB * i * block_w),
                    inp_pixels + (NCH_RGB * (img_w * (y+i) + x)),
                    NCH_RGB * block_w
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
        decode_block = bc1::decode_block;
        break;
    case YCOCG_BC3:
        block_nints = 4; // 128 bits per 4x4 block
        decode_block = ycocg_bc3::decode_block;
        break;
    case ASTC:
        return 0;
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
    return (double)tv.tv_sec + (double)tv.tv_usec * (double)1.0e-6;
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

    // Init & print out format-specific info
    if (ENC_FORMAT == ASTC)
    {
        astc::init_astc(12, 12, 8, 5);
        printf("WARNING: ASTC format decoding is not supported. Instead,"
               " encoded images are saved as .astc files in the output"
               " directory.\n");
    }

    if (ENC_FORMAT == BC1)
    {
        printf("INFO: BC1 is compiled with SELECT_DIAG=%d\n", BC1_SELECT_DIAG);
    }

    // Error code
    int err = 0;

    double total_enc_duration = 0.0;
    double total_duration = 0.0;
    uint64_t total_num_enc_pixels = 0;
    int num_enc_images = 0;

    // Loop through images one by one
    for (int i = 1; i < argc-1; ++i)
    {
        double start_time = get_time();

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

        // Pad image to fit the image dimensions being divisible by block size
        std::vector<uint8_t> padded_img;
        int pad_w = -1;
        int pad_h = -1;
        pad_image(inp_pixels, inp_w, inp_h, padded_img, pad_w, pad_h);

        if ( (inp_w != pad_w) || (inp_h != pad_h) )
        {
            printf("-- Image size %dx%d, after padding %dx%d\n", inp_w, inp_h, pad_w, pad_h);
        }

        double enc_start_time = get_time();

        // Encode it into enc_data
        std::vector<uint32_t> enc_data;
        if (encode_image(padded_img.data(), pad_w, pad_h, enc_data))
        {
            printf("-- Error encoding image\n");
            stbi_image_free(inp_pixels);
            continue;
        }

        total_enc_duration += ( get_time() - enc_start_time );
        total_num_enc_pixels += ( pad_w * pad_h );
        num_enc_images += 1;

        // Dump encoded data to file as raw bits
        // std::string dump_name = out_dir
        //     / (fs::path(inp_name).filename().replace_extension(".bin"));
        // dump_enc_data(enc_data, dump_name);

        // Decode encoded data
        std::vector<uint8_t> dec_image(pad_w*pad_h*NCH_RGB);
        if (decode_image(enc_data, pad_w, pad_h, dec_image))
        {
            printf("-- Error decoding image\n");
            stbi_image_free(inp_pixels);
            continue;
        }

        // Save result file into out_dir
        if (ENC_FORMAT == ASTC)
        {
            std::string out_name = out_dir
                / fs::path(inp_name).filename().replace_extension(".astc");
            printf("-- Saving encoded image to '%s'\n", out_name.data());
            int ret = astc::store_astc_image(
                (uint8_t*)enc_data.data(),
                4*enc_data.size(),
                inp_w,
                inp_h,
                out_name.c_str()
            );
            if (ret != 0)
            {
                printf("-- Error saving .astc file\n");
                stbi_image_free(inp_pixels);
                continue;
            }
        }
        else
        {
            // Trim the decoded image back into the original size before padding
            std::vector<uint8_t> out_pixels(inp_w*inp_h*NCH_RGB);
            trim_image(dec_image, pad_w, out_pixels, inp_w, inp_h);

            std::string out_name = out_dir
                / fs::path(inp_name).filename().replace_extension(".png");
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
        }

        stbi_image_free(inp_pixels);

        total_duration += ( get_time() - start_time );
    }

    printf("\n");
    printf("Encoded images                : %9d\n", num_enc_images);
    printf("Average encoding time (sec)   : %9.5f\n", total_enc_duration / num_enc_images);
    printf("Average encoding rate (Mpx/s) : %9.5f\n", total_num_enc_pixels / total_enc_duration / (double)1e6);
    printf("Average total time (sec)      : %9.5f\n", total_duration / num_enc_images);

    return err;
}
