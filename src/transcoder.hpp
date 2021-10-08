#ifndef TRANSCODER_H
#define TRANSCODER_H

#include <vector>
#include <string>

namespace simple {

/* Possible encoding formats */
typedef enum enc_format_t {
    BC1,
    YCOCG_BC3,
    YCOCG,
    ASTC,
    ASTC_INT,
} enc_format_t;

/* Define encoding format here */
#ifndef ENC_FORMAT_DEF
#define ENC_FORMAT_DEF ASTC_INT
#endif
constexpr enc_format_t ENC_FORMAT = ENC_FORMAT_DEF;

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
);

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
);

/* Encode the whole image according to ENC_FORMAT
 *
 * Returns zero on success, non-zero on failure.
 */
int encode_image(
    const uint8_t *inp_pixels,
    int img_w,
    int img_h,
    std::vector<uint32_t>& enc_data
);

/* Decode the whole image according to ENC_FORMAT
 *
 * Returns zero on success, non-zero on failure.
 */
int decode_image(
    const std::vector<uint32_t>& enc_data,
    int img_w,
    int img_h,
    std::vector<uint8_t>& out_pixels
);

/* Dump raw encoded data into a file (for debugging) */
void dump_enc_data(
    const std::vector<uint32_t>& enc_data,
    const std::string& filename
);

/* Read the current time and return seconds (taken from astcenc) */
double get_time();

/* Transcoder entry point.
 *
 * Returns exit code (non-zero => fail)
 */
int transcoder_entry(
    std::vector<std::string> inp_images,
    std::string out_dir
);

} // namespace simple

#endif // TRANSCODER_H
