# Simple BCn

Simple scalar CPU implementation of BC1 and YCoCg-BC3 texture compression algorithms.
It is a direct translation of the [reference real-time GPU implementation by Waveren and Castaño](https://developer.download.nvidia.com/whitepapers/2007/Real-Time-YCoCg-DXT-Compression/Real-Time%20YCoCg-DXT%20Compression.pdf).

## Files & Notes

* `simple_bc1.h`: Implementation of the BC1 algorithm (low quality, 6:1 compression)
  * Setting `SELECT_DIAG` to `0` disables an optional refinement step and results in slightly lower quality but faster runtime
  * Functions for encoding/decoding one 4x4 pixel block
* `simple_ycocg_bc3.h`: Implementation of the YCoCg-BC3 algorithm (low quality, 3:1 compression)
  * Functions for encoding/decoding one 4x4 pixel block
* `simple_bcn_common.h`: Common definitions for the two files above
* `transcoder.cpp`: An example CLI app for encoding, decoding and saving the results
  * Select encoding algorithm with the `ENC_FORMAT` constant
  * Accepts images to be encoded as CLI arguments, the last argument is an output directory
    * Outputs decoded images into the output directory
    * (optional: uncomment if desired) Outputs raw encoded bytes into the output directory (with a `.bin` extension)
  * Example functions for encoding/decoding the whole image (must have width/height divisible by 4 due to the 4x4 block size)
* `Makefile`: Requires C++17 standard due to the use of `std::filesystem` for convenience in `transcoder.cpp`, `simple_*.h` should be C++11-compatible
* `stb/`: Image reading/writing (https://github.com/nothings/stb)

## Quality Comparison

The quality is slightly different between the CPU and GPU version due to the differences of floating point precision on both platforms.

### Average PSNR (dB) vs. uncompressed:

(decoded with https://github.com/ifeherva/bcndecode, not this repository)

| dataset |  algorithm |  CPU(this) | GPU(OpenCL) |
|:------- |:---------- | ----------:| -----------:|
| kodim   | bc1_nosel  | 34.9096    | 35.6300     |
| kodim   | bc1_sel    | 35.3600    | 35.6200     |
| kodim   | ycocg_bc3  | 41.1189    | 41.1392     |
| rgb8bit | bc1_nosel  | 37.6707    | 37.6845     |
| rgb8bit | bc1_sel    | 37.9579    | 37.9723     |
| rgb8bit | ycocg_bc3  | 42.3804    | 42.4597     |

* bc1_nosel : using `simple_bc1.h` with `#define SELECT_DIAG 0`
* bc1_sel   : using `simple_bc1.h` with `#define SELECT_DIAG 1`
* ycocg_bc3 : using `ycocg_bc3.h`

Datasets:
* kodim: http://r0k.us/graphics/kodak/index.html
* rgb8bit: http://imagecompression.info/test_images
