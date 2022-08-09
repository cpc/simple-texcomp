# Simple Texcomp

Simple CPU implementations of various texture compression encoders:

* A direct translation of the [reference real-time GPU implementation by Waveren and Castaño](https://developer.download.nvidia.com/whitepapers/2007/Real-Time-YCoCg-DXT-Compression/Real-Time%20YCoCg-DXT%20Compression.pdf).
  * `src/simple_bc1.h`: Implementation of the BC1 algorithm (low quality, 6:1 compression)
  * `src/simple_ycocg_bc3.h`: Implementation of the YCoCg-BC3 algorithm (high quality, 3:1 compression)
* Various implementations of [ASTC format](https://en.wikipedia.org/wiki/Adaptive_scalable_texture_compression) encoders with heavily pruned encoding configurations
  * `src/simple_astc.cpp`: Floating point implementation of the pruned ASTC encoder
  * `src/simple_astc_int_<target>.cpp`: Integer implementations fo the pruned ASTC encoder, where `<target>` denotes the target platform:
    * `arm`: Developed for armv8a ISA (Samsung S10 smartphone)
    * `x86`: For a general PC / laptop
    * `tce`: For [TCE](https://github.com/cpc/tce)
* RGB <-> YCoCg color space conversion
  * `src/simple_ycocg.cpp`

The whole project is wrapped under the `src/simple_texcomp.hpp` header as a library.
A `transcoder` command line application is built using this library (`src/transcoder_*` files) and can be used for encoding (in some cases also decoding) groups of images.

## Building

This project uses meson build system, clang++ compiler and requires C++17.

## Android App

There is an app for testing the integer ASTC encoder on a smartphone called "TexcompApp" under `android/`.
This app also includes a standalone JPEG encoder using the turbojpeg library.
The JPEG encoder can also be separately compiled as a command line utility (search for `build_tjenc.sh`).
