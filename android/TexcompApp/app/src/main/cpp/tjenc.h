#ifndef TJENC_H
#define TJENC_H

#include <string>
#include <vector>

int compress_jpeg_image(
    const std::string &inp_name,
    const std::string &out_name,
    int q,
    int *img_w,
    int *img_h,
    double *enc_time
);

int compress_jpeg_images(
    const std::vector<std::string> &inp_images,
    const std::string &out_dir,
    int q
);

int decompress_jpeg_image(
    const std::string &inp_name,
    const std::string &out_name,
    int *img_w,
    int *img_h,
    double *dec_time
);

int decompress_jpeg_images(
    const std::string &out_dir
);

#endif // TJENC_H
