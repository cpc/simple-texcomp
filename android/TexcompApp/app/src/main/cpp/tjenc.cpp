#include "tjenc.h"

#include <cstring>

#include <filesystem.hpp>
namespace fs = ghc::filesystem;

#include <platform.hpp>
#include <transcoder.hpp>
#include <stb_image.h>
#include <stb_image_write.h>

#include "turbojpeg.h"
#define THROW(action, message) { \
  LOGI("ERROR in line %d while %s:\n%s\n", __LINE__, action, message); \
  retval = -1;  goto bailout; \
}
#define THROW_TJ(action)  THROW(action, tjGetErrorStr2(tjInstance))
#define THROW_UNIX(action)  THROW(action, strerror(errno))
#define DEFAULT_SUBSAMP  TJSAMP_420
//#define DEFAULT_QUALITY  45

static inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

int compress_jpeg_image(
    const std::string &inp_name,
    const std::string &out_name,
    int q,
    int *img_w,
    int *img_h,
    double *enc_time
) {
    const char *subsampName[TJ_NUMSAMP] = {
        "4:4:4", "4:2:2", "4:2:0", "Grayscale", "4:4:0", "4:1:1"
    };

    const char *colorspaceName[TJ_NUMCS] = {
        "RGB", "YCbCr", "GRAY", "CMYK", "YCCK"
    };

    const char *pixelFormatName[TJ_NUMPF] = {
        "RGB", "BGR", "RGBX", "BGRX", "XBGR", "XRGB", "GRAY", "RGBA", "BGRA", "ABGR", "ARGB", "CMYK"
    };

    tjscalingfactor *scalingFactors = NULL;
    int numScalingFactors = 0;

    tjscalingfactor scalingFactor = {1, 1};
    int outSubsamp = -1, outQual = -1;
    tjtransform xform;
    int flags = 0;
    int width, height;
    FILE *jpegFile = NULL;
    unsigned char *imgBuf = NULL, *jpegBuf = NULL;
    int retval = 0, i, pixelFormat = TJPF_UNKNOWN;
    tjhandle tjInstance = NULL;

    std::string inFormat = fs::path(inp_name).extension();  // argv[1]
    std::string outFormat = fs::path(out_name).extension(); // argv[2]

    #if 0
    if (!strcasecmp(inFormat.c_str(), ".jpg") || !strcasecmp(inFormat.c_str(), ".jpeg")) {
        /* Input image is a JPEG image.  Decompress and/or transform it. */
        long size;
        int inSubsamp, inColorspace;
        int doTransform = (xform.op != TJXOP_NONE || xform.options != 0 ||
                           xform.customFilter != NULL);
        unsigned long jpegSize;

        /* Read the JPEG file into memory. */
        if ((jpegFile = fopen(inp_name.c_str(), "rb")) == NULL) THROW_UNIX("opening input file");
        if (fseek(jpegFile, 0, SEEK_END) < 0 || ((size = ftell(jpegFile)) < 0) ||
            fseek(jpegFile, 0, SEEK_SET) < 0) THROW_UNIX("determining input file size");
        if (size == 0) THROW("determining input file size", "Input file contains no data");
        jpegSize = (unsigned long) size;
        if ((jpegBuf = (unsigned char *) tjAlloc(jpegSize)) == NULL) THROW_UNIX(
                "allocating JPEG buffer");
        if (fread(jpegBuf, jpegSize, 1, jpegFile) < 1) THROW_UNIX("reading input file");
        fclose(jpegFile);
        jpegFile = NULL;

        if (doTransform) {
            LOGI("Doing transform\n");
            /* Transform it. */
            unsigned char *dstBuf = NULL;  /* Dynamically allocate the JPEG buffer */
            unsigned long dstSize = 0;

            if ((tjInstance = tjInitTransform()) == NULL) THROW_TJ("initializing transformer");
            xform.options |= TJXOPT_TRIM;
            if (tjTransform(tjInstance, jpegBuf, jpegSize, 1, &dstBuf, &dstSize,
                            &xform, flags) < 0) {
                tjFree(dstBuf);
                THROW_TJ("transforming input image");
            }
            tjFree(jpegBuf);
            jpegBuf = dstBuf;
            jpegSize = dstSize;
        } else {
            if ((tjInstance = tjInitDecompress()) == NULL) THROW_TJ("initializing decompressor");
        }

        if (tjDecompressHeader3(tjInstance, jpegBuf, jpegSize, &width, &height,
                                &inSubsamp, &inColorspace) < 0) THROW_TJ("reading JPEG header");

        LOGI("%s Image:  %d x %d pixels, %s subsampling, %s colorspace\n",
             (doTransform ? "Transformed" : "Input"), width, height,
             subsampName[inSubsamp], colorspaceName[inColorspace]);

        if ((!strcasecmp(outFormat.c_str(), ".jpg") || !strcasecmp(outFormat.c_str(), ".jpeg")) &&
            doTransform &&
            scalingFactor.num == 1 && scalingFactor.denom == 1 && outSubsamp < 0 &&
            outQual < 0) {
            /* Input image has been transformed, and no re-compression options
               have been selected.  Write the transformed image to disk and exit. */
            if ((jpegFile = fopen(out_name.c_str(), "wb")) == NULL) THROW_UNIX(
                    "opening output file");
            if (fwrite(jpegBuf, jpegSize, 1, jpegFile) < 1) THROW_UNIX("writing output file");
            fclose(jpegFile);
            jpegFile = NULL;
            goto bailout;
        }

        /* Scaling and/or a non-JPEG output image format and/or compression options
           have been selected, so we need to decompress the input/transformed
           image. */
        width = TJSCALED(width, scalingFactor);
        height = TJSCALED(height, scalingFactor);
        if (outSubsamp < 0)
            outSubsamp = inSubsamp;

        pixelFormat = TJPF_BGRX;
        if ((imgBuf = (unsigned char *) tjAlloc(width * height *
                                                tjPixelSize[pixelFormat])) == NULL) THROW_UNIX(
                "allocating uncompressed image buffer");

        if (tjDecompress2(tjInstance, jpegBuf, jpegSize, imgBuf, width, 0, height,
                          pixelFormat, flags) < 0) THROW_TJ("decompressing JPEG image");
        tjFree(jpegBuf);
        jpegBuf = NULL;
        tjDestroy(tjInstance);
        tjInstance = NULL;
    } else {
        /* Input image is not a JPEG image.  Load it into memory. */
        if ((imgBuf = tjLoadImage(inp_name.c_str(), &width, 1, &height, &pixelFormat,
                                  0)) == NULL) THROW_TJ("loading input image");
        if (outSubsamp < 0) {
            if (pixelFormat == TJPF_GRAY)
                outSubsamp = TJSAMP_GRAY;
            else
                outSubsamp = TJSAMP_444;
        }
        // LOGI("Input Image (%s):  %d x %d pixels, pixel format: %s\n",
        //      inp_name.c_str(), width, height, pixelFormatName[pixelFormat]);
    }
    #endif

    // Read input image
    int nch;
    imgBuf = stbi_load(
        inp_name.data(),
        &width,
        &height,
        &nch,
        3
    );

    if (imgBuf == NULL)
    {
        LOGI("-- Error opening %s, skipping\n", inp_name.data());
        return -1 ;
    }

    if (nch == 1) {
        pixelFormat = TJPF_GRAY;
    } else {
        pixelFormat = TJPF_RGB;
    }

    LOGI("Input Image (%s):  %d x %d pixels, pixel format: %s\n",
         inp_name.c_str(), width, height, pixelFormatName[pixelFormat]);
    LOGI("Output Image (%s):  %d x %d pixels\n", out_name.c_str(), width, height);

    *img_w = width;
    *img_h = height;

    if ( !strcasecmp(outFormat.c_str(), ".jpg") || !strcasecmp(outFormat.c_str(), ".jpeg") ) {
        /* Output image format is JPEG.  Compress the uncompressed image. */
        unsigned long jpegSize = 0;

        jpegBuf = NULL;  /* Dynamically allocate the JPEG buffer */

        outSubsamp = DEFAULT_SUBSAMP;
        outQual = q;
        LOGI("-- %s subsampling, quality = %d\n", subsampName[outSubsamp], outQual);

        if ((tjInstance = tjInitCompress()) == NULL)
            THROW_TJ("initializing compressor");

        double enc_start_time = simple::get_time();
        if (tjCompress2(tjInstance, imgBuf, width, 0, height, pixelFormat,
                        &jpegBuf, &jpegSize, outSubsamp, outQual, flags) < 0)
            THROW_TJ("compressing image");
        double enc_end_time = simple::get_time();

        tjDestroy(tjInstance);  tjInstance = NULL;

        *enc_time = enc_end_time - enc_start_time;

        /* Write the JPEG image to disk. */
        if ((jpegFile = fopen(out_name.c_str(), "wb")) == NULL)
            THROW_UNIX("opening output file");
        if (fwrite(jpegBuf, jpegSize, 1, jpegFile) < 1)
            THROW_UNIX("writing output file");
        tjDestroy(tjInstance);  tjInstance = NULL;
        fclose(jpegFile);  jpegFile = NULL;
        tjFree(jpegBuf);  jpegBuf = NULL;
    } else {
        /* Output image format is not JPEG.  Save the uncompressed image
           directly to disk. */
        LOGI("\n");
        if (tjSaveImage(out_name.c_str(), imgBuf, width, 0, height, pixelFormat, 0) < 0)
            THROW_TJ("saving output image");
    }

bailout:
    tjFree(imgBuf);
    if (tjInstance) tjDestroy(tjInstance);
    tjFree(jpegBuf);
    if (jpegFile) fclose(jpegFile);
    return retval;
}

int decompress_jpeg_image(
    const std::string &inp_name,
    const std::string &out_name,
    int *img_w,
    int *img_h,
    double *dec_time
) {
    const char *subsampName[TJ_NUMSAMP] = {
        "4:4:4", "4:2:2", "4:2:0", "Grayscale", "4:4:0", "4:1:1"
    };

    const char *colorspaceName[TJ_NUMCS] = {
        "RGB", "YCbCr", "GRAY", "CMYK", "YCCK"
    };

    const char *pixelFormatName[TJ_NUMPF] = {
        "RGB", "BGR", "RGBX", "BGRX", "XBGR", "XRGB", "GRAY", "RGBA", "BGRA", "ABGR", "ARGB", "CMYK"
    };

    tjscalingfactor *scalingFactors = NULL;
    int numScalingFactors = 0;

    tjscalingfactor scalingFactor = {1, 1};
    int outSubsamp = -1, outQual = -1;
    tjtransform xform;
    int flags = 0;
    int width, height;
    FILE *jpegFile = NULL;
    unsigned char *imgBuf = NULL, *jpegBuf = NULL;
    int retval = 0, i, pixelFormat = TJPF_UNKNOWN;
    tjhandle tjInstance = NULL;

    std::string inFormat = fs::path(inp_name).extension();  // argv[1]
    // std::string outFormat = fs::path(out_name).extension(); // argv[2]

    if (!strcasecmp(inFormat.c_str(), ".jpg") || !strcasecmp(inFormat.c_str(), ".jpeg")) {
        /* Input image is a JPEG image.  Decompress and/or transform it. */
        long size;
        int inSubsamp, inColorspace;
        int doTransform = (xform.op != TJXOP_NONE || xform.options != 0 ||
                           xform.customFilter != NULL);
        unsigned long jpegSize;

        if ((scalingFactors = tjGetScalingFactors(&numScalingFactors)) == NULL)
            THROW_TJ("getting scaling factors");
        memset(&xform, 0, sizeof(tjtransform));

        /* Read the JPEG file into memory. */
        if ((jpegFile = fopen(inp_name.c_str(), "rb")) == NULL) THROW_UNIX("opening input file");
        if (fseek(jpegFile, 0, SEEK_END) < 0 || ((size = ftell(jpegFile)) < 0) ||
            fseek(jpegFile, 0, SEEK_SET) < 0) THROW_UNIX("determining input file size");
        if (size == 0) THROW("determining input file size", "Input file contains no data");
        jpegSize = (unsigned long) size;
        if ((jpegBuf = (unsigned char *) tjAlloc(jpegSize)) == NULL) THROW_UNIX(
                "allocating JPEG buffer");
        if (fread(jpegBuf, jpegSize, 1, jpegFile) < 1) THROW_UNIX("reading input file");
        fclose(jpegFile);
        jpegFile = NULL;

        double tmp_time = 0.0;
        if ((tjInstance = tjInitDecompress()) == NULL) THROW_TJ("initializing decompressor");

        double dec_start_time = simple::get_time();
        if (tjDecompressHeader3(tjInstance, jpegBuf, jpegSize, &width, &height,
                                &inSubsamp, &inColorspace) < 0) THROW_TJ("reading JPEG header");
        tmp_time += simple::get_time() - dec_start_time;

        // if ((!strcasecmp(outFormat.c_str(), ".jpg") || !strcasecmp(outFormat.c_str(), ".jpeg")) &&
        //     doTransform &&
        //     scalingFactor.num == 1 && scalingFactor.denom == 1 && outSubsamp < 0 &&
        //     outQual < 0) {
        //     [> Input image has been transformed, and no re-compression options
        //        have been selected.  Write the transformed image to disk and exit. */
        //     if ((jpegFile = fopen(out_name.c_str(), "wb")) == NULL) THROW_UNIX(
        //             "opening output file");
        //     if (fwrite(jpegBuf, jpegSize, 1, jpegFile) < 1) THROW_UNIX("writing output file");
        //     fclose(jpegFile);
        //     jpegFile = NULL;
        //     goto bailout;
        // }

        /* Scaling and/or a non-JPEG output image format and/or compression options
           have been selected, so we need to decompress the input/transformed
           image. */
        width = TJSCALED(width, scalingFactor);
        height = TJSCALED(height, scalingFactor);
        if (outSubsamp < 0)
            outSubsamp = inSubsamp;

        // pixelFormat = TJPF_BGRX;
        pixelFormat = TJPF_RGB;
        if ((imgBuf = (unsigned char *) tjAlloc(width * height *
                                                tjPixelSize[pixelFormat])) == NULL) THROW_UNIX(
                "allocating uncompressed image buffer");

        dec_start_time = simple::get_time();
        if (tjDecompress2(tjInstance, jpegBuf, jpegSize, imgBuf, width, 0, height,
                          pixelFormat, flags) < 0) THROW_TJ("decompressing JPEG image");
        double dec_end_time = simple::get_time();
        tmp_time += simple::get_time() - dec_start_time;
        *dec_time = tmp_time;

        LOGI("%s Image:  %d x %d pixels, %s subsampling, %s colorspace\n",
             (doTransform ? "Transformed" : "Input"), width, height,
             subsampName[inSubsamp], colorspaceName[inColorspace]);

        LOGI("-- Saving decoded image to '%s'\n", out_name.data());
        int ret = stbi_write_png(
            out_name.data(),
            width,
            height,
            3,
            imgBuf,
            width*3
        );
        if (ret == 0)
        {
            // err_exit("Can't save output image");
            LOGE("Can't save output image\n");
            return 1;
        }

        tjFree(jpegBuf);
        jpegBuf = NULL;
        tjDestroy(tjInstance);
        tjInstance = NULL;
    } else {
        /* Input image is not a JPEG image. Error. */
        LOGI("Not a JPEG image: %s\n", inp_name.c_str());
    }

    // LOGI("Decomp %s -> %s", inp_name.c_str(), out_name.c_str());
    // *enc_time = 0.0;

bailout:
    *img_w = width;
    *img_h = height;
    tjFree(imgBuf);
    if (tjInstance) tjDestroy(tjInstance);
    tjFree(jpegBuf);
    if (jpegFile) fclose(jpegFile);
    return retval;
}

int compress_jpeg_images(const std::vector<std::string> &inp_images, const std::string &out_dir, int q)
{
    int ret = 0;

    double total_enc_duration = 0.0;
    double total_duration = 0.0;
    uint64_t total_num_enc_pixels = 0;
    int num_enc_images = 0;

    for (auto& inp_name : inp_images)
    {
        double start_time = simple::get_time();

        std::string out_name = out_dir / fs::path(inp_name).filename().replace_extension(".jpg");

        int img_w, img_h;
        double enc_time;
        ret = compress_jpeg_image(inp_name, out_name, q, &img_w, &img_h, &enc_time);

        total_enc_duration += enc_time;
        total_num_enc_pixels += ( img_w * img_h );
        ++num_enc_images;

        double enc_tp = img_w * img_h / enc_time / 1E6;
        LOGI("Encoding time: %.3f ms  (%12.6f Mpx/s)  [%s]\n", enc_time * 1E3, enc_tp, out_name.c_str());

        total_duration += ( simple::get_time() - start_time );
    }

    LOGI("Encoded images                : %12d\n", num_enc_images);
    LOGI("Average encoding time (sec)   : %12.6f\n", total_enc_duration / num_enc_images);
    LOGI("Average encoding rate (Mpx/s) : %12.6f\n", total_num_enc_pixels / total_enc_duration / 1E6);
    LOGI("Average total time (sec)      : %12.6f\n", total_duration / num_enc_images);

    return ret;
}

int decompress_jpeg_images(const std::string &inp_dir)
{
    int ret = 0;

    const fs::directory_iterator inp_dir_iter(inp_dir);
    std::vector<std::string> inp_images;
    for (auto& f : inp_dir_iter) {
        auto fname = f.path().string();
        transform(fname.begin(), fname.end(), fname.begin(), ::tolower);

        if ( ends_with(fname, ".jpg") || ends_with(fname, ".jpeg") )
        {
            inp_images.push_back(f.path().string());
        }
    }

    if (inp_images.size() == 0)
    {
        LOGI("No JPEG images to encode");
        return -1;
    }

    const std::string out_dir = inp_dir + "_dec";
    fs::create_directories(out_dir);

    double total_dec_duration = 0.0;
    double total_duration = 0.0;
    uint64_t total_num_dec_pixels = 0;
    int num_dec_images = 0;

    for (auto& inp_name : inp_images)
    {
        double start_time = simple::get_time();

        std::string out_name = out_dir / fs::path(inp_name).filename().replace_extension(".png");

        int img_w, img_h;
        double dec_time;
        ret = decompress_jpeg_image(inp_name, out_name, &img_w, &img_h, &dec_time);

        total_dec_duration += dec_time;
        total_num_dec_pixels += ( img_w * img_h );
        ++num_dec_images;

        double dec_tp = img_w * img_h / dec_time / 1E6;
        LOGI("Decoding time: %.3f ms  (%12.6f Mpx/s)\n", dec_time * 1E3, dec_tp);

        total_duration += ( simple::get_time() - start_time );
    }


    LOGI("Decoded images                : %12d\n", num_dec_images);
    LOGI("Average decoding time (sec)   : %12.6f\n", total_dec_duration / num_dec_images);
    LOGI("Average decoding rate (Mpx/s) : %12.6f\n", total_num_dec_pixels / total_dec_duration / 1E6);
    LOGI("Average total time (sec)      : %12.6f\n", total_duration / num_dec_images);

    return ret;
}

/* Exit with error, optionally printing usage */
void err_exit(const std::string& err_msg, bool print_usage=false)
{
    fprintf(stderr, "ERROR: %s\n", err_msg.data());
    if (print_usage)
    {
        fprintf(stderr,
            "Usage: tjenc input_image [<input_image> ...] output_folder q_value\n"
        );
    }
    exit(1);
}

int main(int argc, char **argv)
{
    // Check command line arguments
    if (argc < 4)
    {
        err_exit(
            "Provide at least one image, an output folder and a Q parameter",
            true);
    }

    // Last argument is the quality parameter
    std::string q_str = argv[argc-1];
    int q = std::stoi(q_str);

    // Second last argument is the output directory
    std::string out_dir = argv[argc-2];
    if ( !(fs::is_directory(out_dir) && fs::exists(out_dir)) )
    {
        err_exit("Can't open output directory", true);
    }

    // First arguments are input files
    std::vector<std::string> inp_images(argc-3);
    for (int i = 0; i < argc-3; ++i)
    {
        inp_images.at(i) = std::string(argv[i+1]);
    }

    LOGI("Output directory: '%s'\n", out_dir.data());
    if (compress_jpeg_images(inp_images, out_dir, q) != 0)
    {
        err_exit("Error encoding images", false);
    }
}
