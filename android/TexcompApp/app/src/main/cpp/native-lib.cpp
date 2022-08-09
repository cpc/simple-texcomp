#include <jni.h>
#include <cstdint>
#include <string>
#include <algorithm>
#include <cstdio>
#include <fstream>

// std image library
//#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

// 3rd party std::filesystem equivalent
#include <filesystem.hpp>
namespace fs = ghc::filesystem;

// simple-texcomp
#include <transcoder.hpp>
#include <platform.hpp>

// TurboJPEG stuff
#include "tjenc.h"

// this is in C++20 but we use 17 here
static inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

extern "C" JNIEXPORT jstring JNICALL
Java_fi_tuni_cpc_texcompapp_MainActivity_getImgInfo(
        JNIEnv* env,
        jobject /* this */,
        jstring inp_name
) {
    const char* name = env->GetStringUTFChars(inp_name, nullptr);

    // Read input image
    int inp_w, inp_h, nch;
    uint8_t *inp_pixels = stbi_load(
            name,
            &inp_w,
            &inp_h,
            &nch,
            3
    );
    std::string ret = std::to_string(inp_w) + "x" + std::to_string(inp_h);
    if (inp_pixels == nullptr)
    {
        ret = "Could not open file: ";
        ret.append(name);
    }

    env->ReleaseStringUTFChars(inp_name, name);
    return env->NewStringUTF(ret.c_str());
}


static unsigned long long read_id_aa64isar0()
{
    unsigned long long id_aa64isar0;

    __asm ("MRS %x0, ID_AA64ISAR0_EL1 \n" : "=r" (id_aa64isar0) );

    return (id_aa64isar0);
}

static bool dot_product_supported()
{
    if (read_id_aa64isar0() & 0x0000100000000000ULL)
        return true;
    else
        return false;
}

extern "C" JNIEXPORT jstring JNICALL
Java_fi_tuni_cpc_texcompapp_MainActivity_encodeDirectory(
        JNIEnv* env,
        jobject /* this */,
        jstring inp_name,
        jstring format
) {
#ifdef TRACY_ENABLE
    LOGW("WARNING: Profiling enabled!");
#endif

    const std::string stem = fs::path(env->GetStringUTFChars(inp_name, nullptr)).filename();
    const std::string enc_format = env->GetStringUTFChars(format, nullptr);

    int q = 96; // TODO: Allow setting it manually

    std::string out_dir;
    if (enc_format == "jpeg") {
        out_dir = "/sdcard/texcomp/out_jpeg/" + stem;
    } else if (enc_format == "astc") {
        out_dir = "/sdcard/texcomp/out_astc/" + stem;
    } else {
        out_dir = "/sdcard/texcomp/out/" + stem;
    }
    fs::create_directories(out_dir);

    const fs::directory_iterator inp_dir(env->GetStringUTFChars(inp_name, nullptr));
    std::vector<std::string> inp_images;
    for (auto& f : inp_dir) {
        auto fname = f.path().string();
        transform(fname.begin(), fname.end(), fname.begin(), ::tolower);

        if ( ends_with(fname, ".png")
            || ends_with(fname, ".jpg")
            || ends_with(fname, ".jpeg")
            || ends_with(fname, ".bmp")
            || ends_with(fname, ".ppm") )
        {
            inp_images.push_back(f.path().string());
        }
    }

    if (dot_product_supported()) {
        LOGI("Dot product (UDOT instruction) SUPPORTED");
    } else {
        LOGI("Dot product (UDOT instruction) NOT SUPPORTED");
    }

    int ret = -1;
    if (enc_format == "jpeg") {
        LOGI("Encoding JPEG");
        ret = compress_jpeg_images(inp_images, out_dir, q);
        if (ret == 0) {
            ret = decompress_jpeg_images(out_dir);
        }
    } else if (enc_format == "astc") {
        LOGI("Encoding ASTC");
        ret = simple::transcoder_entry(inp_images, out_dir);
    } else {
        LOGI("Unknown encoding format");
    }
    return env->NewStringUTF(std::to_string(ret).c_str());
}
