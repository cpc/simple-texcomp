#!/usr/bin/env bash
#
# Just a quick and dirty script that builds the JPEG encoder for PC

clang++ \
    ../../../../../../src/transcoder.cpp \
    ../../../../../../src/simple_astc_common.cpp \
    ../../../../../../src/simple_astc_int_x86.cpp \
    ./tjenc.cpp \
    -I. \
    -I../../../../../../src \
    -Ilibturbojpeg/include \
    -I../../../../../../extern/stb \
    -I../../../../../../extern/tracy \
    -I../../../../../../extern/ghc \
    -Llibturbojpeg/x86_64 \
    -l:libturbojpeg.a \
    -std=c++17 \
    -Ofast \
    -ffast-math \
    -march=native \
    -o tjenc
