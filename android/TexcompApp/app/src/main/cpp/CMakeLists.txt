# For more information about using CMake with Android Studio, read the
# documentation: https://d.android.com/studio/projects/add-native-code.html

# Sets the minimum version of CMake required to build the native library.

cmake_minimum_required(VERSION 3.10.2)

# Declares and names the project.

project("texcompapp")

# poor man's build config
set(COMMON_FLAGS "-Wextra -ffast-math -std=c++17")
set(COMMON_FLAGS "${COMMON_FLAGS} -DBC1_SELECT_DIAG=0")
set(COMMON_FLAGS "${COMMON_FLAGS} -DASTC_SELECT_DIAG=0")
set(COMMON_FLAGS "${COMMON_FLAGS} -DASTC_TRIM_ENDPOINTS=0")
set(FLOAT_PRECISION "16")
set(COMMON_FLAGS "${COMMON_FLAGS} -DFLOAT_PRECISION=${FLOAT_PRECISION}")  # 16, 32 or 64
set(COMMON_FLAGS "${COMMON_FLAGS} -DENC_FORMAT_DEF=ASTC_INT")
#set(COMMON_FLAGS "${COMMON_FLAGS} -DTRACY_ENABLE")

if (${FLOAT_PRECISION} EQUAL "16")
    set(COMMON_FLAGS "${COMMON_FLAGS} -march=armv8.4-a+simd+sve+fp16")
else()
    set(COMMON_FLAGS "${COMMON_FLAGS} -march=armv8.4-a+simd+sve")
endif()

set(SCALAR_FLAGS "-fno-slp-vectorize -fno-vectorize -fno-tree-vectorize")
set(VECTOR_FLAGS "-Rpass=loop -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize  -fsave-optimization-record")

set(DEBUG_FLAGS "-DDEBUG -g")
set(RELEASE_FLAGS "-DNDEBUG -Ofast")

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)  # enable LTO

set(CMAKE_CXX_FLAGS "-DANDROID")
set(CMAKE_CXX_FLAGS_DEBUG "${COMMON_FLAGS} ${VECTOR_FLAGS} ${RELEASE_FLAGS}")

add_library(stb STATIC
        ../../../../../../extern/stb/stb_image.h
        ../../../../../../extern/stb/stb_image_write.h
        )
set_target_properties(stb PROPERTIES LINKER_LANGUAGE C)

add_library(ghc STATIC
        ../../../../../../extern/ghc/filesystem.hpp
        )
set_target_properties(ghc PROPERTIES LINKER_LANGUAGE CXX)

# Creates and names a library, sets it as either STATIC
# or SHARED, and provides the relative paths to its source code.
# You can define multiple libraries, and CMake builds them for you.
# Gradle automatically packages shared libraries with your APK.
add_library(
         # Sets the name of the library.
         native-lib

         # Sets the library as a shared library.
         SHARED

         # Provides a relative path to your source file(s).
         native-lib.cpp
         tjenc.cpp
)
target_include_directories(native-lib PUBLIC libturbojpeg/include)
target_include_directories(native-lib PUBLIC ../../../../../../extern/stb)
target_include_directories(native-lib PUBLIC ../../../../../../extern/ghc)
target_include_directories(native-lib PUBLIC ../../../../../../src)
set_target_properties(native-lib PROPERTIES LINKER_LANGUAGE CXX)

add_library(transcoder STATIC
        ../../../../../../src/platform.hpp
        ../../../../../../src/simple_mathlib.hpp
        ../../../../../../src/simple_texcomp.hpp
        ../../../../../../src/transcoder.hpp
        ../../../../../../src/simple_astc_int_platform.hpp

        ../../../../../../src/transcoder.cpp
        ../../../../../../src/simple_bilinear.cpp
        ../../../../../../src/simple_astc_common.cpp
        ../../../../../../src/simple_astc.cpp
        ../../../../../../src/simple_astc_int_arm.cpp
        ../../../../../../src/simple_bc1.cpp
        ../../../../../../src/simple_ycocg_bc3.cpp
        ../../../../../../src/simple_ycocg.cpp

        ../../../../../../extern/tracy/TracyClient.cpp
)
target_include_directories(transcoder PUBLIC ../../../../../../src)
target_include_directories(transcoder PUBLIC ../../../../../../extern/stb)
target_include_directories(transcoder PUBLIC ../../../../../../extern/ghc)
target_include_directories(transcoder PUBLIC ../../../../../../extern/tracy)
set_target_properties(transcoder PROPERTIES LINKER_LANGUAGE CXX)

# Searches for a specified prebuilt library and stores the path as a
# variable. Because CMake includes system libraries in the search path by
# default, you only need to specify the name of the public NDK library
# you want to add. CMake verifies that the library exists before
# completing its build.

find_library( # Sets the name of the path variable.
              log-lib

              # Specifies the name of the NDK library that
              # you want CMake to locate.
              log )

message("Android ABI: ${CMAKE_ANDROID_ARCH_ABI}")

#find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("Found OpenMP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    target_link_libraries(
            native-lib

            stb
            ghc

            ${log-lib}

            transcoder
            ${CMAKE_SOURCE_DIR}/libturbojpeg/${CMAKE_ANDROID_ARCH_ABI}/libturbojpeg.a

            OpenMP::OpenMP_CXX
    )
else()
    message(WARNING "OpenMP not found!")
    target_link_libraries(
            native-lib

            stb
            ghc

            ${log-lib}
            ${CMAKE_SOURCE_DIR}/libturbojpeg/${CMAKE_ANDROID_ARCH_ABI}/libturbojpeg.a

            transcoder
    )
endif()
