# To change compiler, set CXX environment variable (to e.g. clang++ or g++)
# To generate compile_commands.json, run `bear -- meson compile ...` (useful
# for clang LSP)
project(
  'simple-texcomp',
  'cpp',
  default_options: [
    'cpp_std=c++17', # due to std::filesystem, otherwise could be lower
    'buildtype=release',
  ]
)

includes = include_directories('extern/stb')

executable('downsampler', 'src/utils/downsampler.cpp', include_directories: includes)
executable('quantizer', 'src/utils/quantizer.cpp', include_directories: includes)
executable('transcoder', 'src/transcoder.cpp', include_directories: includes)
