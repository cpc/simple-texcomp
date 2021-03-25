# To change compiler, set CXX environment variable (to e.g. clang++ or g++)
# To generate compile_commands.json, run `bear -- meson compile ...` (useful
# for clang LSP)
project(
  'simple-texcomp',
  'cpp',
  default_options: [
    'cpp_std=c++17', # due to std::filesystem, otherwise could be lower
    'buildtype=debugoptimized',
    'warning_level=3',
  ]
)

add_project_arguments(
  ['-Wdouble-promotion', '-Wfloat-conversion'],
  language: 'cpp'
)

includes = include_directories('src')
stb_includes = include_directories('extern/stb', is_system: true)

src_downsampler = [
  'src/utils/downsampler.cpp',
]

src_transcoder = [
  'src/transcoder.cpp',
  'src/simple_astc.cpp',
  'src/simple_bc1.cpp',
  'src/simple_ycocg_bc3.cpp',
]

executable(
  'downsampler',
  src_downsampler,
  include_directories: [ includes, stb_includes ],
)

executable(
  'transcoder',
  src_transcoder,
  include_directories: [ includes, stb_includes ],
)
