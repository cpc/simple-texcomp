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
    'b_ndebug=if-release',
  ]
)

add_project_arguments(
  ['-Wdouble-promotion', '-Wfloat-conversion', '-Wnarrowing'],
  language: 'cpp'
)

if get_option('bc1_sel')
  add_project_arguments('-DBC1_SELECT_DIAG=1', language: 'cpp')
else
  add_project_arguments('-DBC1_SELECT_DIAG=0', language: 'cpp')
endif

if get_option('use_double')
  add_project_arguments('-DUSE_DOUBLE=1', language: 'cpp')
else
  add_project_arguments('-DUSE_DOUBLE=0', language: 'cpp')
endif

includes = include_directories('src')
stb_includes = include_directories('extern/stb', is_system: true)

src_downsampler = [
  'src/simple_bilinear.cpp',
  'src/utils/downsampler.cpp',
]

src_transcoder = [
  'src/transcoder.cpp',
  'src/simple_bilinear.cpp',
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
