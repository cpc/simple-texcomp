# To change compiler, set CXX environment variable (to e.g. clang++ or g++)
# To generate compile_commands.json, run `bear -- meson compile ...` (useful
# for clang LSP)
project(
  'simple-texcomp',
  'cpp',
  default_options: [
    'cpp_std=c++17', # c++17 for std::filesystem
    'buildtype=release',
    'warning_level=3',
    'b_ndebug=if-release',
    'b_lto=true'
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

if get_option('astc_sel')
  add_project_arguments('-DASTC_SELECT_DIAG=1', language: 'cpp')
else
  add_project_arguments('-DASTC_SELECT_DIAG=0', language: 'cpp')
endif

if get_option('astc_trim_endpoints')
  add_project_arguments('-DASTC_TRIM_ENDPOINTS=1', language: 'cpp')
else
  add_project_arguments('-DASTC_TRIM_ENDPOINTS=0', language: 'cpp')
endif

if get_option('use_double')
  add_project_arguments('-DUSE_DOUBLE=1', language: 'cpp')
else
  add_project_arguments('-DUSE_DOUBLE=0', language: 'cpp')
endif

if get_option('fast_math')
  add_project_arguments('-ffast-math', language: 'cpp')
endif

scalar_args = [ '-fno-slp-vectorize', '-fno-vectorize',  '-fno-tree-vectorize', ]
vector_args = [ '-march=native' ]

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
  cpp_args: vector_args,
)

executable(
  'transcoder',
  src_transcoder,
  include_directories: [ includes, stb_includes ],
  cpp_args: vector_args,
)

executable(
  'transcoder-scalar',
  src_transcoder,
  include_directories: [ includes, stb_includes ],
  cpp_args: scalar_args,
)
