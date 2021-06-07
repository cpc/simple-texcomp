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

float_precision = get_option('float_precision')
if float_precision == 'double'
  add_project_arguments('-DFLOAT_PRECISION=64', language: 'cpp')
elif float_precision == 'half'
  add_project_arguments('-DFLOAT_PRECISION=16', language: 'cpp')
else
  add_project_arguments('-DFLOAT_PRECISION=32', language: 'cpp')
endif

if get_option('fast_math')
  add_project_arguments('-ffast-math', language: 'cpp')
endif

tracy_deps = [ ]
if get_option('tracy_enable')
  add_project_arguments('-DTRACY_ENABLE', language: 'cpp')
  compiler = meson.get_compiler('cpp')
  libdl = compiler.find_library('dl', required: true)
  libpthread = compiler.find_library('pthread', required: true)
  tracy_deps = [ libdl, libpthread ]
endif

scalar_args = [ '-fno-slp-vectorize', '-fno-vectorize',  '-fno-tree-vectorize', ]
vector_args = [ '-march=native' ]
if get_option('diagnostics')
  vector_args = [ vector_args, '-Rpass=loop', '-Rpass-missed=loop-vectorize', '-Rpass-analysis=loop-vectorize',  '-fsave-optimization-record' ]
endif

includes = include_directories('src')
stb_includes = include_directories('extern/stb', is_system: true)
ghc_includes = include_directories('extern/ghc', is_system: true)
tracy_includes = include_directories('extern/tracy', is_system: true)

src_downsampler = [
  'src/simple_bilinear.cpp',
  'src/utils/downsampler.cpp',
  'extern/tracy/TracyClient.cpp',
]

src_transcoder = [
  'src/transcoder_main.cpp',
  'src/transcoder.cpp',
  'src/simple_bilinear.cpp',
  'src/simple_astc.cpp',
  'src/simple_bc1.cpp',
  'src/simple_ycocg_bc3.cpp',
  'extern/tracy/TracyClient.cpp',
]

executable(
  'downsampler',
  src_downsampler,
  include_directories: [ includes, stb_includes, ghc_includes, tracy_includes ],
  cpp_args: vector_args,
  dependencies: tracy_deps,
)

executable(
  'transcoder',
  src_transcoder,
  include_directories: [ includes, stb_includes, ghc_includes, tracy_includes ],
  cpp_args: vector_args,
  dependencies: tracy_deps,
)

executable(
  'transcoder-scalar',
  src_transcoder,
  include_directories: [ includes, stb_includes, ghc_includes, tracy_includes ],
  cpp_args: scalar_args,
  dependencies: tracy_deps,
)
