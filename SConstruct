# file SConstruct
#
# A scons based build system for simuPOP. It is better for development
# usage because it is a real multi-thread, dependency-based build system
# that can build part of the simuPOP modules..
#
# Usage:
#    scons [scons-options] [options] [targets]
# where:
#    scons-options:
#      standard scons options like -j (number of threads)
#
#    options:
#      prefix=/path/to/prefix:  
#          prefix of installation, equivalent to the --prefix option of python setup.py
#      include-dirs=/path/to/include:/path/to/includes2
#          extra include directories, usually to boost library. The path separator is ; under windows.
#      library-dirs=/path/to/lib:/path/to/lib2
#          extra library directories, usually to boost library. The path separator is ; under windows.
#      
#   targets: one of more of
#      std op la laop ba baop mpi opmpi lampi laopmpi bampi baopmpi: individual module
#      all = all modules
#      default to all standard modules (non-MPI version)
#      install: install specified targets
#
#
#  BUG: This script does *not* work under solaris because a needed change in _wrap.cpp files.
#       This change is handled in setup.py, but not in SConstruct
#
import os, sys

# do not update version for this development version to avoid rebuild
SIMUPOP_REV = '9990'
SIMUPOP_VER = '9.9.9'
all_modu = ['std', 'op', 'la', 'laop', 'ba', 'baop']
mpi_modu = ['mpi', 'opmpi', 'lampi', 'laopmpi', 'bampi', 'baopmpi']

if not os.path.isfile('SConstruct'):
    print 'Please run scons from simuPOP source directory.'
    Exit(1)

# Needs scons 0.96.93
from SCons import __version__
version = map(int, __version__.split('.'))
if version[0] != 0 or version[1] != 96 or version[2] != 93:
    print "Scons version 0.96.93 is required."
    Exit(1)

# load all the module information from setup.py
from setup import *

# get information from python distutils
import distutils.sysconfig, os
vars = distutils.sysconfig.get_config_vars('CC', 'CXX', 'OPT', 'BASECFLAGS',
    'CCSHARED', 'LDSHARED', 'SO', 'LIBDEST', 'INCLUDEPY')
for i in range(len(vars)):
    if vars[i] is None:
        vars[i] = ""
(cc, cxx, opt, basicflags, ccshared, ldshared, so_ext, lib_dest, python_inc_dir) = vars
#
# get compiler and hack options
from distutils.ccompiler import new_compiler
comp = new_compiler()
if comp.__dict__.has_key('initialized'):
    comp.initialize()
if not comp.__dict__.has_key('ldflags_shared'):
    comp.ldflags_shared = ''
if not comp.__dict__.has_key('compile_options'):
    comp.compile_options = []

opts = Options()
opts.AddOptions(
    PathOption('prefix', 'Where to install. see "python setup.py install --prefix"', None),
    PathOption('include-dirs', 'Extra include directories, see "python setup.py build_ext --help"', None),
    PathOption('library-dirs', 'Extra library directories, see "python setup.py build_ext --help"', None),
)
    
env = Environment(
    options=opts,
	# pass all environment variables because MSVC needs INCLUDE and LIB
	# but they may not exist on other platforms
    ENV=os.environ,
    tools=['default', 'swig'])
# try to use the compiler used to build python
if cc != "":
    env['CC'] = cc
if cxx != "":
    env['CXX'] = cxx
if ldshared != '':
    env['SHLINK'] = ldshared

if env.has_key('prefix') and env['prefix'] is not None:
    dest_dir  = distutils.sysconfig.get_python_lib(plat_specific=1, prefix=env['prefix'])
    print "Installing to", dest_dir
else:
    dest_dir  = os.path.join(lib_dest, 'site-packages')

if env.has_key('include-dirs') and env['include-dirs'] is not None:
    boost_inc_search_paths.extend(env['include-dirs'].split(os.pathsep))
if env.has_key('library-dirs') and env['library-dirs'] is not None:
    boost_lib_search_paths.extend(env['library-dirs'].split(os.pathsep))

build_dir = 'build'
env.BuildDir(build_dir, 'src', duplicate = 0)
env['build_dir'] = build_dir

if os.name == 'nt':
    extra_lib_path = [os.path.join(sys.exec_prefix, 'libs'), 'win32']
else:
    extra_lib_path = []

# for swig 1.3.30, but do not use -outdir src, use -outdir . instead
# because the scons/swig module requires .py in the current directory (a bug, I would say).
env['SWIGFLAGS'] = SWIG_FLAGS.replace('-outdir src', '-outdir .')
env.Command('$build_dir/swigpyrun.h', None, ['swig %s $TARGET' % SWIG_RUNTIME_FLAGS.replace('-outdir src', '-outdir .')])

# this lib may contain gsl, boost/iostreams and boost/serialization
extra_lib = env.StaticLibrary(
    target = '$build_dir/extra_libs',
    source = LIB_FILES,
    CCFLAGS = ' '.join([opt, ccshared]),
    CPPPATH = ['.', ModuInfo('std')['include_dirs']],
)

targets = []
for key in all_modu + mpi_modu:
    if key in BUILD_TARGETS:
        targets.append(key)
if targets == []:
    # by default, do not build mpi version
    targets = all_modu 
if 'all' in BUILD_TARGETS:
    # if all is specified, build all
    targets = all_modu + mpi_modu

def mod_src(file, mod):
    return file.replace('src', '$build_dir').replace('.cpp', '_%s.cpp' % mod)

def convert_def(defines):
    new_list = []
    for d in defines:
        if d[1] is not None:
            new_list.append(d)
        else:
            new_list.append(d[0])
    return new_list

for mod in targets:
    info = ModuInfo(mod)
    for file in SOURCE_FILES:
        env.Command(mod_src(file, mod), file, [Copy('$TARGET', '$SOURCE')])
	module_source = [mod_src(x, mod) for x in SOURCE_FILES]
    lib = env.SharedLibrary(
        target = '$build_dir/_simuPOP_%s%s' % (mod, so_ext),
        source = module_source + ['$build_dir/simuPOP_%s.i' % mod],
        LIBS = info['libraries'] + [extra_lib],
        SHLIBPREFIX = "",
        SHLIBSUFFIX = so_ext,
        SHLINKFLAGS = comp.ldflags_shared,
        LIBPATH = info['library_dirs'] + extra_lib_path,
        CPPPATH = [python_inc_dir, '.', 'src'] + info['include_dirs'],
        CPPDEFINES = convert_def(info['define_macros']),
        CCFLAGS = info['extra_compile_args'] + comp.compile_options,
        CPPFLAGS = ' '.join([basicflags, ccshared, opt])
    )
    env.Depends(module_source, '$build_dir/swigpyrun.h')
    env.Depends(['$build_dir/simuPOP_%s_wrap$CXXFILESUFFIX' % mod, lib],
        ['src/simuPOP_cfg.h', 'src/simuPOP_common.i', 'src/simuPOP_doc.i'] + \
        HEADER_FILES)
    env.Depends('src/utility_%s.cpp' % mod, 'src/arraymodule.c')
    #
    Alias(mod, lib)
    Alias('all', lib)
    dp1 = env.InstallAs(os.path.join(dest_dir, 'simuPOP_%s.py' % mod),
        'simuPOP_%s.py' % mod)
    dp2 = env.InstallAs(os.path.join(dest_dir, '_simuPOP_%s%s' % (mod, so_ext)),
        lib[0])
    env.Depends(dp1, dp2)
    Alias('install', dp1)


for pyfile in SIMUPOP_FILES:
    env.Install(dest_dir, 'src/%s.py' % pyfile)
    Alias('install', dest_dir)

# install to share directory, later.
Default('install')



