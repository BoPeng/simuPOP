# file SConstruct
#
# A scons based build system for simuPOP. It is better for development
# usage because it is a real multi-thread, dependency-based build system
# that can build part of the simuPOP modules..
#
#
import os, sys
import SCons.Defaults
import SCons.Tool

# do not update version for this development version to avoid rebuild
SIMUPOP_REV = '9990'
SIMUPOP_VER = '9.9.9'
all_modu = ['std', 'op', 'la', 'laop', 'ba', 'baop']
mpi_modu = ['mpi', 'opmpi', 'lampi', 'laopmpi', 'bampi', 'baopmpi']

# load all the module information from setup.py
from setup import *

# get information from python distutils
import distutils.sysconfig, os
vars = distutils.sysconfig.get_config_vars('CC', 'CXX', 'OPT', 'BASECFLAGS',
    'CCSHARED', 'LDSHARED', 'SO', 'LIBDEST')
for i in range(len(vars)):
    if vars[i] is None:
        vars[i] = ""
(cc, cxx, opt, basicflags, ccshared, ldshared, so_ext, lib_dest) = vars

env = Environment(ENV={'PATH':os.environ['PATH']},
    tools=['default', 'swig'])
# try to use the compiler used to build python
if cc != "":
    env['CC'] = cc
if cxx != "":
    env['CXX'] = cxx
if ldshared != '':
    env['SHLINK'] = ldshared

dest_dir  = os.path.join(lib_dest, 'site-packages')
build_dir = 'build'
env.BuildDir(build_dir, 'src', duplicate = 0)
env['build_dir'] = build_dir

# for swig 1.3.30, but do not use -outdir src, use -outdir . instead
# because the scons/swig module requires .py in the current directory (a bug, I would say).
env['SWIGFLAGS'] = SWIG_FLAGS.replace('-outdir src', '-outdir .')
env.Command('build_dir/swigpyrun.h', None, ['swig %s $TARGET' % SWIG_RUNTIME_FLAGS.replace('-outdir src', '-outdir .')])

gsl = env.StaticLibrary(
    target = '$build_dir/gsl',
    source = GSL_FILES,
    CCFLAGS = ' '.join([opt, ccshared]),
    CPPPATH = ['.'],
)

targets = []
for key in all_modu + mpi_modu:
    if key in BUILD_TARGETS:
        targets.append(key)
if targets == [] or 'all' in BUILD_TARGETS:
    # by default, do not build mpi version
    targets = all_modu

def mod_src(file, mod):
    return file.replace('src', '$build_dir').replace('.cpp', '_%s.cpp' % mod)

for mod in targets:
    info = ModuInfo(mod)
    for file in SOURCE_FILES:
        env.Command(mod_src(file, mod), file, [Copy('$TARGET', '$SOURCE')])
    lib = env.SharedLibrary(
        target = '$build_dir/_simuPOP_%s%s' % (mod, so_ext),
        source = ['$build_dir/simuPOP_%s.i' % mod] + [mod_src(x, mod) for x in SOURCE_FILES],
        LIBS = info['libraries'] + [gsl],
        SHLIBPREFIX = "",
        SHLIBSUFFIX = so_ext,
        LIBPATH = info['library_dirs'],
        CPPPATH = [distutils.sysconfig.get_python_inc(), '.', 'src'] + info['include_dirs'],
        CPPDEFINES = info['define_macros'],
        CCFLAGS = info['extra_compile_args'],
        CPPFLAGS = ' '.join([basicflags, ccshared, opt])
    )
    env.Depends(['$build_dir/simuPOP_%s_wrap$CXXFILESUFFIX' % mod, lib],
        ['src/simupop_cfg.h', 'src/simuPOP_common.i', 'src/simuPOP_doc.i'] + \
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



