#
# $File: SConstruct $
# $LastChangedDate: 2009-02-24 15:50:38 -0600 (Tue, 24 Feb 2009) $
# $Rev: 2491 $
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit http://simupop.sourceforge.net
# for details.
#
# Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

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
#      std op la laop ba baop mu muop: individual module
#      all = all modules
#      default to all modules
#      install: install specified targets
#      std_exe op_exe ....: individual executable for debugging purposes.
#      executable: std_exe, op_exe etc
#
#
import os, sys
import distutils.sysconfig

# load all the module information from setup.py
from setup import *
# do not update version for this development version to avoid rebuild
SIMUPOP_VER, SIMUPOP_REV = simuPOP_version()
print 'Using swig version %d.%d.%d' % tuple(swig_version())
# get information from python distutils
vars = distutils.sysconfig.get_config_vars('CC', 'CXX', 'OPT', 'BASECFLAGS',
    'CCSHARED', 'LDSHARED', 'SO', 'LIBDEST', 'INCLUDEPY')
for i in range(len(vars)):
    if vars[i] is None:
        vars[i] = ""
(cc, cxx, opt, basicflags, ccshared, ldshared, so_ext, lib_dest, python_inc_dir) = vars
# C++ does not need this option, remove it to avoid annoying warning messages.
opt = opt.replace('-Wstrict-prototypes', '').replace('-Wall', '')

# get compiler and hack options
from distutils.ccompiler import new_compiler
comp = new_compiler()
if comp.__dict__.has_key('initialized'):
    comp.initialize()
if not comp.__dict__.has_key('ldflags_shared'):
    comp.ldflags_shared = ''
if not comp.__dict__.has_key('compile_options'):
    comp.compile_options = []

opts = Variables()
opts.AddVariables(
    PathVariable('prefix', 'Where to install. see "python setup.py install --prefix"', None),
    PathVariable('include-dirs', 'Extra include directories, see "python setup.py build_ext --help"', None),
    PathVariable('library-dirs', 'Extra library directories, see "python setup.py build_ext --help"', None),
)

# common environment
env = Environment(
    options=opts,
    # pass all environment variables because MSVC needs INCLUDE and LIB
    # but they may not exist on other platforms
    ENV=os.environ,
    tools=['default', 'swig'])

# try to use the compiler used to build python
if cc != '':
    env['CC'] = cc
if cxx != '':
    env['CXX'] = cxx
if ldshared != '':
    env['SHLINK'] = ldshared

if env.has_key('prefix') and env['prefix'] is not None:
    pylib_dir = distutils.sysconfig.get_python_lib(plat_specific=1, prefix=env['prefix'])
    dest_dir = os.path.join(pylib_dir, 'simuPOP')
    prefix = env['prefix']
    print "Installing to", dest_dir
else:
    pylib_dir = distutils.sysconfig.get_python_lib(plat_specific=1, prefix=None)
    dest_dir = os.path.join(pylib_dir, 'simuPOP')
    prefix = sys.prefix

if os.name == 'nt':
    wrapper_flags = []
    extra_path = [os.path.join(sys.exec_prefix, 'libs'), 'win32']
    # treat warning as error. This is only set for development builds (scons)
    # because setup.py builds might be built by imcompatible version of compilers
    extra_flags = ['/Zi', '/WX']
else:
    wrapper_flags = ['-Wno-unused-parameter', '-Wno-missing-field-initializers']
    extra_path = [os.path.join(os.path.split(lib_dest)[:-1])]
    # during development, output symbols (for profiling and debugging)
    # treat warning as error... (although an unused parameter warning from boost has to be allowed)
    extra_flags = ['-g', '-Wno-unused-parameter', '-Wno-char-subscripts', '-Wno-unused-variable', 
        '-Wno-unused-function']
    #
    # GNU gcc:
    #
    # -Wconversion will produce a lot of warnings about type conversion. It is
    # useful but I cannot always turn it on because many warnings from boost will
    # overwhelm the output, and stop compiling if '-Werror' is defined. Because
    # MSVC does a reasonable job in detecting conversion problems and we have /WX
    # turned on for MSVC, we can turn on this option (and disable Werror) from
    # time to time to pick out warnings not detected by MSVC.
    #
    # Intel C++:
    #
    # -Wcheck and -Weffc++ give more warnings about coding style etc. They
    # also should not be turned on by default.
    #
    # -vec-report2 reports loops that are automatically vectorized, and remarks
    # for loops that cannot be vectorized.
    #
    #if USE_ICC:
    #    extra_flags.extend(['-Wcheck', '-Weffc++', '-vec-report2'])
    #else:
    #    extra_flags.append('-Wconversion')
    #

def convert_def(defines):
    new_list = []
    for d in defines:
        if d[1] is not None:
            new_list.append(d)
        else:
            new_list.append(d[0])
    return new_list
#
# Building library gsl
#
gsl_env = env.Clone()
gsl_env.VariantDir('build/gsl', '.')
gsl_env['SWIGOUTDIR'] = 'build/gsl/src'
gsl_env['SWIGFLAGS'] = SWIG_CC_FLAGS 
gsl = gsl_env.SharedLibrary(
    target = 'build/gsl/_gsl%s' % so_ext,
    source = ['build/gsl/' + x for x in GSL_FILES] + ['build/gsl/src/gsl.i'],
    SHLIBPREFIX = "",
    SHLIBSUFFIX = so_ext,
    SHLINKFLAGS = comp.ldflags_shared,
    LIBPATH = extra_path,
    # build for config.h
    CPPPATH = ['.', 'gsl', 'gsl/cdf', 'gsl/specfunc', 'build', python_inc_dir],
    CCFLAGS = comp.compile_options,
    CPPFLAGS = ' '.join([basicflags, ccshared, opt])
)
Alias('gsl', gsl)
Alias('all', gsl)
Alias('install', gsl_env.InstallAs(os.path.join(dest_dir, '_gsl%s' % so_ext), gsl[0]))
Alias('install', gsl_env.InstallAs(os.path.join(dest_dir, 'gsl.py'), 'build/gsl/src/gsl.py'))
#
#
# Building a library for common files
#
# Note that we suppress warnings when building these external source files
# using -wd options
std_common_env = env.Clone()
std_common_env.VariantDir('build/std_common', '.')
std_common_lib = std_common_env.StaticLibrary(
    target = 'build/std_common/common',
    source = ['build/std_common/' + x for x in LIB_FILES],
    CCFLAGS = ModuInfo('std', SIMUPOP_VER, SIMUPOP_REV)['extra_compile_args'] +\
        comp.compile_options + (['-wd1125', '-wd2259'] if USE_ICC else []),
    CPPPATH = ['.', 'gsl', 'build', ModuInfo('std', SIMUPOP_VER, SIMUPOP_REV)['include_dirs']],
    CPPDEFINES = convert_def(ModuInfo('std', SIMUPOP_VER, SIMUPOP_REV)['define_macros']),
    CPPFLAGS = ccshared + ' ' + opt,
)
op_common_env = env.Clone()
op_common_env.VariantDir('build/op_common', '.')
op_common_lib = op_common_env.StaticLibrary(
    target = 'build/op_common/common',
    source = ['build/op_common/' + x for x in LIB_FILES],
    CCFLAGS = ModuInfo('op', SIMUPOP_VER, SIMUPOP_REV)['extra_compile_args'] + \
        comp.compile_options + (['-wd1125'] if USE_ICC else []),
    CPPPATH = ['.', 'gsl', 'build', ModuInfo('op', SIMUPOP_VER, SIMUPOP_REV)['include_dirs']],
    CPPDEFINES = convert_def(ModuInfo('op', SIMUPOP_VER, SIMUPOP_REV)['define_macros']),
    CPPFLAGS = ccshared + ' ' + opt,
)
Alias('common', (std_common_lib, op_common_lib))

#
# Building modules
# 
targets = []
exe_targets = []
all_modu = ['std', 'op', 'la', 'laop', 'ba', 'baop', 'mu', 'muop', 'lin', 'linop']
all_exe = [x + '_exe' for x in all_modu]
for key in all_modu:
    if key in BUILD_TARGETS:
        targets.append(key)
for key in all_exe:
    if key in BUILD_TARGETS:
        exe_targets.append(key)

if targets == [] and 'gsl' not in BUILD_TARGETS and len(exe_targets) == 0:
    targets = all_modu
if 'all' in BUILD_TARGETS:
    targets = all_modu

for mod in targets:
    mod_env = env.Clone()
    mod_env.VariantDir('build/' + mod, '.')
    mod_env['SWIGFLAGS'] = SWIG_CPP_FLAGS + ' -Isrc'  # -Isrc for %include interface files under src
    mod_env['SWIGOUTDIR'] = 'build/%s/src' % mod
    info = ModuInfo(mod, SIMUPOP_VER, SIMUPOP_REV)
    mod_env.Command('build/%s/src/swigpyrun.h' % mod, None, ['swig %s $TARGET' % SWIG_RUNTIME_FLAGS])
    if 'op' in mod:
        common_lib = op_common_lib
    else:
        common_lib = std_common_lib
    mod_wrapper = mod_env.StaticLibrary(
        target = 'build/%s/src/_simuPOP_wrap_%s' % (mod, mod),
        source = ['build/%s/src/simuPOP_%s.i' % (mod, mod)],
        SHLIBPREFIX = "",
        SHLIBSUFFIX = so_ext,
        SHLINKFLAGS = comp.ldflags_shared,
        LIBPATH = info['library_dirs'] + extra_path,
        CPPPATH = [python_inc_dir, '.', 'src', 'build'] + info['include_dirs'],
        CPPDEFINES = convert_def(info['define_macros']),
        # No extra-flags because it can have -Werror and we cannot guarantee that there is
        # no warning coming out of the wrapper code.
        CCFLAGS = opt.split() + info['extra_compile_args'] + comp.compile_options + wrapper_flags,
        CPPFLAGS = ' '.join([basicflags, ccshared])
    )
    mod_lib = mod_env.SharedLibrary(
        target = 'build/%s/_simuPOP_%s' % (mod, mod),
        source = ['build/%s/src/%s' % (mod, x) for x in SOURCE_FILES],
        LIBS = info['libraries'] + [mod_wrapper, common_lib],
        SHLIBPREFIX = "",
        SHLIBSUFFIX = so_ext,
        SHLINKFLAGS = comp.ldflags_shared,
        LIBPATH = info['library_dirs'] + extra_path,
        CPPPATH = [python_inc_dir, '.', 'src', 'build'] + info['include_dirs'],
        CPPDEFINES = convert_def(info['define_macros']),
        CCFLAGS = info['extra_compile_args'] + comp.compile_options + extra_flags,
        CPPFLAGS = ' '.join([basicflags, ccshared, opt])
    )
    env.Depends('build/%s/src/swigpyrun.h' % mod, 'build/%s/src/utility.cpp' % mod)
    Alias(mod, mod_lib)
    Alias('all', mod_lib)
    pyInstall = env.InstallAs(os.path.join(dest_dir, 'simuPOP_%s.py' % mod),
        'build/%s/src/simuPOP_%s.py' % (mod, mod))
    libInstall = env.InstallAs(os.path.join(dest_dir, '_simuPOP_%s%s' % (mod, so_ext)),
        mod_lib[0])
    mod_env.Depends(pyInstall, libInstall)
    Alias('install', libInstall)
    Alias('install', pyInstall)

# module executable, for testing purposes.
if os.name == 'nt':
    pylibs = ['python%d%d' % (sys.version_info[0], sys.version_info[1])]
else:
    pylibs = ['python%d.%d' % (sys.version_info[0], sys.version_info[1]), 'util']

for mod in exe_targets:
    mod_name = mod.split('_')[0]
    exe_env = env.Clone()
    exe_env.VariantDir('build/' + mod, '.')
    exe_env['SWIGFLAGS'] = SWIG_CPP_FLAGS + ' -Isrc'  # -Isrc for %include interface files under src
    exe_env['SWIGOUTDIR'] = 'build/%s/src' % mod
    info = ModuInfo(mod_name, SIMUPOP_VER, SIMUPOP_REV)
    if 'op' in mod:
        common_lib = op_common_lib
    else:
        common_lib = std_common_lib
    mod_executable = exe_env.Program(
        target = 'build/simuPOP_%s' % mod,
        source = ['build/%s/src/%s' % (mod, x) for x in SOURCE_FILES],
        LIBS = info['libraries'] + [common_lib] + pylibs,
        LIBPATH = info['library_dirs'] + extra_path,
        CPPPATH = [python_inc_dir, '.', 'src', 'build'] + info['include_dirs'],
        CPPDEFINES = convert_def(info['define_macros'] + [('STANDALONE_EXECUTABLE', None)]),
        CCFLAGS = info['extra_compile_args'] + comp.compile_options + extra_flags,
        CPPFLAGS = ' '.join([basicflags, ccshared, opt])
    )
    Alias(mod, mod_executable)
    Alias('debug', mod_executable)


env.Install(pylib_dir, 'simuOpt.py')
Alias('install', pylib_dir)
for pyfile in ['__init__.py', 'utils.py', 'plotter.py', 'sampling.py', 'sandbox.py', 'demography.py']:
    env.Install(dest_dir, 'src/%s' % pyfile)
    Alias('install', dest_dir)

for datafile in PACKAGE_DATA:
    env.Install(dest_dir, 'src/%s' % datafile)
    Alias('install', dest_dir)

Default('install')
