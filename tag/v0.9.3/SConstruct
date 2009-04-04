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
#

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
#      std op la laop ba baop: individual module
#      all = all modules
#      default to all modules
#      install: install specified targets
#
#
import os, sys

all_modu = ['std', 'op', 'la', 'laop', 'ba', 'baop']

if not os.path.isfile('SConstruct'):
    print 'Please run scons from simuPOP source directory.'
    Exit(1)

# load all the module information from setup.py
from setup import *

# do not update version for this development version to avoid rebuild
SIMUPOP_VER, SIMUPOP_REV = simuPOP_version()

# get information from python distutils
import distutils.sysconfig, os
vars = distutils.sysconfig.get_config_vars('CC', 'CXX', 'OPT', 'BASECFLAGS',
    'CCSHARED', 'LDSHARED', 'SO', 'LIBDEST', 'INCLUDEPY')
for i in range(len(vars)):
    if vars[i] is None:
        vars[i] = ""
(cc, cxx, opt, basicflags, ccshared, ldshared, so_ext, lib_dest, python_inc_dir) = vars
# C++ does not need this option, remove it to avoid annoying warning messages.
opt = opt.replace('-Wstrict-prototypes', '')

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
	# force the use of MS Visual Studio .NET 2003 if under windows
	# This does not work right when there are multiple versions of MSVS.
	MSVS_VERSION = '7.1',
    #MSVS_IGNORE_IDE_PATHS = 1,
    tools=['default', 'swig'])
# try to use the compiler used to build python
if cc != "":
    env['CC'] = cc
if cxx != "":
    env['CXX'] = cxx
if ldshared != '':
    env['SHLINK'] = ldshared

# Users can set a few options in a local file scons_cfg.py
# Only prefix is supported now.
if os.path.isfile('scons_cfg.py'):
    import scons_cfg as cfg
    try:
        if not env.has_key('prefix'):
            print 'Setting prefix to %s from scons_cfg.py' % cfg.prefix
            env['prefix'] = cfg.prefix
    except:
        pass

if env.has_key('prefix') and env['prefix'] is not None:
    dest_dir = distutils.sysconfig.get_python_lib(plat_specific=1, prefix=env['prefix'])
    prefix = env['prefix']
    print "Installing to", dest_dir
else:
    dest_dir = os.path.join(lib_dest, 'site-packages')
    prefix = sys.prefix

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
env['SWIGFLAGS'] = SWIG_FLAGS # .replace('-outdir src', '-outdir .')
env['SWIGOUTDIR'] = '$build_dir'
env.Command('$build_dir/swigpyrun.h', None, ['swig %s $TARGET' % SWIG_RUNTIME_FLAGS])

# this lib may contain gsl, boost/iostreams and boost/serialization
extra_lib = env.StaticLibrary(
    target = '$build_dir/extra_libs',
    source = LIB_FILES,
    CCFLAGS = ModuInfo('std', SIMUPOP_VER, SIMUPOP_REV)['extra_compile_args'] + comp.compile_options,
    CPPPATH = ['.', ModuInfo('std', SIMUPOP_VER, SIMUPOP_REV)['include_dirs']],
	CPPFLAGS = ccshared + ' ' + opt,
)

targets = []
for key in all_modu:
    if key in BUILD_TARGETS:
        targets.append(key)
if targets == []:
    targets = all_modu
if 'all' in BUILD_TARGETS:
    targets = all_modu

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
    info = ModuInfo(mod, SIMUPOP_VER, SIMUPOP_REV)
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
        '$build_dir/simuPOP_%s.py' % mod)
    dp2 = env.InstallAs(os.path.join(dest_dir, '_simuPOP_%s%s' % (mod, so_ext)),
        lib[0])
    env.Depends(dp1, dp2)
    Alias('install', dp1)


for pyfile in SIMUPOP_FILES:
    env.Install(dest_dir, 'src/%s.py' % pyfile)
    Alias('install', dest_dir)

for data in DATA_FILES:
    dest = data[0]
    for file in data[1]:
        env.Install(os.path.join(prefix, dest), file)
        Alias('install', os.path.join(prefix, dest))

Default('install')
