#!/usr/bin/env python

#
# $File: setup.py $
# $LastChangedDate$
# $Rev$
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

"""
simuPOP installer: 

A standard simuPOP package has boost and needed SWIG-generated wrapper files so
you should be able to run this script and install simuPOP using:
  > python setup.py install

If your copy of simuPOP is checked out from the subversion server, you will
need to download a supported version of boost (see variable boost_versions
below) and uncompress it under the simuPOP source directory. You also need to
install swig >= 1.3.35 for the generation of Python wrapper files. Please see
http://simupop.sourceforge.net/main/GetInvolved for details.

"""
import os, sys, shutil, glob, re, tempfile

# simuPOP works with these boost versions. Newer versions will be used if these
# versions are not available, and will most likely work just fine.
boost_versions = ['1_35_0', '1_36_0', '1_37_0', '1_38_0', '1_39_0', '1_40_0']
invalid_boost_versions = ['1_41_0']

included_version = [x for x in boost_versions if os.path.isdir('boost_' + x)]
invalid_version = [x for x in invalid_boost_versions if os.path.isdir('boost_' + x)]
unsupported_version = [x for x in glob.glob('boost_*') if os.path.isdir(x) \
    and x[6:] not in invalid_version + included_version]

boost_dir = ''
if len(included_version) > 0:
    boost_dir = 'boost_' + included_version[-1]  # use the latest version
elif len(invalid_version) > 0:
    print 'This version of boost is known to cause problems to simuPOP: ' + ', '.join(invalid_version)

if boost_dir == '' and len(unsupported_version) > 0:
    print 'This version of boost is not tested. It may or may not work: ' + ', '.join(unsupported_version)
    boost_dir = unsupported_version[-1]  # use the latest version

if boost_dir == '':
    print 'Cannot find an useful version of boost header files.'
    print 'Please download boost from http://www.boost.org and unpack it under the simuPOP directory'
    sys.exit(1)
elif len(included_version + unsupported_version) > 1:
    print 'Using boost version: %s' % boost_dir

boost_include_dir = boost_dir
boost_serialization_dir = os.path.join(boost_dir, 'libs', 'serialization', 'src')
boost_iostreams_dir = os.path.join(boost_dir, 'libs', 'iostreams', 'src')
boost_regex_dir = os.path.join(boost_dir, 'libs', 'regex', 'src')

# if you need to use full path name for swig, change it here.
SWIG = 'swig'

############################################################################
#
# THE FOLLOWING IS NOT SUPPOSED TO BE MODIFIED
#
############################################################################

from distutils.core import setup, Extension
from distutils.sysconfig import get_config_var

def swig_version():
    ''' get the version of swig '''
    fout = os.popen(SWIG + ' -version')
    #
    try:
        version = re.match('SWIG Version\s*(\d+).(\d+).(\d+).*', fout.readlines()[1]).groups()
    except:
        print 'Can not obtain swig version, please install swig'
        sys.exit(1)
    return map(int, version)


def simuPOP_version():
    import simuPOP_version
    SIMUPOP_VER = simuPOP_version.SIMUPOP_VER
    SIMUPOP_REV = simuPOP_version.SIMUPOP_REV
    if SIMUPOP_VER.endswith('svn'):
        rev = SIMUPOP_REV
        try:
            rev = os.popen('svnversion .').readline().strip()
            if ':' in rev:
                rev = rev.split(':')[1]
            rev = rev.rstrip('M')
        except:
            pass
        # if 'svnversion' exists and the revision has changed
        if rev != '':
            SIMUPOP_REV = rev
    return SIMUPOP_VER, SIMUPOP_REV


def replaceIntHeader(file):
    ''' Replace "#include <stdint.h>" with "#include <inttypes.h>"
        Try to keep time stamp unchanged '''
    # save file modification time
    time = os.path.getmtime(file)
    # create a temp file
    tmp, name = tempfile.mkstemp()
    input = open(file)
    output = open(name, 'w')
    for line in input.readlines():
        if '#include <stdint.h>' in line:
            print >> output, '#include <inttypes.h>  /* no stdint.h is found so we use inttypes.h instead */'
        else:
            print >> output, line,
    output.close()
    input.close()
    # replace file with temp file
    shutil.copyfile(name, file)
    os.remove(name)    
    # restore old file modification time
    try:
        os.utime(file, (-1, time))
    except:
        pass

#
# SOURCE FILES
#

HEADER_FILES = [
    'src/simuPOP_cfg.h',
    'src/utility.h',
    'src/genoStru.h',
    'src/individual.h',
    'src/population.h',
    'src/simulator.h',
    'src/mating.h',
    'src/operator.h',
    'src/initializer.h',
    'src/migrator.h',
    'src/outputer.h',
    'src/selector.h',
    'src/penetrance.h',
    'src/qtrait.h',
    'src/stator.h',
    'src/mutator.h',
    'src/transmitter.h',
    'src/tagger.h',
    'src/pedigree.h',
    'src/virtualSubPop.h',
]

SOURCE_FILES = [
    'src/utility.cpp',
    'src/genoStru.cpp',
    'src/individual.cpp',
    'src/population.cpp',
    'src/simulator.cpp',
    'src/mating.cpp',
    'src/operator.cpp',
    'src/initializer.cpp',
    'src/migrator.cpp',
    'src/outputer.cpp',
    'src/selector.cpp',
    'src/penetrance.cpp',
    'src/qtrait.cpp',
    'src/stator.cpp',
    'src/mutator.cpp',
    'src/transmitter.cpp',
    'src/tagger.cpp',
    'src/pedigree.cpp',
    'src/virtualSubPop.cpp',
]

# since it is troublesome to link to external gsl library,
# I embed some GSL files with simuPOP. 
LIB_FILES = [ 
    'gsl/sys/infnan.c',
    'gsl/sys/coerce.c',
    'gsl/sys/fdiv.c',
    'gsl/sys/pow_int.c',
    'gsl/sys/fcmp.c',
    'gsl/sys/log1p.c',
    'gsl/sys/invhyp.c',
    'gsl/complex/math.c',
    'gsl/specfunc/psi.c',
    'gsl/specfunc/trig.c',
    'gsl/specfunc/exp.c',
    'gsl/specfunc/expint.c',
    'gsl/specfunc/log.c',
    'gsl/specfunc/erfc.c',
    'gsl/specfunc/zeta.c',
    'gsl/specfunc/elementary.c',
    'gsl/specfunc/gamma.c',
    'gsl/specfunc/gamma_inc.c',
    'gsl/rng/borosh13.c',
    'gsl/rng/fishman2x.c',
    'gsl/rng/mt.c',
    'gsl/rng/rand.c',
    'gsl/rng/ranmar.c',
    'gsl/rng/types.c',
    'gsl/rng/cmrg.c',
    'gsl/rng/gfsr4.c',
    'gsl/rng/r250.c',
    'gsl/rng/random.c',
    'gsl/rng/rng.c',
    'gsl/rng/uni32.c',
    'gsl/rng/coveyou.c',
    'gsl/rng/knuthran2.c',
    'gsl/rng/ran0.c',
    'gsl/rng/randu.c',
    'gsl/rng/slatec.c',
    'gsl/rng/uni.c',
    'gsl/rng/default.c',
    'gsl/rng/knuthran.c',
    'gsl/rng/ran1.c',
    'gsl/rng/ranf.c',
    'gsl/rng/taus113.c',
    'gsl/rng/vax.c',
    'gsl/rng/file.c',
    'gsl/rng/lecuyer21.c',
    'gsl/rng/ran2.c',
    'gsl/rng/ranlux.c',
    'gsl/rng/taus.c',
    'gsl/rng/waterman14.c',
    'gsl/rng/fishman18.c',
    'gsl/rng/minstd.c',
    'gsl/rng/ran3.c',
    'gsl/rng/ranlxd.c',
    'gsl/rng/transputer.c',
    'gsl/rng/zuf.c',
    'gsl/rng/fishman20.c',
    'gsl/rng/mrg.c',
    'gsl/rng/rand48.c',
    'gsl/rng/ranlxs.c',
    'gsl/rng/tt.c',
    'gsl/rng/knuthran2002.c',
    'gsl/randist/binomial_tpe.c',
    'gsl/randist/beta.c',
    'gsl/randist/exponential.c',
    'gsl/randist/geometric.c',
    'gsl/randist/nbinomial.c',
    'gsl/randist/poisson.c',
    'gsl/randist/multinomial.c',
    'gsl/randist/chisq.c',
    'gsl/randist/gauss.c',
    'gsl/randist/gausszig.c',
    'gsl/randist/gamma.c',
    'gsl/cdf/chisq.c',
    'gsl/cdf/gamma.c',
    'gsl/error.c' 
] + [x for x in glob.glob(os.path.join(boost_serialization_dir, '*.cpp')) if 'xml' not in x and 'binary' not in x]\
  + [x for x in glob.glob(os.path.join(boost_iostreams_dir, '*.cpp')) if 'bzip' not in x]\
  + glob.glob(os.path.join(boost_regex_dir, '*.cpp'))


GSL_FILES = [
    'gsl/error.c', 
    'gsl/sys/infnan.c',
    'gsl/sys/coerce.c',
    'gsl/sys/fdiv.c',
    'gsl/sys/pow_int.c',
    'gsl/sys/fcmp.c',
    'gsl/sys/log1p.c',
    'gsl/sys/invhyp.c',
    'gsl/sys/expm1.c',
    'gsl/complex/math.c',
    'gsl/specfunc/elementary.c',
    'gsl/specfunc/erfc.c',
    'gsl/specfunc/exp.c',
    'gsl/specfunc/expint.c',
    'gsl/specfunc/log.c',
    'gsl/specfunc/psi.c',
    'gsl/specfunc/gamma.c',
    'gsl/specfunc/gamma_inc.c',
    'gsl/specfunc/trig.c',
    'gsl/specfunc/zeta.c',
    'gsl/cdf/gauss.c',
    'gsl/cdf/gaussinv.c',
    'gsl/cdf/exponential.c',
    'gsl/cdf/exponentialinv.c',
    # function gsl_ran_gamma_pdf is copied to gsl.i to avoid subsequent
    # inclusion of RNG related functions.
    #'gsl/randist/gamma.c',
    'gsl/cdf/gamma.c',
    'gsl/cdf/gammainv.c',
    'gsl/cdf/chisq.c',
    'gsl/cdf/chisqinv.c',
]

# build zlib from source for windows system to avoid distributing zlib1.dll
# along with simuPOP.
if os.name == 'nt':
    LIB_FILES.extend([os.path.join('win32', 'zlib-1.2.3', x) for x in [
        'adler32.c',
        'compress.c',
        'crc32.c',
        'gzio.c',
        'uncompr.c',
        'deflate.c',
        'trees.c',
        'zutil.c',
        'inflate.c',
        'infback.c',
        'inftrees.c',
        'inffast.c'
        ]
    ])



SWIG_CPP_FLAGS = '-O -templatereduce -shadow -python -c++ -keyword -nodefaultctor -w-503,-312,-511,-362,-383,-384,-389,-315,-509,-525'
SWIG_CC_FLAGS = '-python -keyword'
SWIG_RUNTIME_FLAGS = '-python -external-runtime'
# python setup.py reads py_modules from src so we have to produce simuPOP_std.py
# etc to this directory.
SWIG_OUTDIR = 'src'
                
MACROS = {
    'std':    [('SIMUPOP_MODULE', 'simuPOP_std')],
    'op':     [('SIMUPOP_MODULE', 'simuPOP_op'), ('OPTIMIZED', None)],
    'la':     [('SIMUPOP_MODULE', 'simuPOP_la'), ('LONGALLELE', None)],
    'laop':   [('SIMUPOP_MODULE', 'simuPOP_laop'), ('LONGALLELE', None), ('OPTIMIZED', None)],
    'ba':     [('SIMUPOP_MODULE', 'simuPOP_ba'), ('BINARYALLELE', None) ],
    'baop':   [('SIMUPOP_MODULE', 'simuPOP_baop'), ('BINARYALLELE', None), ('OPTIMIZED', None)],
}
 
WRAP_INFO = {
    'std':    ['src/simuPOP_std_wrap.cpp', 'src/simuPOP_std.i', ''],
    'op':     ['src/simuPOP_op_wrap.cpp', 'src/simuPOP_op.i', '-DOPTIMIZED'],
    'la':     ['src/simuPOP_la_wrap.cpp', 'src/simuPOP_la.i', '-DLONGALLELE'],
    'laop':   ['src/simuPOP_laop_wrap.cpp', 'src/simuPOP_laop.i', '-DLONGALLELE -DOPTIMIZED'],
    'ba':     ['src/simuPOP_ba_wrap.cpp', 'src/simuPOP_ba.i', '-DBINARYALLELE'],
    'baop':   ['src/simuPOP_baop_wrap.cpp', 'src/simuPOP_baop.i', '-DBINARYALLELE -DOPTIMIZED'],
}

DESCRIPTION = """
simuPOP is a forward-time population genetics simulation environment.
The core of simuPOP is a scripting language (Python) that provides 
a large number of objects and functions to manipulate populations, 
and a mechanism to evolve populations forward in time. Using this 
R/Splus-like environment, users can create, manipulate and evolve 
populations interactively, or write a script and run it as a batch 
file. Owing to its flexible and extensible design, simuPOP can simulate
large and complex evolutionary processes with ease. At a more 
user-friendly level, simuPOP provides an increasing number of built-in
scripts that perform simulations ranging from implementation of basic 
population genetics models to generating datasets under complex 
evolutionary scenarios.
"""
            

#
DATA_FILES = [
    ('share/simuPOP', ['README', 'INSTALL', 'ChangeLog', 'AUTHORS', 
        'COPYING', 'TODO', 'simuPOP_version.py']), 
    ('share/simuPOP/doc', ['doc/userGuide.pdf', 'doc/userGuide.py', 'doc/refManual.pdf']), 
    ('share/simuPOP/test', glob.glob('test/test_*.py') + ['test/run_tests.py'])
]


def ModuInfo(modu, SIMUPOP_VER, SIMUPOP_REV):
    #
    boost_inc_path = boost_include_dir
    boost_lib_names = []
    boost_lib_path = None
    res = {}
    res['src'] =  ['src/simuPOP_' + modu + '_wrap.cpp']
    for src in SOURCE_FILES:
        res['src'].append(src[:-4] + '_' + modu + '.cpp')
    res['src'].extend(LIB_FILES)
    # lib
    if os.name == 'nt':    # Windows, build zlib from source
        res['libraries'] = []
    else:
        res['libraries'] = ['stdc++', 'z']
    res['libraries'].extend(boost_lib_names)
    res['include_dirs'] = ['.', boost_inc_path]
    res['library_dirs'] = ['build']
    if os.name == 'nt':
        # I have a portable stdint.h for msvc, to avoid distributing
        # zdll1.dll, I also build zlib from source
        res['include_dirs'].extend(['win32', 'win32/zlib-1.2.3'])
        # zdll.lib is under win32
        res['library_dirs'].append('win32')
    if os.name == 'nt':
        # msvc does not have O3 option, /GR is to fix a C4541 warning
        # /EHsc is for VC exception handling,
		# /wd4819 is used to disable a warning for non-unicode character in boost/uitlity/enable_if.hpp
        res['extra_compile_args'] = ['/O2', '/GR', '/EHsc', '/wd4819']
    else:
        res['extra_compile_args'] = ['-O3', '-Wall']
    # if Intel ICC is used, turn off remark 981
    if os.getenv('CC', '').endswith('icc') or \
        (get_config_var('CC') is not None and 'icc' in get_config_var('CC')):
        res['extra_compile_args'].append('-wd981')
    # define_macros (deep copy)
    res['define_macros'] = [x for x in MACROS[modu]]
    res['define_macros'].extend([('SIMUPOP_VER', SIMUPOP_VER), ('SIMUPOP_REV', SIMUPOP_REV)])
    if os.name == 'nt':
        res['define_macros'].extend([('BOOST_ALL_NO_LIB', None),
            ('NO_ZLIB', 0), ('NO_BZIP' , 1),
			# this one disables a lot of warnings about VC Checked iterators. Might not be a good idea.
			('_SCL_SECURE_NO_WARNINGS', None)
			])
    res['undef_macros'] = []
    return res


############################################################################
#
# Build extensions
#
############################################################################
# checking os type and copy configuration files
if os.name == 'nt':    # Windows
    shutil.copy('config_win32.h', 'config.h')
elif os.name == 'posix':
    if sys.platform == 'linux2':     # Linux
        shutil.copy('config_linux.h', 'config.h')
    elif sys.platform == 'sunos5': # Solaris
        shutil.copy('config_solaris.h', 'config.h')
    elif sys.platform == 'darwin':    # MacOS
        shutil.copy('config_mac.h', 'config.h')
    else: # HPUX?
        shutil.copy('config_linux.h', 'config.h')
else:
    # otherwise, assume a posix system
    shutil.copy('config_linux.h', 'config.h')

if __name__ == '__main__':
    SIMUPOP_VER, SIMUPOP_REV = simuPOP_version()
    # create source file for each module
    MODULES = ['std', 'op', 'la', 'laop', 'ba', 'baop']
    #
    # Generate Wrapping files
    #
    # if any of the wrap files does not exist
    # or if the wrap files are older than any of the source files.
    if (False in [os.path.isfile(WRAP_INFO[x][0]) for x in MODULES]) or \
        (max( [os.path.getmtime(x) for x in HEADER_FILES] ) > \
         min( [os.path.getmtime(WRAP_INFO[x][0]) for x in MODULES])):
        (v1, v2, v3) = swig_version()
        if (v1, v2, v3) < (1, 3, 35):
            print 'Swig >= 1.3.35 is required, please upgrade it.'
            sys.exit(1)
        # generate header file 
        print "Generating external runtime header file..."
        os.system('swig %s swigpyrun.h' % SWIG_RUNTIME_FLAGS)
        # try the first option set with the first library
        for lib in MODULES:
            print "Generating wrapper file " + WRAP_INFO[lib][0]
            if os.system('%s %s -outdir %s %s -o %s %s' % (SWIG, SWIG_CPP_FLAGS, \
                SWIG_OUTDIR, WRAP_INFO[lib][2], WRAP_INFO[lib][0], WRAP_INFO[lib][1])) != 0:
                print "Calling swig failed. Please check your swig version."
                sys.exit(1)
        print "Generating wrapper file src/gsl_wrap.cpp"
        if os.system('%s %s -outdir %s %s -o %s %s' % (SWIG, SWIG_CC_FLAGS, \
            SWIG_OUTDIR, '', 'src/gsl_wrap.cpp', 'src/gsl.i')) != 0:
            print "Calling swig failed. Please check your swig version."
            sys.exit(1)
        print
        print "All wrap files are generated successfully."
        print
    # under solaris, there is no stdint.h so I need to replace stdint.h
    # in the wrap files with inttypes.h
    if sys.platform == 'sunos5':
        for lib in MODULES:
            replaceIntHeader(WRAP_INFO[lib][0])
    # copy needed files
    copied_files = []
    for modu in MODULES:
        for src in SOURCE_FILES:
            mod_src = src[:-4] + '_' + modu + '.cpp'
            shutil.copy(src, mod_src)
            copied_files.append(mod_src)
    # build
    EXT_MODULES = []
    for modu in MODULES:
        info = ModuInfo(modu, SIMUPOP_VER=SIMUPOP_VER, SIMUPOP_REV=SIMUPOP_REV)
        EXT_MODULES.append(
            Extension('simuPOP._simuPOP_%s' % modu,
                sources = info['src'],
                extra_compile_args = info['extra_compile_args'],
                include_dirs = info['include_dirs'],
                library_dirs = info['library_dirs'],
                libraries = info['libraries'],
                define_macros = info['define_macros'],
                undef_macros = info['undef_macros'],
            )
        )
    # For module simuPOP.gsl
    EXT_MODULES.append(
        Extension('simuPOP.gsl',
            sources = GSL_FILES,
            extra_compile_args = info['extra_compile_args'],
            include_dirs = info['include_dirs'],
        )
    )
    setup(
        name = "simuPOP",
        version = SIMUPOP_VER,
        author = "Bo Peng",
        author_email = "bpeng@mdanderson.org",
        maintainer = "Bo Peng",
        maintainer_email = "bpeng@mdanderson.org",
        url = "http://simupop.sourceforge.net",
        description = "Forward-time population genetics simulation environment",
        long_description = DESCRIPTION, 
        download_url = 'http://simupop.sourceforge.net/Main/Download',
        classifiers = [
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Education',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Natural Language :: English',
            'Operating System :: OS Independent',
            'Programming Language :: C++',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        platforms = ['all'],
        #
        packages = ['simuPOP'],
        package_dir = {'simuPOP': 'src'}, 
        py_modules = [
            'simuOpt', 
            'simuPOP.__init__',
            'simuPOP.utils', 
            'simuPOP.plotter', 
            'simuPOP.sampling', 
        ] + ['simuPOP.simuPOP_%s' % x for x in MODULES],
        ext_modules = EXT_MODULES,
        data_files = DATA_FILES,
    )
    # remove copied files
    for file in copied_files:
        os.remove(file)


