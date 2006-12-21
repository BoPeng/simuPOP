"""
simuPOP installer 

In case that you have modified the C++ code, or you are
checking out simuPOP from svn repository, you need to 
install swig >= 1.3.27 to generate the wrap files.

"""
import os
#

# this is the toolset to build boost. It is needed because
# boost libraries have form libboost_xxx-TOOLSET.a
TOOLSET = 'gcc'

############################################################################
#
# THE FOLLOWING IS NOT SUPPOSED TO BE MODIFIED
#
############################################################################

from distutils.core import setup, Extension
from distutils.ccompiler import new_compiler

import shutil, sys, glob, re

# this is borrowed from pypar, since I do not have access to other
# MPI implementation, I use mpicc only.
def getMPIFlags():
    ''' get and parse result of mpiCC -showme or mpicc -show'''
    fin, fout, ferr = os.popen3('mpiCC -show')
    output = fout.read()
    fout.close()
    #
    if output == '':
        fin, fout, ferr = os.popen3('mpiCC -showme')
        output = fout.read()
        fout.close()
    if output == '':
        return {'mpi': False, 'inc_dirs': [], 'lib_dirs': [], 'libs': [],
            'def_macros': [], 'undef_macros': []}
    # now get the include, library dirs and the libs to link with.
    flags = output.split()
    inc_dirs = []
    lib_dirs = []
    libs = []
    def_macros = []
    undef_macros = []
    for f in flags:
        if f[:2] == '-I':
            inc_dirs.append(f[2:])
        elif f[:2] == '-L':
            lib_dirs.append(f[2:])
        elif f[:2] == '-l':
            libs.append(f[2:])
        elif f[:2] == '-U':
            undef_macros.append(f[2:])
        elif f[:2] == '-D':
            tmp = string.split(f[2:], '=')
            if len(tmp) == 1:
                def_macros.append((tmp[0], None))
            else:
                def_macros.append(tuple(tmp))
    return {'mpi': True, 'inc_dirs': inc_dirs, 'lib_dirs': lib_dirs, 'libs':libs,
            'def_macros': def_macros, 'undef_macros': undef_macros}


def swig_version():
    ''' get the version of swig '''
    fout = os.popen('swig -version')
    #
    try:
        version = re.match('SWIG Version\s*(\d+).(\d+).(\d+).*', fout.readlines()[1]).groups()
    except:
        print 'Can not obtain swig version, please install swig'
        sys.exit(1)
    return map(int, version)


############################################################################
#
# SOURCE FILES
#
############################################################################

HEADER_FILES = [
    'src/simupop_cfg.h',
    'src/utility.h',
    'src/individual.h',
    'src/population.h',
    'src/simulator.h',
    'src/mating.h',
    'src/operator.h',
    'src/initializer.h',
    'src/migrator.h',
    'src/outputer.h',
    'src/selector.h',
    'src/stator.h',
    'src/terminator.h',
    'src/mutator.h',
    'src/recombinator.h',
    'src/tagger.h',
]

SOURCE_FILES = [
    'src/utility.cpp',
    'src/individual.cpp',
    'src/population.cpp',
    'src/simulator.cpp',
    'src/mating.cpp',
    'src/operator.cpp',
    'src/initializer.cpp',
    'src/migrator.cpp',
    'src/outputer.cpp',
    'src/selector.cpp',
    'src/stator.cpp',
    'src/terminator.cpp',
    'src/mutator.cpp',
    'src/recombinator.cpp',
    'src/tagger.cpp',
]

# since it is troublesome to link to external gsl library,
# I embed some GSL files with simuPOP. 
GSL_FILES = [ 
    'gsl/sys/infnan.c',
    'gsl/sys/coerce.c',
    'gsl/sys/fdiv.c',
    'gsl/sys/pow_int.c',
    'gsl/sys/fcmp.c',
    'gsl/specfunc/psi.c',
    'gsl/specfunc/trig.c',
    'gsl/specfunc/exp.c',
    'gsl/specfunc/log.c',
    'gsl/specfunc/erfc.c',
    'gsl/specfunc/zeta.c',
    'gsl/specfunc/elementary.c',
    'gsl/specfunc/gamma.c',
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
    'gsl/randist/nbinomial.c',
    'gsl/randist/beta.c',
    'gsl/randist/exponential.c',
    'gsl/randist/geometric.c',
    'gsl/randist/binomial.c',
    'gsl/randist/poisson.c',
    'gsl/randist/rdgamma.c',
    'gsl/randist/multinomial.c',
    'gsl/randist/chisq.c',
    'gsl/randist/gauss.c',
    'gsl/error.c' 
]


SIMUPOP_FILES = [
    'simuPOP', 
    'simuOpt', 
    'simuUtil', 
    'simuRPy', 
    'simuViewPop'
]

############################################################################
#
# MODULE SETTINGS
#
############################################################################

# for every official release, there will be a file recording release info
# Othersise, SIMUPOP_VER and SIMUPOP_REV will be provided as environmental
# variables.
if os.environ.has_key('SIMUPOP_VER') and os.environ.has_key('SIMUPOP_REV'):
        SIMUPOP_VER = os.environ['SIMUPOP_VER']
        SIMUPOP_REV = os.environ['SIMUPOP_REV']
else:
        execfile('simuPOP.release')

# explore availability of mpi library
MPIFlags = getMPIFlags()

# create source file for each module
MODULES = ['std', 'op', 'la', 'laop', 'ba', 'baop']
#if MPIFlags['mpi']:
if False:
    # feel like compiling mpi version first :-)
    MODULES = ['mpi', 'opmpi', 'lampi', 'laopmpi', 'bampi', 'baopmpi'] + MODULES
MODU_INFO = {}

SIMUPOP_FILES += ['simuPOP_%s' % x for x in MODULES]

#
MACROS = {
    'std':    [('SIMUPOP_MODULE', 'simuPOP_std')],
    'op':     [('SIMUPOP_MODULE', 'simuPOP_op'), ('OPTIMIZED', None)],
    'la':     [('SIMUPOP_MODULE', 'simuPOP_la'), ('LONGALLELE', None)],
    'laop':   [('SIMUPOP_MODULE', 'simuPOP_laop'), ('LONGALLELE', None), ('OPTIMIZED', None)],
    'ba':     [('SIMUPOP_MODULE', 'simuPOP_ba'), ('BINARYALLELE', None) ],
    'baop':   [('SIMUPOP_MODULE', 'simuPOP_baop'), ('BINARYALLELE', None), ('OPTIMIZED', None)],
    'mpi':    [('SIMUPOP_MODULE', 'simuPOP_mpi'), ('SIMUMPI', None)] + MPIFlags['def_macros'],
    'opmpi':  [('SIMUPOP_MODULE', 'simuPOP_opmpi'), ('OPTIMIZED', None), ('SIMUMPI', None)] + MPIFlags['def_macros'], 
    'lampi':  [('SIMUPOP_MODULE', 'simuPOP_lampi'), ('LONGALLELE', None), ('SIMUMPI', None)] + MPIFlags['def_macros'],
    'laopmpi':[('SIMUPOP_MODULE', 'simuPOP_laopmpi'), ('LONGALLELE', None), ('OPTIMIZED', None), ('SIMUMPI', None)] + MPIFlags['def_macros'],
    'bampi':  [('SIMUPOP_MODULE', 'simuPOP_bampi'), ('BINARYALLELE', None), ('SIMUMPI', None)] + MPIFlags['def_macros'],
    'baopmpi':[('SIMUPOP_MODULE', 'simuPOP_baopmpi'), ('BINARYALLELE', None), ('OPTIMIZED', None), ('SIMUMPI', None)] + MPIFlags['def_macros'],
}
 
WRAP_INFO = {
    'std':    ['src/simuPOP_std_wrap.cpp', 'src/simuPOP_std.i', ''],
    'op':     ['src/simuPOP_op_wrap.cpp', 'src/simuPOP_op.i', '-DOPTIMIZED'],
    'la':     ['src/simuPOP_la_wrap.cpp', 'src/simuPOP_la.i', '-DLONGALLELE'],
    'laop':   ['src/simuPOP_laop_wrap.cpp', 'src/simuPOP_laop.i', '-DLONGALLELE -DOPTIMIZED'],
    'ba':     ['src/simuPOP_ba_wrap.cpp', 'src/simuPOP_ba.i', '-DBINARYALLELE'],
    'baop':   ['src/simuPOP_baop_wrap.cpp', 'src/simuPOP_baop.i', '-DBINARYALLELE -DOPTIMIZED'],
    'mpi':    ['src/simuPOP_mpi_wrap.cpp', 'src/simuPOP_mpi.i', '-DSIMUMPI'],
    'opmpi':  ['src/simuPOP_opmpi_wrap.cpp', 'src/simuPOP_opmpi.i', '-DOPTIMIZED -DSIMUMPI'],
    'lampi':  ['src/simuPOP_lampi_wrap.cpp', 'src/simuPOP_lampi.i', '-DLONGALLELE -DSIMUMPI'],
    'laopmpi':['src/simuPOP_laopmpi_wrap.cpp', 'src/simuPOP_laopmpi.i', '-DLONGALLELE -DOPTIMIZED -DSIMUMPI'],
    'bampi':  ['src/simuPOP_bampi_wrap.cpp', 'src/simuPOP_bampi.i', '-DBINARYALLELE -DSIMUMPI'],
    'baopmpi':['src/simuPOP_baopmpi_wrap.cpp', 'src/simuPOP_baopmpi.i', '-DBINARYALLELE -DOPTIMIZED -DSIMUMPI'],
}

for modu in MODULES:
    # source files
    MODU_INFO[modu] = {}
    MODU_INFO[modu]['src'] = []
    for src in SOURCE_FILES:
        mod_src = src[:-4] + '_' + modu + '.cpp'
        shutil.copy(src, mod_src)
        MODU_INFO[modu]['src'].append(mod_src)
    MODU_INFO[modu]['src'].append('src/simuPOP_' + modu + '_wrap.cpp' )
    MODU_INFO[modu]['src'].extend(GSL_FILES)
    # lib
    if os.name == 'nt':    # Windows
        MODU_INFO[modu]['libraries'] = ['libboost_serialization-mgw-mt-s-1_33_1', 'libboost_iostreams-mgw-mt-s-1_33_1', 'stdc++']
        MODU_INFO[modu]['libraries'].append('zdll')
    else:
        MODU_INFO[modu]['libraries'] = ['boost_serialization-%s' % TOOLSET, 'boost_iostreams-%s' % TOOLSET, 'stdc++']
        MODU_INFO[modu]['libraries'].append('z')
    MODU_INFO[modu]['include_dirs'] = ['.']
    #
    MODU_INFO[modu]['library_dirs'] = ['build']
    if os.name == 'nt':
        MODU_INFO[modu]['library_dirs'].append('win32')
    MODU_INFO[modu]['extra_compile_args'] = ['-O3']
    # define_macros
    MODU_INFO[modu]['define_macros'] = MACROS[modu]
    MODU_INFO[modu]['define_macros'].extend([('SIMUPOP_VER', SIMUPOP_VER), ('SIMUPOP_REV', SIMUPOP_REV)])
    MODU_INFO[modu]['undef_macros'] = []
    if 'mpi' in modu:
        MODU_INFO[modu]['include_dirs'].extend(MPIFlags['inc_dirs'])
        MODU_INFO[modu]['library_dirs'].extend(MPIFlags['lib_dirs'])
        MODU_INFO[modu]['libraries'].extend(MPIFlags['libs'])
        MODU_INFO[modu]['define_macros'].extend(MPIFlags['def_macros'])
        MODU_INFO[modu]['undef_macros'].extend(MPIFlags['undef_macros'])


if sys.argv[1] not in ['sdist']:
    # checking os type and copy configuration files
    if os.name == 'nt':    # Windows
        shutil.copy('config_win32.h', 'config.h')
    elif os.name == 'posix':
        sysName = os.uname()[0]
        if sysName == 'Linux':     # Linux
            shutil.copy('config_linux.h', 'config.h')
        elif sysName == 'SunOS': # Solaris
            shutil.copy('config_solaris.h', 'config.h')
        elif sysName == 'Darwin':    # MacOS
            shutil.copy('config_mac.h', 'config.h')
    else:
        try:
            open('config.h')
            close('config.h')
            print "Warning: Unknown system type. Using existing config.h"
        except IOError:
            print "Error: Unknown system type. Use configure to generate config.h."
            sys.exit()


        

############################################################################
#
# Generate Wrapping files
#
############################################################################

# if any of the wrap files does not exist
# or if the wrap files are older than any of the source files.
if (False in [os.path.isfile(WRAP_INFO[x][0]) for x in MODULES]) or \
    (max( [os.path.getmtime(x) for x in HEADER_FILES] ) > \
     min( [os.path.getmtime(WRAP_INFO[x][0]) for x in MODULES])):
    (v1, v2, v3) = swig_version()
    if (v1, v2, v3) >= (1, 3, 28):
        # for swig >= 1.3.28
        SWIG = 'swig -O -templatereduce -shadow -python -c++ -keyword -nodefaultctor -w-503,-312,-511,-362,-383,-384,-389,-315,-525'
    elif (v1, v2, v3) >= (1, 3, 25):
        # for swig from 1.3.25 to 1.3.27
        SWIG = 'swig -shadow -c++ -python -keyword -w-312,-401,-503,-511,-362,-383,-384,-389,-315,-525'
    else:
        print 'Swig >= 1.3.25 is required, please upgrade it.'
        sys.exit(1)
    # generate header file 
    print "Generating external runtime header file..."
    os.system( 'swig -python -external-runtime swigpyrun.h' )
    # try the first option set with the first library
    for lib in MODULES:
        print "Generating wrap file " + WRAP_INFO[lib][0]
        if os.system('%s %s -o %s %s' % (SWIG, WRAP_INFO[lib][2], WRAP_INFO[lib][0], WRAP_INFO[lib][1])) != 0:
            print "Calling swig failed. Please check your swig version."
            sys.exit(1)
    print
    print "All wrap files are generated successfully."
    print



############################################################################
#
# Build extensions
#
############################################################################


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
evolutionary scenarios. simuPOP is currently bundled with a Python
binding of coaSim.
"""
            

# find all test files
DATA_FILES = [
    ('share/simuPOP', ['README', 'INSTALL', 'ChangeLog', 'AUTHORS', 
        'COPYING', 'TODO', 'simuPOP.release']), 
    ('share/simuPOP/doc', ['doc/userGuide.pdf', 'doc/userGuide.py', 'doc/refManual.pdf']), 
    ('share/simuPOP/test', glob.glob('test/test_*.py')),
    ('share/simuPOP/misc', ['misc/README', 'misc/python-mode.el', 'misc/emacs-python.el']),
    ('share/simuPOP/scripts', glob.glob('scripts/*.py'))
]

EXT_MODULES = []
for modu in MODULES:
    EXT_MODULES.append(
        Extension('_simuPOP_%s' % modu,
            sources = MODU_INFO[modu]['src'],
            extra_compile_args = MODU_INFO[modu]['extra_compile_args'],
            include_dirs = MODU_INFO[modu]['include_dirs'],
            library_dirs = MODU_INFO[modu]['library_dirs'],
            libraries = MODU_INFO[modu]['libraries'],
            define_macros = MODU_INFO[modu]['define_macros'],
            undef_macros = MODU_INFO[modu]['undef_macros'],
        )
    )

setup(
    name = "simuPOP",
    version = SIMUPOP_VER,
    author = "Bo Peng",
    author_email = "bpeng@rice.edu",
    description = "Forward-time population genetics simulation environment",
    long_description = DESCRIPTION, 
    url = "http://simupop.sourceforge.net",
    package_dir = {'': 'src' }, 
    py_modules = SIMUPOP_FILES,
    ext_modules = EXT_MODULES,
    data_files = DATA_FILES
)


