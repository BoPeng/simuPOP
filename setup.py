"""
simuPOP installer 

In case that you have modified the C++ code, or you are
checking out simuPOP from svn repository, you need to 
install swig >= 1.3.23 (and remove src/*wrap.cpp if they
exist) to generate these wrap files.
"""

#
# XML_SUPPORT will be disabled for mac machines due to a bug
# in mac/gcc. You will not be able to save population in xml 
# format under a mac. Binary and txt formats are still supported
# and should suffice most applications.
#
# If you do need XML_SUPPORT under mac, you can go to line 42
# set XML_SUPPORT=True and try to compile.
#
XML_SUPPORT = True

#########################################################

from distutils.core import setup, Extension

import os, shutil, sys, glob, re

# handle version and revision
# 1. use environmental variable if possible
#
# for every official release, there will be a file recording release info
execfile('simuPOP.release')
std_macro = [('SIMUPOP_VER', SIMUPOP_VER), 
             ('SIMUPOP_REV', SIMUPOP_REV) ]
  
if sys.argv[1] not in ['sdist']:
  shutil.copy('src/utility.cpp', 'src/utility_std.cpp')
  shutil.copy('src/utility.cpp', 'src/utility_op.cpp')
  shutil.copy('src/utility.cpp', 'src/utility_la.cpp')
  shutil.copy('src/utility.cpp', 'src/utility_laop.cpp')
  shutil.copy('src/utility.cpp', 'src/utility_ba.cpp')
  shutil.copy('src/utility.cpp', 'src/utility_baop.cpp')
  # checking os type and copy configuration files
  if os.name == 'nt':  # Windows
    shutil.copy('config_win32.h', 'config.h')
  elif os.name == 'posix':
    sysName = os.uname()[0]
    if sysName == 'Linux':   # Linux
      shutil.copy('config_linux.h', 'config.h')
    elif sysName == 'SunOS': # Solaris
      shutil.copy('config_solaris.h', 'config.h')
    elif sysName == 'Darwin':  # MacOS
      shutil.copy('config_mac.h', 'config.h')
      XML_SUPPORT = False
  else:
    try:
      open('config.h')
      close('config.h')
      print "Warning: Unknown system type. Using existing config.h"
    except IOError:
      print "Error: Unknown system type. Use configure to generate config.h."
      sys.exit()
    
# run swig, modify generated wrap file.
# determine if the source if newer than _wrap files.
SOURCE_FILES = [
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
  'src/utility.cpp'
]

WRAP_INFO = [
  ['src/simuPOP_std_wrap.cpp', 'src/simuPOP_std.i', ''],
  ['src/simuPOP_op_wrap.cpp', 'src/simuPOP_op.i', '-DOPTIMIZED'],
  ['src/simuPOP_la_wrap.cpp', 'src/simuPOP_la.i', '-DLONGALLELE'],
  ['src/simuPOP_laop_wrap.cpp', 'src/simuPOP_laop.i', '-DLONGALLELE -DOPTIMIZED'],
  ['src/simuPOP_ba_wrap.cpp', 'src/simuPOP_ba.i', '-DBINARYALLELE'],
  ['src/simuPOP_baop_wrap.cpp', 'src/simuPOP_baop.i', '-DBINARYALLELE -DOPTIMIZED']
]  

# if any of the wrap files does not exist
# or if the wrap files are older than any of the source files.
if (False in [os.path.isfile(WRAP_INFO[x][0]) for x in range(len(WRAP_INFO))]) or \
  (max( [os.path.getmtime(x) for x in SOURCE_FILES] ) > \
   min( [os.path.getmtime(WRAP_INFO[x][0]) for x in range(len(WRAP_INFO))])):
  try:
    # for swig >= 1.3.28
    SWIG1 = 'swig -O -templatereduce -shadow -python -c++ -keyword -nodefaultctor -w-503,-312,-511,-362,-383,-384,-389,-315,-525'
    # for swig <= 1.3.27, >= 1.3.25
    SWIG2 = 'swig -shadow -c++ -python -keyword -w-312,-401,-503,-511,-362,-383,-384,-389,-315,-525'
    # generate header file 
    print "Generating external runtime header file..."
    os.system( 'swig  -python -external-runtime swigpyrun.h' )
    # try the first option set with the first library
    lib = 0
    print "Generating wrap file " + WRAP_INFO[lib][0]
    if os.system('%s %s -o %s %s' % (SWIG1, WRAP_INFO[lib][2], WRAP_INFO[lib][0], WRAP_INFO[lib][1])) != 0:
      print "Your swig version is not up to date. Trying options with swig <= 1.3.27"
      # all libraries
      for lib in range(len(WRAP_INFO)):
        print "Generating wrap file " + WRAP_INFO[lib][0]
        if os.system('%s %s -o %s %s' % (SWIG2, WRAP_INFO[lib][2], WRAP_INFO[lib][0], WRAP_INFO[lib][1])) != 0:
          print "None of the swig option sets works, please check if you have SWIG >= 1.3.25 installed"
          sys.exit(1)
    # if OK, generate the rest of them, and I do not expect any error
    else:
      for lib in range(1,len(WRAP_INFO)):
        print "Generating wrap file " + WRAP_INFO[lib][0]
        os.system('%s %s -o %s %s' % (SWIG2, WRAP_INFO[lib][2], WRAP_INFO[lib][0], WRAP_INFO[lib][1]))
  except:
    print "Can not generate wrap files. Please check your swig installation."
    raise

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
  'gsl/randist/chisq.c',
  'gsl/randist/gauss.c',
  'gsl/error.c' 
]

SERIAL_FILES = [
  'src/serialization/basic_archive.cpp',
  'src/serialization/basic_iarchive.cpp',
  'src/serialization/basic_oarchive.cpp',
  'src/serialization/basic_serializer_map.cpp',
  'src/serialization/basic_text_iprimitive.cpp',
  'src/serialization/basic_text_oprimitive.cpp',
  'src/serialization/binary_iarchive.cpp',
  'src/serialization/binary_oarchive.cpp',
  'src/serialization/extended_type_info.cpp',
  'src/serialization/extended_type_info_no_rtti.cpp',
  'src/serialization/extended_type_info_typeid.cpp',
  'src/serialization/text_iarchive.cpp',
  'src/serialization/text_oarchive.cpp',
  'src/serialization/void_cast.cpp',
  'src/serialization/polymorphic_iarchive.cpp',
  'src/serialization/polymorphic_oarchive.cpp'
]
if XML_SUPPORT: 
  SERIAL_FILES.extend( [
    'src/serialization/basic_xml_archive.cpp',
    # 'src/serialization/basic_xml_grammar.ipp',
    'src/serialization/xml_grammar.cpp',
    'src/serialization/xml_iarchive.cpp',
    'src/serialization/xml_oarchive.cpp'
    ]
  )

# find all test files
DATA_FILES =  [
  ('share/simuPOP', ['README', 'INSTALL', 'ChangeLog', 'AUTHORS', 
    'COPYING', 'TODO', 'simuPOP.release']), 
  ('share/simuPOP/doc', ['doc/userGuide.pdf', 'doc/userGuide.py', 'doc/refManual.pdf']), 
  ('share/simuPOP/test', glob.glob('test/test_*.py')),
  ('share/simuPOP/misc', ['misc/README', 'misc/python-mode.el', 'misc/emacs-python.el']),
  ('share/simuPOP/scripts', glob.glob('scripts/*.py')),
]

if not XML_SUPPORT: # using serialization library with xml support
  std_macro.append( ('__NO_XML_SUPPORT__', None) )
  
setup(
  name = "simuPOP",
  version = SIMUPOP_VER,
  author = "Bo Peng",
  author_email = "bpeng@rice.edu",
  description = "Forward-based population genetics simulation framework",
  long_description = DESCRIPTION, 
  url = "http://bp6.stat.rice.edu:8080/simuPOP",
  package_dir = {'': 'src'}, 
  py_modules = ['simuPOP', 'simuOpt', 'simuPOP_std', 'simuPOP_op', 'simuPOP_la', 'simuPOP_laop', 
    'simuUtil', 'simuSciPy', 'simuMatPlt', 'simuRPy', 'simuViewPop'],
  ext_modules = [
    Extension('_simuPOP_std',
      extra_compile_args=['-O3'],
      include_dirs = ["."],
      libraries = ['stdc++'],
      define_macros = [ ('SIMUPOP_MODULE', 'simuPOP_std')] + std_macro,
      sources = GSL_FILES + SERIAL_FILES + [
        'src/simuPOP_std_wrap.cpp',
        'src/utility_std.cpp'] 
    ),
    Extension('_simuPOP_op',
      extra_compile_args=['-O3'],
      include_dirs = ["."],
      libraries = ['stdc++'],
      define_macros = [ ('SIMUPOP_MODULE', 'simuPOP_op'), ('OPTIMIZED', None)] + std_macro,
      sources = GSL_FILES + SERIAL_FILES + [
        'src/simuPOP_op_wrap.cpp',
        'src/utility_op.cpp'] 
    ),
    Extension('_simuPOP_la',
      extra_compile_args=['-O3'],
      include_dirs = ["."],
      libraries = ['stdc++'],
      define_macros = [ ('SIMUPOP_MODULE', 'simuPOP_la'), ('LONGALLELE', None) ] + std_macro,
      sources = GSL_FILES + SERIAL_FILES + [
        'src/simuPOP_la_wrap.cpp',
        'src/utility_la.cpp'] 
    ),
    Extension('_simuPOP_laop',
      extra_compile_args=['-O3'],
      include_dirs = ["."],
      libraries = ['stdc++'],
      define_macros = [ ('SIMUPOP_MODULE', 'simuPOP_laop'), ('LONGALLELE', None), 
                        ('OPTIMIZED', None) ] + std_macro,
      sources = GSL_FILES + SERIAL_FILES + [
        'src/simuPOP_laop_wrap.cpp',
        'src/utility_laop.cpp'] 
    ),
    Extension('_simuPOP_ba',
      extra_compile_args=['-O3'],
      include_dirs = ["."],
      libraries = ['stdc++'],
      define_macros = [ ('SIMUPOP_MODULE', 'simuPOP_ba'), ('BINARYALLELE', None) ] + std_macro,
      sources = GSL_FILES + SERIAL_FILES + [
        'src/simuPOP_ba_wrap.cpp',
        'src/utility_ba.cpp'] 
    ),
    Extension('_simuPOP_baop',
      extra_compile_args=['-O3'],
      include_dirs = ["."],
      libraries = ['stdc++'],
      define_macros = [ ('SIMUPOP_MODULE', 'simuPOP_baop'), ('BINARYALLELE', None), 
                        ('OPTIMIZED', None) ] + std_macro,
      sources = GSL_FILES + SERIAL_FILES + [
        'src/simuPOP_baop_wrap.cpp',
        'src/utility_baop.cpp'] 
    )
  ],
  data_files = DATA_FILES
)
