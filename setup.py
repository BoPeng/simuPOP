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

SIMUPOP_VER="snapshot"
if sys.argv[1] not in ['sdist']:
  shutil.copy('src/utility.cpp', 'src/utility_std.cpp')
  shutil.copy('src/utility.cpp', 'src/utility_op.cpp')
  shutil.copy('src/utility.cpp', 'src/utility_la.cpp')
  shutil.copy('src/utility.cpp', 'src/utility_laop.cpp')
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
    
def addCarrayEntry(file):
  ' add a line at the wrap file for carray type definition'
  shutil.copy(file, file+'tmp')
  ofile = open(file, 'w')
  ifile = open(file+'tmp', 'r')
  for line in ifile.readlines():
    if line.find('static PyMethodDef SwigMethods[] = ') != -1:
      ofile.write('''static PyMethodDef SwigMethods[] = {
        /* add carray entries */ 
        { (char*)"carray", a_array, METH_VARARGS, a_array_doc}, 
        ''')
    else:
      ofile.write(line)
  ifile.close()
  ofile.close()
  os.remove(file+'tmp')

# run swig, modify generated wrap file.
if not ( os.path.isfile('src/simuPOP_std_wrap.cpp') and \
  os.path.isfile('src/simuPOP_op_wrap.cpp') and \
  os.path.isfile('src/simuPOP_la_wrap.cpp') and \
  os.path.isfile('src/simuPOP_laop_wrap.cpp') ):
  SWIG = 'swig  -shadow -python -keyword -w-503,-312,-511,-362,-383,-384,-389,-315,-525 -nodefault -c++ '
  print "Generating wrap file for standard library..."
  os.system(SWIG + ' -o src/simuPOP_std_wrap.cpp src/simuPOP_std.i')
  addCarrayEntry('src/simuPOP_std_wrap.cpp')
  # for optimized library
  print "Generating wrap file for optimzied library..."
  os.system(SWIG + ' -DOPTIMIZED -o src/simuPOP_op_wrap.cpp src/simuPOP_op.i')
  addCarrayEntry('src/simuPOP_op_wrap.cpp')
  # for long allele library
  print "Generating wrap file for long allele library..."
  os.system(SWIG + ' -DLONGALLELE -o src/simuPOP_la_wrap.cpp src/simuPOP_la.i')
  addCarrayEntry('src/simuPOP_la_wrap.cpp')
  # for long allele optimized library
  print "Generating wrap file for long allele library..."
  os.system(SWIG + ' -DLONGALLELE -DOPTIMIZED -o src/simuPOP_laop_wrap.cpp src/simuPOP_laop.i')
  addCarrayEntry('src/simuPOP_laop_wrap.cpp')

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
    'COPYING', 'TODO']), 
  ('share/simuPOP/doc', ['doc/userGuide.pdf', 'doc/userGuide.py', 'doc/refManual.pdf']), 
  ('share/simuPOP/test', glob.glob('test/test_*.py')),
  ('share/simuPOP/misc', ['misc/README', 'misc/python-mode.el', 'misc/emacs-python.el']),
  ('share/simuPOP/scripts', ['scripts/simuComplexDisease.py',  
   'scripts/simuLDDecay.py', 'scripts/simuCDCV.py', 'scripts/simuRecHotSpot.py']),
]

if XML_SUPPORT: # using serialization library with xml support
  serial_macro = []
else:
  serial_macro = [ ('__NO_XML_SUPPORT__', None) ]
  
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
      extra_compile_args=['-O2'],
      include_dirs = ["."],
      libraries = ['stdc++'],
      define_macros = serial_macro,
      sources= GSL_FILES + SERIAL_FILES + [
        'src/simuPOP_std_wrap.cpp',
        'src/utility_std.cpp'] 
    ),
    Extension('_simuPOP_op',
      extra_compile_args=['-O2'],
      include_dirs = ["."],
      libraries = ['stdc++'],
      define_macros = [ ('OPTIMIZED', None)] + serial_macro,
      sources= GSL_FILES + SERIAL_FILES + [
        'src/simuPOP_op_wrap.cpp',
        'src/utility_op.cpp'] 
    ),
    Extension('_simuPOP_la',
      extra_compile_args=['-O2'],
      include_dirs = ["."],
      libraries = ['stdc++'],
      define_macros = [ ('LONGALLELE', None) ] + serial_macro,
      sources= GSL_FILES + SERIAL_FILES + [
        'src/simuPOP_la_wrap.cpp',
        'src/utility_la.cpp'] 
    ),
    Extension('_simuPOP_laop',
      extra_compile_args=['-O2'],
      include_dirs = ["."],
      libraries = ['stdc++'],
      define_macros = [ ('LONGALLELE', None), ('OPTIMIZED', None) ] + serial_macro,
      sources= GSL_FILES + SERIAL_FILES + [
        'src/simuPOP_laop_wrap.cpp',
        'src/utility_laop.cpp'] 
    )
  ],
  data_files = DATA_FILES
 )
