"""
CoaSim/Python installer 
"""

from distutils.core import setup, Extension
import os, sys 

HEADER_FILES = [
  'all_markers.hh',
  'configuration.hh',
  'marker.hh',
  'simulator.hh',
  'builder_events.hh',
  'descender.hh',
  'micro_satellite_marker.hh',
  'snp_marker.hh',
  'builder.hh',
  'dist_funcs.hh',
  'monitor.hh',
  'compile_options.hh',
  'epochs.hh',
  'node.hh',
  'interval.hh',
  'retired_interval.hh',
  'trait_marker.hh'
]

SOURCE_FILES = [
  'configuration.cc',
  'marker.cc',
  'simulator.cc',
  'builder_events.cc',
  'descender.cc',
  'micro_satellite_marker.cc',
  'snp_marker.cc',
  'builder.cc',
  'dist_funcs.cc',
  'epochs.cc',
  'node.cc',
  'interval.cc',
  'retired_interval.cc',
  'trait_marker.cc'
]

WRAP_FILE = 'coaSim_wrap.cpp'
INTERFACE_FILE = 'coaSim.i'

# if the wrap files does not exist
# or if the wrap files are older than any of the source files.
if not os.path.isfile(WRAP_FILE) or \
  max( [os.path.getmtime(x) for x in HEADER_FILES + SOURCE_FILES + [INTERFACE_FILE]] ) > \
   os.path.getmtime(WRAP_FILE):
  try:
    # for swig >= 1.3.28
    SWIG1 = 'swig -O -templatereduce -shadow -python -c++ -keyword -nodefaultctor -w-503,-312,-511,-362,-383,-384,-389,-315,-525'
    # for swig <= 1.3.27, >= 1.3.25
    SWIG2 = 'swig -shadow -c++ -python -keyword -nodefault -w-312,-401,-503,-511,-362,-383,-384,-389,-315,-525'
    # generate header file 
    # try the first option set with the first library
    print "Generating wrap file coasim_wrap.cpp"
    if os.system('%s -o %s %s' % (SWIG1, WRAP_FILE, INTERFACE_FILE)) != 0:
      print "Your swig version is not up to date. Trying options with swig <= 1.3.27"
      if os.system('%s -o %s %s' % (SWIG2, WRAP_FILE, INTERFACE_FILE)) != 0:
        print "None of the swig option sets works, please check if you have SWIG >= 1.3.25 installed"
        sys.exit(1)
    print
    print WRAP_FILE, "is generated successfully."
    print
  except:
    print "Can not generate wrap files. Please check your swig installation."
    raise

DESCRIPTION = """
CoaSim is a tool for simulating the coalescent process with recombination
and geneconversion under various demographic models. It effectively 
constructs the ancestral recombination graph for a given number of 
individuals and uses this to simulate samples of SNP, micro-satellite, 
and other haplotypes/genotypes. The generated sample can afterwards be 
separated in cases and controls, depending on states of selected individual
markers. The tool can accordingly also be used to construct cases and 
control data sets for association studies.
"""

setup(
  name = "CoaSim",
  version = "4.0.5",
  author = "T Mailund, wrapped by Bo Peng",
  author_email = "mailund@birc.au.dk, bpeng@rice.edu",
  description = "Coalescent-based population genetics simulation program",
  long_description = DESCRIPTION, 
  url = "http://www.birc.dk/Software/CoaSim/",
  py_modules = ['coaSim'],
  ext_modules = [
    Extension('_coaSim',
      include_dirs = ["."],
      libraries = ['stdc++'],
      sources = SOURCE_FILES + [WRAP_FILE]
    ),
  ],
)
