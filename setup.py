"""
simuPOP installer 

In case that you have modified the C++ code, or you are
checking out simuPOP from svn repository, you need to 
install swig >= 1.3.27 to generate the wrap files.

"""
import os, sys, shutil, glob, re, tempfile

# simuPOP works with these boost versions.
boost_versions = ['1_33_1', '1_34_0', '1_34_1', '1_35_0', '1_36_0']

# if the source package is districuted with boost 1.33.1, use
# this bundled version.
included_version = [x for x in boost_versions if os.path.isdir('boost_' + x)]
included_boost_dir = 'boost_1_36_0'
if len(included_version) > 0:
    included_boost = True
    included_boost_dir = 'boost_' + included_version[0]
    included_boost_include_dir = included_boost_dir
    included_boost_serialization_dir = os.path.join(included_boost_dir, 'libs', 'serialization', 'src')
    included_boost_iostreams_dir = os.path.join(included_boost_dir, 'libs', 'iostreams', 'src')
else:
    included_boost = False

# If setup.py can not find boost libraries, change boost_lib_seach_paths
# and/or boost_inc_search_paths. 
#
# use_vc is used only once to indicate if a portable stdint.h
# need to be used. (msvc does not ship with stdint.h)
if os.name == 'nt':
    use_vc = True
    # under windows, boost/iostreams/gzip decompressor seems
    # to be broken. has to be disabled by now.
    disable_compression = True
    # win32 is default since 1.33.1 libraries are bundled with simuPOP windows
    # distribution
    boost_lib_search_paths = [r'win32', r'c:\boost\lib', r'c:\program files\boost\lib']
    boost_inc_search_paths = [included_boost_dir, r'c:\boost', r'c:\program files\boost']
    boost_lib_prefix = ''
    boost_lib_suffix = '.lib'
else:    
    use_vc = False
    disable_compression = False
    boost_lib_search_paths = ['/usr/lib', '/usr/local/lib']
    boost_inc_search_paths = [included_boost_dir, '/usr/include', '/usr/local/include']
    home = os.environ.get('HOME', None)
    if home is not None:
        dirs = glob.glob(os.path.join(home, 'boost*'))
        for dir in dirs:
            if os.path.isdir(os.path.join(dir, 'lib')):
                boost_lib_search_paths.append(os.path.join(dir, 'lib'))
            if os.path.isdir(os.path.join(dir, 'include')):
                boost_inc_search_paths.append(os.path.join(dir, 'include'))
            else:
                boost_inc_search_paths.append(dir)
    boost_lib_prefix = 'lib'
    boost_lib_suffix = '.a'

# if you need to use full path name for swig, change it here.
SWIG = 'swig'

############################################################################
#
# THE FOLLOWING IS NOT SUPPOSED TO BE MODIFIED
#
############################################################################

from distutils.core import setup, Extension


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


def getBoostLibraries(libs, lib_paths, lib_prefix, lib_suffix, inc_paths, versions):
    ''' look for boost libraries
      libs: library names
      lib_paths: try these paths for boost libraries
      inc_paths: try these paths for boost headers
      versions:   supported boost versions
    '''
    found_lib = False
    found_inc = False
    lib_names = []
    lib_path = None
    inc_path = None
    for path in lib_paths:
        if path is None:
            continue
        for lib in libs:
            # get all the libs, then filter for the right library
            files = glob.glob(os.path.join(path, '%sboost_%s*.*' % (lib_prefix, lib)))
            # check things like libboost_iostreams-gcc-mt-d-1_33_1.a
            if len(files) > 0:
                # runtime code includes s,g,y,d,p,n, where we should look for
                # d,g,y for debug, s,p,n for release
                lib_files = []
                for ver in versions:
                    lib_files += filter(lambda x: re.search('%sboost_%s-\w+-mt-([^dgy]+-)*%s%s' \
                        % (lib_prefix, lib, ver, lib_suffix), x), files)
                if len(lib_files) == 0:
                    # use alternative libraries
                    for ver in versions:
                        lib_files += filter(lambda x: re.search('%sboost_%s[-\w]*[-%s]*%s' \
                            % (lib_prefix, lib, ver, lib_suffix), x), files)
                if len(lib_files) > 0:
                    # get xxx-gcc-1_33_1 from /usr/local/lib/libboost_xxx-gcc-1_33_1.a
                    name = lib_files[0].split(os.sep)[-1][len(lib_prefix):]
                    lib_names.append(name.split('.')[0])
                else:
                    break
        if len(lib_names) == len(libs):
            found_lib = True
            lib_path = path
            break
    if not found_lib:
        print "Can not find boost libraries. Please read the front part"
        print "of setup.py for instructions."
        print "Library search paths: ", lib_paths
        print "Include search paths: ", inc_paths
        sys.exit(1)
    # check version number in boost/version.hpp
    def isValidBoostDir(dir):
        version_file = os.path.join(dir, 'boost', 'version.hpp')
        if not os.path.isfile(version_file):
            return False
        version_file_content = open(version_file).read()
        version_strings = ['#define BOOST_LIB_VERSION "%s"' % ver for ver in versions]
        return True in [x in version_file_content for x in version_strings]
    # check for boost header file
    for path in inc_paths:
        if path is None:
            continue
        if isValidBoostDir(path):
            inc_path = path
            found_inc = True
        else:   # check path/boost_1_xx_x/boost
            dirs = glob.glob(os.path.join(path, 'boost-*')) + glob.glob(os.path.join(path, 'boost_*'))
            if len(dirs) > 0 and isValidBoostDir(dirs[0]):
                inc_path = dirs[0]
                found_inc = True
    # return result
    if found_inc:
        return (lib_names, lib_path, inc_path)
    else:
        print "Can not find boost libraries. Please read the front part"
        print "of setup.py for instructions."
        sys.exit(1)
        
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
    'src/sampler.h',
    'src/qtrait.h',
    'src/stator.h',
    'src/terminator.h',
    'src/mutator.h',
    'src/recombinator.h',
    'src/tagger.h',
    'src/misc.h',
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
    'src/sampler.cpp',
    'src/qtrait.cpp',
    'src/stator.cpp',
    'src/terminator.cpp',
    'src/mutator.cpp',
    'src/recombinator.cpp',
    'src/tagger.cpp',
    'src/misc.cpp',
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
]

if included_boost and os.name != 'nt':
    LIB_FILES.extend([os.path.join(included_boost_serialization_dir, x) for x in [
        'basic_archive.cpp',
        'basic_iarchive.cpp',
        'basic_oarchive.cpp',
        'basic_serializer_map.cpp',
        'basic_text_iprimitive.cpp',
        'basic_text_oprimitive.cpp',
        'binary_iarchive.cpp',
        'binary_oarchive.cpp',
        'extended_type_info.cpp',
        'extended_type_info_no_rtti.cpp',
        'extended_type_info_typeid.cpp',
        'text_iarchive.cpp',
        'text_oarchive.cpp',
        'void_cast.cpp',
        'polymorphic_iarchive.cpp',
        'polymorphic_oarchive.cpp',
        'stl_port.cpp',
        'basic_pointer_iserializer.cpp',
        'basic_iserializer.cpp',
        'basic_oserializer.cpp',
        'basic_pointer_oserializer.cpp',
        'basic_xml_archive.cpp',
        'xml_grammar.cpp',
        'xml_iarchive.cpp',
        'xml_oarchive.cpp'
        ]
    ])


if included_boost and os.name != 'nt':
    LIB_FILES.extend([os.path.join(included_boost_iostreams_dir, x) for x in [
        'mapped_file.cpp',
        'file_descriptor.cpp',
        'zlib.cpp'
        ]
    ])


SIMUPOP_FILES = [
    'simuPOP', 
    'simuOpt', 
    'simuUtil', 
    'hapMapUtil',
    'simuRPy', 
]


#
# DETECT BASIC SYSTEM SETTINGS
#


# I use -O option for 'Optimization' and part of my extension code actually
# depend on this feature (XX_swiginit calls). However, this makes __doc__
# messages of member functions disappear, as discussed in SWIG user mailing list.
# 
# I am using a patched SVN swig version to fix this (as of Apr. 2007). Users 
# of later version of swig (>1.3.21) may not need this. To test if __doc__ works, 
# try
#     help(population.absIndIndex)
#
# The patch can be found here:
#     http://sf.net/tracker/index.php?func=detail&aid=1700146&group_id=1645&atid=301645
SWIG_FLAGS = '-O -templatereduce -shadow -python -c++ -keyword -nodefaultctor -w-503,-312,-511,-362,-383,-384,-389,-315,-509,-525'
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
evolutionary scenarios. simuPOP is currently bundled with a Python
binding of coaSim.
"""
            

# find all test files
DATA_FILES = [
    ('share/simuPOP', ['README', 'INSTALL', 'ChangeLog', 'AUTHORS', 
        'COPYING', 'TODO', 'simuPOP.release']), 
    ('share/simuPOP/doc', ['doc/userGuide.pdf', 'doc/userGuide.py', 'doc/refManual.pdf']), 
    ('share/simuPOP/doc/cookbook', ['doc/cookbook/%s' % x for x in 
        ('Mating_assortativeMating.py', 'Mating_overlappingGeneration.py',
         'Mating_pyParentsChooser.py', 'Mating_selfing.py', 'Mating_consanguineous.py',
         'Mating_pyMating_cpp.py', 'Mating_pyMating_cpp.i', 'Mating_pyMating_cpp.h',
         'Operator_pyOperator.py')]),
    ('share/simuPOP/test', glob.glob('test/test_*.py') + \
        ['test/run_tests.py', 'test/run_tests.sh']),
    ('share/simuPOP/misc', ['misc/README', 'misc/python-mode.el', 'misc/emacs-python.el']),
    ('share/simuPOP/scripts', glob.glob('scripts/*.py'))
]

if os.name == 'nt':
    DATA_FILES += [('Lib/site-packages', ['win32/zlib1.dll'])]


def ModuInfo(modu, SIMUPOP_VER='9.9.9', SIMUPOP_REV='9999'):
    if included_boost and os.name != 'nt':
        boost_inc_path = included_boost_include_dir
        boost_lib_names = []
        boost_lib_path = None
    else:
        (boost_lib_names, boost_lib_path, boost_inc_path) = getBoostLibraries(
            ['iostreams', 'serialization'], boost_lib_search_paths,
            boost_lib_prefix, boost_lib_suffix,
            boost_inc_search_paths, boost_versions)
    res = {}
    res['src'] =  ['src/simuPOP_' + modu + '_wrap.cpp']
    for src in SOURCE_FILES:
        res['src'].append(src[:-4] + '_' + modu + '.cpp')
    res['src'].extend(LIB_FILES)
    # lib
    if os.name == 'nt':    # Windows
        res['libraries'] = ['zdll']
    else:
        res['libraries'] = ['stdc++', 'z']
    res['libraries'].extend(boost_lib_names)
    res['include_dirs'] = ['.', boost_inc_path]
    if use_vc:
        # I have a portable stdint.h for msvc
        res['include_dirs'].append('win32')
    #
    if included_boost:
        res['library_dirs'] = ['build']
    else:
        res['library_dirs'] = ['build', boost_lib_path]
    if use_vc:
        res['library_dirs'].append('win32')
    if os.name == 'nt':
        # msvc does not have O3 option
        res['extra_compile_args'] = ['/O2']
    else:
        # force the use of static boost libraries because I do not
        # want to bundle boost libraries with simuPOP distributions.
        res['extra_compile_args'] = ['-O3', '-Wall']
    # define_macros (deep copy)
    res['define_macros'] = [x for x in MACROS[modu]]
    res['define_macros'].extend([('SIMUPOP_VER', SIMUPOP_VER), ('SIMUPOP_REV', SIMUPOP_REV)])
    if disable_compression:
        res['define_macros'].extend([('DISABLE_COMPRESSION', None)])
    if os.name == 'nt':
        res['define_macros'].extend([('BOOST_ALL_NO_LIB', None)])
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
    # for every official release, there will be a file recording release info
    # Othersise, SIMUPOP_VER and SIMUPOP_REV will be provided as environmental
    # variables.
    if os.environ.has_key('SIMUPOP_VER') and os.environ.has_key('SIMUPOP_REV'):
        SIMUPOP_VER = os.environ['SIMUPOP_VER']
        SIMUPOP_REV = os.environ['SIMUPOP_REV']
    else:
        execfile('simuPOP.release')
    # create source file for each module
    MODULES = ['std', 'op', 'la', 'laop', 'ba', 'baop']
    SIMUPOP_FILES += ['simuPOP_%s' % x for x in MODULES]
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
            print "Generating wrap file " + WRAP_INFO[lib][0]
            if os.system('%s %s -outdir %s %s -o %s %s' % (SWIG, SWIG_FLAGS, \
                SWIG_OUTDIR, WRAP_INFO[lib][2], WRAP_INFO[lib][0], WRAP_INFO[lib][1])) != 0:
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
            Extension('_simuPOP_%s' % modu,
                sources = info['src'],
                extra_compile_args = info['extra_compile_args'],
                include_dirs = info['include_dirs'],
                library_dirs = info['library_dirs'],
                libraries = info['libraries'],
                define_macros = info['define_macros'],
                undef_macros = info['undef_macros'],
            )
        )
    #
    setup(
        name = "simuPOP",
        version = SIMUPOP_VER,
        author = "Bo Peng",
        author_email = "bpeng@mdanderson.org",
        description = "Forward-time population genetics simulation environment",
        long_description = DESCRIPTION, 
        url = "http://simupop.sourceforge.net",
        package_dir = {'': 'src'}, 
        py_modules = SIMUPOP_FILES,
        ext_modules = EXT_MODULES,
        data_files = DATA_FILES,
    )
    # remove copied files
    for file in copied_files:
        os.remove(file)


