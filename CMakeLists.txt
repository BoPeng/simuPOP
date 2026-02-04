cmake_minimum_required(VERSION 3.20)
project(simuPOP LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Options
option(USE_OPENMP "Enable OpenMP support" ON)
option(BUILD_TESTING "Build tests" OFF)

# Find Python
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

# Find SWIG
find_package(SWIG 4.0 REQUIRED COMPONENTS python)
include(UseSWIG)

# Find OpenMP (optional)
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_FLAGS}")
    else()
        message(STATUS "OpenMP not found, building without OpenMP support")
        set(USE_OPENMP OFF)
    endif()
endif()

# Find zlib
find_package(ZLIB REQUIRED)

#------------------------------------------------------------------------------
# Platform and Compiler Feature Detection (replaces config.h)
#------------------------------------------------------------------------------
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckFunctionExists)
include(CheckSymbolExists)
include(CheckCXXSourceCompiles)

# Check for header files
check_include_file("stdint.h" HAVE_STDINT_H)
check_include_file("inttypes.h" HAVE_INTTYPES_H)
check_include_file("float.h" HAVE_FLOAT_H)
check_include_file("dlfcn.h" HAVE_DLFCN_H)
check_include_file("unistd.h" HAVE_UNISTD_H)
check_include_file("ieeefp.h" HAVE_IEEEFP_H)

# Check for math functions
set(CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(acosh HAVE_ACOSH)
check_function_exists(asinh HAVE_ASINH)
check_function_exists(atanh HAVE_ATANH)
check_function_exists(expm1 HAVE_EXPM1)
check_function_exists(log1p HAVE_LOG1P)
check_function_exists(hypot HAVE_HYPOT)

# Check for declarations in math.h
check_symbol_exists(isfinite "math.h" HAVE_DECL_ISFINITE)
check_symbol_exists(finite "math.h" HAVE_DECL_FINITE)
check_symbol_exists(isinf "math.h" HAVE_DECL_ISINF)
check_symbol_exists(isnan "math.h" HAVE_DECL_ISNAN)
check_symbol_exists(acosh "math.h" HAVE_DECL_ACOSH)
check_symbol_exists(asinh "math.h" HAVE_DECL_ASINH)
check_symbol_exists(atanh "math.h" HAVE_DECL_ATANH)
check_symbol_exists(expm1 "math.h" HAVE_DECL_EXPM1)
check_symbol_exists(log1p "math.h" HAVE_DECL_LOG1P)
check_symbol_exists(frexp "math.h" HAVE_DECL_FREXP)
check_symbol_exists(ldexp "math.h" HAVE_DECL_LDEXP)
check_symbol_exists(hypot "math.h" HAVE_DECL_HYPOT)

# Check for IEEE comparisons
check_cxx_source_compiles("
    #include <cmath>
    int main() {
        double x = 0.0/0.0;
        return (x != x) ? 0 : 1;
    }
" HAVE_IEEE_COMPARISONS)

# Check for unordered_map (TR1 support)
check_cxx_source_compiles("
    #include <unordered_map>
    int main() {
        std::unordered_map<int, int> m;
        return 0;
    }
" HAVE_UNORDERED_MAP)

check_cxx_source_compiles("
    #include <tr1/unordered_map>
    int main() {
        std::tr1::unordered_map<int, int> m;
        return 0;
    }
" HAVE_TR1_UNORDERED_MAP)

if(HAVE_UNORDERED_MAP)
    set(TR1_SUPPORT 1)
elseif(HAVE_TR1_UNORDERED_MAP)
    set(TR1_SUPPORT 2)
else()
    set(TR1_SUPPORT 0)
endif()

# Detect libc++ vs libstdc++
check_cxx_source_compiles("
    #include <vector>
    #include <iostream>
    int main() {
        std::cout << std::_S_word_bit;
        return 0;
    }
" HAVE_LIBSTDCPP)

if(NOT HAVE_LIBSTDCPP)
    set(USE_LIBCPP 1)
else()
    set(USE_LIBCPP 0)
endif()

# Generate config headers
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/SimuPOPConfig.h.in
    ${CMAKE_BINARY_DIR}/generated/simuPOP_config.h
    @ONLY
)

# Also generate config.h for GSL compatibility (GSL files use #include <config.h>)
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/SimuPOPConfig.h.in
    ${CMAKE_BINARY_DIR}/generated/config.h
    @ONLY
)

#------------------------------------------------------------------------------
# Boost handling
#------------------------------------------------------------------------------
# Look for existing boost directories
file(GLOB BOOST_DIRS "${CMAKE_SOURCE_DIR}/boost_*")
list(FILTER BOOST_DIRS EXCLUDE REGEX ".*\\.tar\\.gz$")
list(SORT BOOST_DIRS ORDER DESCENDING)

if(BOOST_DIRS)
    list(GET BOOST_DIRS 0 BOOST_DIR)
    message(STATUS "Using Boost from: ${BOOST_DIR}")
else()
    # Download boost if not found
    set(BOOST_VERSION "1_86_0")
    set(BOOST_DOT_VERSION "1.86.0")
    set(BOOST_DIR "${CMAKE_SOURCE_DIR}/boost_${BOOST_VERSION}")

    if(NOT EXISTS ${BOOST_DIR})
        message(STATUS "Downloading Boost ${BOOST_DOT_VERSION}...")
        set(BOOST_URL "https://boostorg.jfrog.io/artifactory/main/release/${BOOST_DOT_VERSION}/source/boost_${BOOST_VERSION}.tar.gz")
        file(DOWNLOAD ${BOOST_URL} "${CMAKE_SOURCE_DIR}/boost_${BOOST_VERSION}.tar.gz"
             SHOW_PROGRESS
             STATUS DOWNLOAD_STATUS)
        list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
        if(NOT STATUS_CODE EQUAL 0)
            message(FATAL_ERROR "Failed to download Boost")
        endif()

        message(STATUS "Extracting Boost...")
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xzf "boost_${BOOST_VERSION}.tar.gz"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
    endif()
endif()

set(BOOST_INCLUDE_DIR ${BOOST_DIR})
set(BOOST_SERIALIZATION_DIR ${BOOST_DIR}/libs/serialization/src)
set(BOOST_IOSTREAMS_DIR ${BOOST_DIR}/libs/iostreams/src)
set(BOOST_REGEX_DIR ${BOOST_DIR}/libs/regex/src)

# Collect Boost source files
file(GLOB BOOST_SERIALIZATION_SOURCES "${BOOST_SERIALIZATION_DIR}/*.cpp")
list(FILTER BOOST_SERIALIZATION_SOURCES EXCLUDE REGEX ".*xml.*")
list(FILTER BOOST_SERIALIZATION_SOURCES EXCLUDE REGEX ".*binary.*")

file(GLOB BOOST_IOSTREAMS_SOURCES "${BOOST_IOSTREAMS_DIR}/*.cpp")
list(FILTER BOOST_IOSTREAMS_SOURCES EXCLUDE REGEX ".*bzip.*")
list(FILTER BOOST_IOSTREAMS_SOURCES EXCLUDE REGEX ".*zstd.*")
list(FILTER BOOST_IOSTREAMS_SOURCES EXCLUDE REGEX ".*lzma.*")

file(GLOB BOOST_REGEX_SOURCES "${BOOST_REGEX_DIR}/*.cpp")

#------------------------------------------------------------------------------
# GSL sources
#------------------------------------------------------------------------------
set(GSL_SOURCES
    gsl/error.c
    gsl/sys/infnan.c
    gsl/sys/coerce.c
    gsl/sys/fdiv.c
    gsl/sys/pow_int.c
    gsl/sys/fcmp.c
    gsl/sys/log1p.c
    gsl/sys/invhyp.c
    gsl/sys/expm1.c
    gsl/complex/math.c
    gsl/specfunc/beta.c
    gsl/specfunc/psi.c
    gsl/specfunc/trig.c
    gsl/specfunc/exp.c
    gsl/specfunc/expint.c
    gsl/specfunc/log.c
    gsl/specfunc/erfc.c
    gsl/specfunc/zeta.c
    gsl/specfunc/elementary.c
    gsl/specfunc/gamma.c
    gsl/specfunc/gamma_inc.c
    gsl/rng/borosh13.c
    gsl/rng/fishman2x.c
    gsl/rng/mt.c
    gsl/rng/rand.c
    gsl/rng/ranmar.c
    gsl/rng/types.c
    gsl/rng/cmrg.c
    gsl/rng/gfsr4.c
    gsl/rng/r250.c
    gsl/rng/random.c
    gsl/rng/rng.c
    gsl/rng/uni32.c
    gsl/rng/coveyou.c
    gsl/rng/knuthran2.c
    gsl/rng/ran0.c
    gsl/rng/randu.c
    gsl/rng/slatec.c
    gsl/rng/uni.c
    gsl/rng/default.c
    gsl/rng/knuthran.c
    gsl/rng/ran1.c
    gsl/rng/ranf.c
    gsl/rng/taus113.c
    gsl/rng/vax.c
    gsl/rng/file.c
    gsl/rng/lecuyer21.c
    gsl/rng/ran2.c
    gsl/rng/ranlux.c
    gsl/rng/taus.c
    gsl/rng/waterman14.c
    gsl/rng/fishman18.c
    gsl/rng/minstd.c
    gsl/rng/ran3.c
    gsl/rng/ranlxd.c
    gsl/rng/transputer.c
    gsl/rng/zuf.c
    gsl/rng/fishman20.c
    gsl/rng/mrg.c
    gsl/rng/rand48.c
    gsl/rng/ranlxs.c
    gsl/rng/tt.c
    gsl/rng/knuthran2002.c
    gsl/randist/binomial.c
    gsl/randist/binomial_tpe.c
    gsl/randist/beta.c
    gsl/randist/exponential.c
    gsl/randist/geometric.c
    gsl/randist/nbinomial.c
    gsl/randist/poisson.c
    gsl/randist/multinomial.c
    gsl/randist/chisq.c
    gsl/randist/gauss.c
    gsl/randist/gausszig.c
    gsl/randist/gamma.c
    gsl/cdf/gammainv.c
    gsl/cdf/binomial.c
    gsl/cdf/beta.c
    gsl/cdf/betainv.c
    gsl/cdf/gauss.c
    gsl/cdf/gaussinv.c
    gsl/cdf/chisq.c
    gsl/cdf/chisqinv.c
    gsl/cdf/gamma.c
    gsl/cdf/poisson.c
    gsl/cdf/exponential.c
    gsl/cdf/exponentialinv.c
)

# GSL files for standalone gsl module (without RNG to avoid conflicts)
set(GSL_MODULE_SOURCES
    gsl/error.c
    gsl/sys/infnan.c
    gsl/sys/coerce.c
    gsl/sys/fdiv.c
    gsl/sys/pow_int.c
    gsl/sys/fcmp.c
    gsl/sys/log1p.c
    gsl/sys/invhyp.c
    gsl/sys/expm1.c
    gsl/complex/math.c
    gsl/specfunc/beta.c
    gsl/specfunc/elementary.c
    gsl/specfunc/erfc.c
    gsl/specfunc/exp.c
    gsl/specfunc/expint.c
    gsl/specfunc/log.c
    gsl/specfunc/psi.c
    gsl/specfunc/gamma.c
    gsl/specfunc/gamma_inc.c
    gsl/specfunc/trig.c
    gsl/specfunc/zeta.c
    gsl/cdf/beta.c
    gsl/cdf/betainv.c
    gsl/cdf/binomial.c
    gsl/cdf/gauss.c
    gsl/cdf/gaussinv.c
    gsl/cdf/exponential.c
    gsl/cdf/exponentialinv.c
    gsl/cdf/gamma.c
    gsl/cdf/gammainv.c
    gsl/cdf/chisq.c
    gsl/cdf/chisqinv.c
    gsl/cdf/poisson.c
)

#------------------------------------------------------------------------------
# simuPOP source files
#------------------------------------------------------------------------------
set(SIMUPOP_SOURCES
    src/simuPOP/utility.cpp
    src/simuPOP/genoStru.cpp
    src/simuPOP/individual.cpp
    src/simuPOP/population.cpp
    src/simuPOP/simulator.cpp
    src/simuPOP/mating.cpp
    src/simuPOP/operator.cpp
    src/simuPOP/initializer.cpp
    src/simuPOP/migrator.cpp
    src/simuPOP/outputer.cpp
    src/simuPOP/selector.cpp
    src/simuPOP/penetrance.cpp
    src/simuPOP/qtrait.cpp
    src/simuPOP/stator.cpp
    src/simuPOP/mutator.cpp
    src/simuPOP/transmitter.cpp
    src/simuPOP/tagger.cpp
    src/simuPOP/pedigree.cpp
    src/simuPOP/virtualSubPop.cpp
)

#------------------------------------------------------------------------------
# Common compile definitions
#------------------------------------------------------------------------------
set(COMMON_DEFINITIONS
    BOOST_UBLAS_NDEBUG
    BOOST_NO_CXX98_FUNCTION_BASE
    BOOST_COMPUTE_USE_CPP11
    _HAS_ITERATOR_DEBUGGING=0
    BOOST_ALL_NO_LIB
    NO_ZLIB=0
    NO_BZIP=1
    NO_ZSTD=1
    NO_LZMA=1
    TR1_SUPPORT=${TR1_SUPPORT}
)

if(USE_LIBCPP)
    list(APPEND COMMON_DEFINITIONS USE_LIBCPP=1)
endif()

if(USE_OPENMP AND OpenMP_FOUND AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    list(APPEND COMMON_DEFINITIONS _GLIBCXX_PARALLEL)
endif()

# Get version from _version.py
file(READ "${CMAKE_SOURCE_DIR}/src/simuPOP/_version.py" VERSION_FILE)
string(REGEX MATCH "__version__[^'\"]*['\"]([^'\"]+)['\"]" _ ${VERSION_FILE})
set(SIMUPOP_VER ${CMAKE_MATCH_1})
string(REGEX MATCH "__revision__[^'\"]*['\"]([^'\"]+)['\"]" _ ${VERSION_FILE})
set(SIMUPOP_REV ${CMAKE_MATCH_1})

list(APPEND COMMON_DEFINITIONS
    SIMUPOP_VER=${SIMUPOP_VER}
    SIMUPOP_REV=${SIMUPOP_REV}
)

#------------------------------------------------------------------------------
# Common include directories
#------------------------------------------------------------------------------
set(COMMON_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/gsl
    ${CMAKE_SOURCE_DIR}/gsl/specfunc
    ${CMAKE_SOURCE_DIR}/src/simuPOP
    ${CMAKE_BINARY_DIR}/generated
    ${BOOST_INCLUDE_DIR}
)

#------------------------------------------------------------------------------
# Build static support libraries
#------------------------------------------------------------------------------
set(LIB_SOURCES
    ${GSL_SOURCES}
    ${BOOST_SERIALIZATION_SOURCES}
    ${BOOST_IOSTREAMS_SOURCES}
    ${BOOST_REGEX_SOURCES}
)

# Standard library (with debug checks)
add_library(simuPOP_shstd STATIC ${LIB_SOURCES})
target_include_directories(simuPOP_shstd PRIVATE ${COMMON_INCLUDE_DIRS})
target_compile_definitions(simuPOP_shstd PRIVATE
    ${COMMON_DEFINITIONS}
    _SECURE_SCL=1
)
if(MSVC)
    target_compile_options(simuPOP_shstd PRIVATE /w)
else()
    target_compile_options(simuPOP_shstd PRIVATE -w)
endif()
target_link_libraries(simuPOP_shstd PRIVATE ZLIB::ZLIB)

# Optimized library (no debug checks)
add_library(simuPOP_shop STATIC ${LIB_SOURCES})
target_include_directories(simuPOP_shop PRIVATE ${COMMON_INCLUDE_DIRS})
target_compile_definitions(simuPOP_shop PRIVATE
    ${COMMON_DEFINITIONS}
    NDEBUG
    _SECURE_SCL=0
    OPTIMIZED
)
if(MSVC)
    target_compile_options(simuPOP_shop PRIVATE /w)
else()
    target_compile_options(simuPOP_shop PRIVATE -w)
endif()
target_link_libraries(simuPOP_shop PRIVATE ZLIB::ZLIB)

#------------------------------------------------------------------------------
# GSL Python module (using pre-generated wrapper)
#------------------------------------------------------------------------------
# Use pre-generated SWIG wrapper to avoid SWIG version compatibility issues
Python3_add_library(gsl MODULE
    src/simuPOP/gsl_wrap.c
    ${GSL_MODULE_SOURCES}
)
target_include_directories(gsl PRIVATE ${COMMON_INCLUDE_DIRS})
target_compile_definitions(gsl PRIVATE ${COMMON_DEFINITIONS})
if(MSVC)
    target_compile_options(gsl PRIVATE /w)
else()
    target_compile_options(gsl PRIVATE -w)
endif()
set_target_properties(gsl PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/simuPOP
    OUTPUT_NAME "_gsl"
    PREFIX ""
)

# Copy gsl.py to build directory
configure_file(src/simuPOP/gsl.py ${CMAKE_BINARY_DIR}/simuPOP/gsl.py COPYONLY)

#------------------------------------------------------------------------------
# Module variants configuration
#------------------------------------------------------------------------------
# Define all module variant names
set(MODULE_VARIANT_NAMES std op la laop ba baop mu muop lin linop)

# Define macros for each variant (use | as separator to avoid CMake list issues)
set(VARIANT_std_MACROS "")
set(VARIANT_std_LIB "STD")

set(VARIANT_op_MACROS "OPTIMIZED")
set(VARIANT_op_LIB "OPT")

set(VARIANT_la_MACROS "LONGALLELE")
set(VARIANT_la_LIB "STD")

set(VARIANT_laop_MACROS "LONGALLELE|OPTIMIZED")
set(VARIANT_laop_LIB "OPT")

set(VARIANT_ba_MACROS "BINARYALLELE")
set(VARIANT_ba_LIB "STD")

set(VARIANT_baop_MACROS "BINARYALLELE|OPTIMIZED")
set(VARIANT_baop_LIB "OPT")

set(VARIANT_mu_MACROS "MUTANTALLELE")
set(VARIANT_mu_LIB "STD")

set(VARIANT_muop_MACROS "MUTANTALLELE|OPTIMIZED")
set(VARIANT_muop_LIB "OPT")

set(VARIANT_lin_MACROS "LINEAGE")
set(VARIANT_lin_LIB "STD")

set(VARIANT_linop_MACROS "LINEAGE|OPTIMIZED")
set(VARIANT_linop_LIB "OPT")

# Process each module variant
foreach(VARIANT_NAME ${MODULE_VARIANT_NAMES})
    set(VARIANT_MACROS_STR "${VARIANT_${VARIANT_NAME}_MACROS}")
    set(VARIANT_LIB "${VARIANT_${VARIANT_NAME}_LIB}")

    # Convert pipe-separated macros to list
    string(REPLACE "|" ";" VARIANT_MACROS "${VARIANT_MACROS_STR}")

    # Use pre-generated SWIG wrapper file
    set(SWIG_WRAPPER "src/simuPOP/simuPOP_${VARIANT_NAME}_wrap.cpp")

    # Create the Python module using pre-generated wrapper
    Python3_add_library(simuPOP_${VARIANT_NAME} MODULE
        ${SWIG_WRAPPER}
        ${SIMUPOP_SOURCES}
    )

    # Module-specific compile definitions
    set(VARIANT_DEFINITIONS
        ${COMMON_DEFINITIONS}
        SIMUPOP_MODULE=simuPOP_${VARIANT_NAME}
    )

    if(VARIANT_MACROS)
        foreach(MACRO ${VARIANT_MACROS})
            list(APPEND VARIANT_DEFINITIONS ${MACRO})
        endforeach()
    endif()

    # Add STD or OPT specific macros
    if(VARIANT_LIB STREQUAL "OPT")
        list(APPEND VARIANT_DEFINITIONS NDEBUG _SECURE_SCL=0 OPTIMIZED)
        set(LINK_LIB simuPOP_shop)
    else()
        list(APPEND VARIANT_DEFINITIONS _SECURE_SCL=1)
        set(LINK_LIB simuPOP_shstd)
    endif()

    target_include_directories(simuPOP_${VARIANT_NAME} PRIVATE ${COMMON_INCLUDE_DIRS})
    target_compile_definitions(simuPOP_${VARIANT_NAME} PRIVATE ${VARIANT_DEFINITIONS})

    # On macOS, we need to force-load the static library to include all symbols
    if(APPLE)
        target_link_libraries(simuPOP_${VARIANT_NAME} PRIVATE
            -Wl,-force_load,$<TARGET_FILE:${LINK_LIB}>
            ZLIB::ZLIB
        )
    else()
        target_link_libraries(simuPOP_${VARIANT_NAME} PRIVATE
            ${LINK_LIB}
            ZLIB::ZLIB
        )
    endif()

    if(USE_OPENMP AND OpenMP_FOUND)
        target_link_libraries(simuPOP_${VARIANT_NAME} PRIVATE OpenMP::OpenMP_CXX)
    endif()

    set_target_properties(simuPOP_${VARIANT_NAME} PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/simuPOP
        OUTPUT_NAME "_simuPOP_${VARIANT_NAME}"
        PREFIX ""
    )

    # Add dependency on static library
    add_dependencies(simuPOP_${VARIANT_NAME} ${LINK_LIB})

    # Copy the Python wrapper file to build directory
    configure_file(src/simuPOP/simuPOP_${VARIANT_NAME}.py ${CMAKE_BINARY_DIR}/simuPOP/simuPOP_${VARIANT_NAME}.py COPYONLY)
endforeach()

#------------------------------------------------------------------------------
# Install rules
#------------------------------------------------------------------------------
# Install all Python extension modules
install(TARGETS gsl
    LIBRARY DESTINATION simuPOP
)

foreach(VARIANT_NAME ${MODULE_VARIANT_NAMES})
    install(TARGETS simuPOP_${VARIANT_NAME}
        LIBRARY DESTINATION simuPOP
    )
endforeach()

# Install Python source files
install(FILES
    src/simuPOP/__init__.py
    src/simuPOP/utils.py
    src/simuPOP/demography.py
    src/simuPOP/sampling.py
    src/simuPOP/simuOpt.py
    DESTINATION simuPOP
)

# Install Python wrapper files (from source, pre-generated by SWIG)
install(FILES src/simuPOP/gsl.py DESTINATION simuPOP)
foreach(VARIANT_NAME ${MODULE_VARIANT_NAMES})
    install(FILES src/simuPOP/simuPOP_${VARIANT_NAME}.py DESTINATION simuPOP)
endforeach()

# Install simuOpt.py at package root
install(FILES simuOpt.py DESTINATION .)

message(STATUS "=== simuPOP Build Configuration ===")
message(STATUS "Version: ${SIMUPOP_VER}")
message(STATUS "OpenMP: ${USE_OPENMP}")
message(STATUS "TR1 Support: ${TR1_SUPPORT}")
message(STATUS "Using libc++: ${USE_LIBCPP}")
message(STATUS "Boost directory: ${BOOST_DIR}")
message(STATUS "===================================")
