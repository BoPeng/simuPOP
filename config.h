/**
 *  $File: config.h $
 *  $LastChangedDate: 2010-11-17 23:23:30 -0600 (Wed, 17 Nov 2010) $
 *  $Rev: 3896 $
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/** \file config.h
 *  This file includes appropriate configuration file for different
 *  operating systems and compilers.
 */

#ifdef _WIN32

/* define if bool is a built-in type */
#  define HAVE_BOOL

/* define if the compiler supports const_cast<> */
#  define HAVE_CONST_CAST

/* Define to 1 if you have the declaration of `acosh', and to 0 if you don't.
 */
#  define HAVE_DECL_ACOSH 1

/* Define to 1 if you have the declaration of `asinh', and to 0 if you don't.
 */
#  define HAVE_DECL_ASINH 1

/* Define to 1 if you have the declaration of `atanh', and to 0 if you don't.
 */
#  define HAVE_DECL_ATANH 1

/* Define to 1 if you have the declaration of `expm1', and to 0 if you don't.
 */
#  define HAVE_DECL_EXPM1 1

/* Define to 1 if you have the declaration of `feenableexcept', and to 0 if
   you don't. */
/* #undef HAVE_DECL_FEENABLEEXCEPT */

/* Define to 1 if you have the declaration of `fesettrapenable', and to 0 if
   you don't. */
/* #undef HAVE_DECL_FESETTRAPENABLE */

/* Define to 1 if you have the declaration of `finite', and to 0 if you don't.
 */
#  define HAVE_DECL_FINITE 1

/* Define to 1 if you have the declaration of `frexp', and to 0 if you don't.
 */
#  define HAVE_DECL_FREXP 1

/* Define to 1 if you have the declaration of `hypot', and to 0 if you don't.
 */
#  define HAVE_DECL_HYPOT 1

/* Define to 1 if you have the declaration of `isfinite', and to 0 if you
   don't. */
#  define HAVE_DECL_ISFINITE 0

/* Define to 1 if you have the declaration of `isinf', and to 0 if you don't.
 */
#  define HAVE_DECL_ISINF 1

/* Define to 1 if you have the declaration of `isnan', and to 0 if you don't.
 */
#  define HAVE_DECL_ISNAN 1

/* Define to 1 if you have the declaration of `ldexp', and to 0 if you don't.
 */
#  define HAVE_DECL_LDEXP 1

/* Define to 1 if you have the declaration of `log1p', and to 0 if you don't.
 */
/* #undef HAVE_DECL_LOG1P */

/* define if the compiler supports default template parameters */
#  define HAVE_DEFAULT_TEMPLATE_PARAMETERS

/* Define if /dev/null exists */
#  define HAVE_DEV_NULL 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#  define HAVE_DLFCN_H 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* define if the compiler supports dynamic_cast<> */
#  define HAVE_DYNAMIC_CAST

/* define if the compiler supports exceptions */
#  define HAVE_EXCEPTIONS

/* have exi success and failure */
/* #undef HAVE_EXIT_SUCCESS_AND_FAILURE */

/* "HAVE_EXTENDED_PRECISION_REGISTERS" */
/* #undef HAVE_EXTENDED_PRECISION_REGISTERS */

/* Define to 1 if you have the <float.h> header file. */
#  define HAVE_FLOAT_H 1

/* Define to 1 if you have the `floor' function. */
/* #undef HAVE_FLOOR */

/* "Define this is IEEE comparisons work correctly (e.g. NaN != NaN)" */
/* #undef HAVE_IEEE_COMPARISONS */

/* "Define this is IEEE denormalized numbers are available" */
/* #undef HAVE_IEEE_DENORMALS */

/* Define to 1 if you have the <inttypes.h> header file. */
/* #define HAVE_INTTYPES_H 1 */

/* Define to 1 if you have the `iswprint' function. */
#  define HAVE_ISWPRINT 1

/* Define to 1 if you have the `m' library (-lm). */
/* #undef HAVE_LIBM */

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#  define HAVE_MALLOC 1

/* Define to 1 if you have the `memcpy' function. */
/* #undef HAVE_MEMCPY */

/* Define to 1 if you have the `memmove' function. */
/* #undef HAVE_MEMMOVE */

/* Define to 1 if you have the <memory.h> header file. */
#  define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#  define HAVE_MEMSET 1

/* define if the compiler implements namespaces */
#  define HAVE_NAMESPACES

/* define if the compiler accepts the new for scoping rules */
#  define HAVE_NEW_FOR_SCOPING

/* Define to 1 if you have the `pow' function. */
/* #undef HAVE_POW */

/* "Define this if printf can handle %Lf for long double" */
/* #undef HAVE_PRINTF_LONGDOUBLE */

/* Define to 1 if the system has the type `ptrdiff_t'. */
#  define HAVE_PTRDIFF_T 1

/* Define to 1 if you have the `setenv' function. */
#  define HAVE_SETENV 1

/* Define to 1 if you have the `sqrt' function. */
/* #undef HAVE_SQRT */

/* Define to 1 if stdbool.h conforms to C99. */
#  define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stddef.h> header file. */
#  define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
// I have a portable stdint.h for msvc
#  define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#  define HAVE_STDLIB_H 1

/* define if the compiler supports Standard Template Library */
#  define HAVE_STL

/* Define to 1 if you have the `strdup' function. */
/* #undef HAVE_STRDUP */

/* Define to 1 if you have the <strings.h> header file. */
#  define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#  define HAVE_STRING_H 1

/* Define to 1 if you have the `strtol' function. */
/* #undef HAVE_STRTOL */

/* Define to 1 if you have the `strtoul' function. */
#  define HAVE_STRTOUL 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#  define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#  define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
// #define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vprintf' function. */
/* #undef HAVE_VPRINTF */

/* Define to 1 if the system has the type `_Bool'. */
#  define HAVE__BOOL 1

/* Define to the address where bug reports for this package should be sent. */
#  define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#  define PACKAGE_NAME ""

/* Define to the full name and version of this package. */
#  define PACKAGE_STRING ""

/* Define to the one symbol short name of this package. */
#  define PACKAGE_TARNAME ""

/* Define to the version of this package. */
#  define PACKAGE_VERSION ""

/* "Defined if this is an official release" */
/* #undef RELEASED */

/* Define to 1 if you have the ANSI C header files. */
#  define STDC_HEADERS 1
/* config.h.  Generated automatically by configure.  */
/* config.h.in.  Generated automatically from configure.in by autoheader.  */

/* Define to empty if the keyword does not work.  */
/* #undef const */

/* Define if you don't have vprintf but do have _doprnt.  */
/* #undef HAVE_DOPRNT */

/* Define if you have the vprintf function.  */
#  define HAVE_VPRINTF 1

/* Define as __inline if that's what the C compiler calls it.  */
#  define inline __inline

/* Define to `unsigned' if <sys/types.h> doesn't define.  */
/* #undef size_t */

/* Define if you have the ANSI C header files.  */
#  define STDC_HEADERS 1

/* Define if you have the acosh function.  */
/* #undef HAVE_ACOSH */
#define HAVE_ACOSH

/* Define if you have the asinh function.  */
/* #undef HAVE_ASINH */
#define HAVE_ASINH

/* Define if you have the atanh function.  */
/* #undef HAVE_ATANH */
#define HAVE_ATANH

/* Define if you have the expm1 function.  */
/* #undef HAVE_EXPM1 */
#define HAVE_EXPM1

/* Define if you have the finite function.  */
/* #undef HAVE_FINITE */

/* Define if you have the isfinite function.  */
/* #undef HAVE_ISFINITE */

/* Define if you have the isinf function.  */
/* #undef HAVE_ISINF */
#define HAVE_ISINF

/* Define if you have the isnan function.  */
/* #undef HAVE_ISNAN */
#define HAVE_ISNAN

/* Define if you have the log1p function.  */
/* #undef HAVE_LOG1P */
#define HAVE_LOG1P

/* Define if you have the memcpy function.  */
#  define HAVE_MEMCPY 1

/* Define if you have the memmove function.  */
#  define HAVE_MEMMOVE 1

/* Define if you have the strdup function.  */
#  define HAVE_STRDUP 1

/* Define if you have the strtol function.  */
#  define HAVE_STRTOL 1

/* Define if you have the strtoul function.  */
#  define HAVE_STRTOUL 1

/* Define if you have the <dlfcn.h> header file.  */
/* #undef HAVE_DLFCN_H */

/* Define if you have the m library (-lm).  */
/* #undef HAVE_LIBM */

/* Name of package */
#  define PACKAGE "gsl"

/* Version number of package */
#  define VERSION "1.0"

/* Define if you have inline */
#  define HAVE_INLINE 1

/* Define if you need to hide the static definitions of inline functions */
#  define HIDE_INLINE_STATIC 1

/* Define if you have the ansi CLOCKS_PER_SEC clock rate */
#  define HAVE_CLOCKS_PER_SEC 1

/* Defined if configure has guessed a missing ansi CLOCKS_PER_SEC clock rate */
/* #undef HAVE_GUESSED_CLOCKS_PER_SEC */

/* Use configure's best guess for CLOCKS_PER_SEC if it is unknown */
#  ifndef HAVE_CLOCKS_PER_SEC
#    define CLOCKS_PER_SEC HAVE_GUESSED_CLOCKS_PER_SEC
#  endif

/* Defined if you have ansi EXIT_SUCCESS and EXIT_FAILURE in stdlib.h */
#  define HAVE_EXIT_SUCCESS_AND_FAILURE 1

/* Use 0 and 1 for EXIT_SUCCESS and EXIT_FAILURE if we don't have them */
#  ifndef HAVE_EXIT_SUCCESS_AND_FAILURE
#    define EXIT_SUCCESS 0
#    define EXIT_FAILURE 1
#  endif

/* Define one of these if you have a known IEEE arithmetic interface */
/* #undef HAVE_SPARCLINUX_IEEE_INTERFACE */
/* #undef HAVE_M68KLINUX_IEEE_INTERFACE */
/* #undef HAVE_PPCLINUX_IEEE_INTERFACE */
/* #undef HAVE_X86LINUX_IEEE_INTERFACE */
/* #undef HAVE_SUNOS4_IEEE_INTERFACE */
/* #undef HAVE_SOLARIS_IEEE_INTERFACE */
/* #undef HAVE_HPUX11_IEEE_INTERFACE */
/* #undef HAVE_HPUX_IEEE_INTERFACE */
/* #undef HAVE_TRU64_IEEE_INTERFACE */
/* #undef HAVE_IRIX_IEEE_INTERFACE */
/* #undef HAVE_AIX_IEEE_INTERFACE */
/* #undef HAVE_FREEBSD_IEEE_INTERFACE */
/* #undef HAVE_OS2EMX_IEEE_INTERFACE */
/* #undef HAVE_NETBSD_IEEE_INTERFACE */
/* #undef HAVE_OPENBSD_IEEE_INTERFACE */
/* #undef HAVE_DARWIN_IEEE_INTERFACE */

/* Define this if we need to include /usr/include/float.h explicitly
   in order to get FP_RND_RN and related macros.  This is known to be
   a problem on some Compaq Tru64 unix systems when compiled with GCC. */
/* #undef FIND_FP_RND_IN_USR_INCLUDE_FLOAT_H */

/* Define a rounding function which moves extended precision values
   out of registers and rounds them to double-precision. This should
   be used *sparingly*, in places where it is necessary to keep
   double-precision rounding for critical expressions while running in
   extended precision. For example, the following code should ensure
   exact equality, even when extended precision registers are in use,

      double q = GSL_COERCE_DBL(3.0/7.0) ;
      if (q == GSL_COERCE_DBL(3.0/7.0)) { ... } ;

   It carries a penalty even when the program is running in double
   precision mode unless you compile a separate version of the
   library with HAVE_EXTENDED_PRECISION_REGISTERS turned off. */

#  define HAVE_EXTENDED_PRECISION_REGISTERS 1

#  ifdef HAVE_EXTENDED_PRECISION_REGISTERS
#    define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#  else
#    define GSL_COERCE_DBL(x) (x)
#  endif

/* Define this if printf can handle %Lf for long double */
/* #undef HAVE_PRINTF_LONGDOUBLE */

/* Substitute gsl functions for missing system functions */

#  ifndef HAVE_HYPOT
#    ifndef _MSC_VER
#      define hypot gsl_hypot
#    endif
#  endif

#  ifndef HAVE_LOG1P
#    define log1p gsl_log1p
#  endif

#  ifndef HAVE_EXPM1
#    define expm1 gsl_expm1
#  endif

#  ifndef HAVE_ACOSH
#    define acosh gsl_acosh
#  endif

#  ifndef HAVE_ASINH
#    define asinh gsl_asinh
#  endif

#  ifndef HAVE_ATANH
#    define atanh gsl_atanh
#  endif

#  ifndef HAVE_ISINF
#    define isinf gsl_isinf
#  endif

#  ifndef HAVE_ISNAN
#    define isnan gsl_isnan
#  endif

#  ifndef HAVE_FINITE
#    ifdef HAVE_ISFINITE
#      define finite isfinite
#    else
#      define finite gsl_finite
#    endif
#  endif

#  define RETURN_IF_NULL(x) if (!x) { return ; }

#else
#  ifdef MACOSX

/* define if bool is a built-in type */
#    define HAVE_BOOL

/* define if the compiler supports const_cast<> */
#    define HAVE_CONST_CAST

/* Define to 1 if you have the declaration of `acosh', and to 0 if you don't.
 */
#    define HAVE_DECL_ACOSH 1

/* Define to 1 if you have the declaration of `atanh', and to 0 if you don't.
 */
#    define HAVE_DECL_ATANH 1

/* Define to 1 if you have the declaration of `expm1', and to 0 if you don't.
 */
#    define HAVE_DECL_EXPM1 1

/* Define to 1 if you have the declaration of `feenableexcept', and to 0 if
   you don't. */
/* #undef HAVE_DECL_FEENABLEEXCEPT */

/* Define to 1 if you have the declaration of `fesettrapenable', and to 0 if
   you don't. */
/* #undef HAVE_DECL_FESETTRAPENABLE */

/* Define to 1 if you have the declaration of `finite', and to 0 if you don't.
 */
#    define HAVE_DECL_FINITE 1

/* Define to 1 if you have the declaration of `frexp', and to 0 if you don't.
 */
#    define HAVE_DECL_FREXP 1

/* Define to 1 if you have the declaration of `hypot', and to 0 if you don't.
 */
#    define HAVE_DECL_HYPOT 1

/* Define to 1 if you have the declaration of `isfinite', and to 0 if you
   don't. */
#    define HAVE_DECL_ISFINITE 1

/* Define to 1 if you have the declaration of `isinf', and to 0 if you don't.
 */
#    define HAVE_DECL_ISINF 1

/* Define to 1 if you have the declaration of `isnan', and to 0 if you don't.
 */
#    define HAVE_DECL_ISNAN 1

/* Define to 1 if you have the declaration of `ldexp', and to 0 if you don't.
 */
#    define HAVE_DECL_LDEXP 1

/* Define to 1 if you have the declaration of `log1p', and to 0 if you don't.
 */
/* #undef HAVE_DECL_LOG1P */

/* define if the compiler supports default template parameters */
#    define HAVE_DEFAULT_TEMPLATE_PARAMETERS

/* Define if /dev/null exists */
#    define HAVE_DEV_NULL 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#    define HAVE_DLFCN_H 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* define if the compiler supports dynamic_cast<> */
#    define HAVE_DYNAMIC_CAST

/* define if the compiler supports exceptions */
#    define HAVE_EXCEPTIONS

/* have exi success and failure */
/* #undef HAVE_EXIT_SUCCESS_AND_FAILURE */

/* "HAVE_EXTENDED_PRECISION_REGISTERS" */
/* #undef HAVE_EXTENDED_PRECISION_REGISTERS */

/* Define to 1 if you have the <float.h> header file. */
#    define HAVE_FLOAT_H 1

/* Define to 1 if you have the `floor' function. */
#    define HAVE_FLOOR 1

/* "Define this is IEEE comparisons work correctly (e.g. NaN != NaN)" */
/* #undef HAVE_IEEE_COMPARISONS */

/* "Define this is IEEE denormalized numbers are available" */
/* #undef HAVE_IEEE_DENORMALS */

/* "Define if you have inline" */
#    define HAVE_INLINE 1

/* Define to 1 if you have the <inttypes.h> header file. */
#    define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `iswprint' function. */
#    define HAVE_ISWPRINT 1

/* Define to 1 if you have the `m' library (-lm). */
/* #undef HAVE_LIBM */

/* Define to 1 if you have the <limits.h> header file. */
#    define HAVE_LIMITS_H 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#    define HAVE_MALLOC 1

/* Define to 1 if you have the `memcpy' function. */
/* #undef HAVE_MEMCPY */

/* Define to 1 if you have the `memmove' function. */
/* #undef HAVE_MEMMOVE */

/* Define to 1 if you have the <memory.h> header file. */
#    define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#    define HAVE_MEMSET 1

/* define if the compiler implements namespaces */
#    define HAVE_NAMESPACES

/* define if the compiler accepts the new for scoping rules */
#    define HAVE_NEW_FOR_SCOPING

/* Define to 1 if you have the `pow' function. */
#    define HAVE_POW 1

/* "Define this if printf can handle %Lf for long double" */
/* #undef HAVE_PRINTF_LONGDOUBLE */

/* Define to 1 if the system has the type `ptrdiff_t'. */
#    define HAVE_PTRDIFF_T 1

/* Define to 1 if you have the `snprintf' function. */
#    define HAVE_SNPRINTF 1

/* Define to 1 if you have the `sqrt' function. */
#    define HAVE_SQRT 1

/* Define to 1 if stdbool.h conforms to C99. */
#    define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stddef.h> header file. */
#    define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#    define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#    define HAVE_STDLIB_H 1

/* define if the compiler supports Standard Template Library */
#    define HAVE_STL

/* Define to 1 if you have the `strdup' function. */
/* #undef HAVE_STRDUP */

/* Define to 1 if you have the <strings.h> header file. */
#    define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#    define HAVE_STRING_H 1

/* Define to 1 if you have the `strtol' function. */
/* #undef HAVE_STRTOL */

/* Define to 1 if you have the `strtoul' function. */
#    define HAVE_STRTOUL 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#    define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#    define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#    define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vprintf' function. */
/* #undef HAVE_VPRINTF */

/* Define to 1 if the system has the type `_Bool'. */
#    define HAVE__BOOL 1

/* Define to the address where bug reports for this package should be sent. */
#    define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#    define PACKAGE_NAME ""

/* Define to the full name and version of this package. */
#    define PACKAGE_STRING ""

/* Define to the one symbol short name of this package. */
#    define PACKAGE_TARNAME ""

/* Define to the version of this package. */
#    define PACKAGE_VERSION ""

/* "Defined if this is an official release" */
/* #undef RELEASED */

/* Define to 1 if you have the ANSI C header files. */
#    define STDC_HEADERS 1

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#    ifndef __cplusplus
/* #undef inline */
#    endif

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to `unsigned' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */

#    if HAVE_EXTENDED_PRECISION_REGISTERS
#      define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#    else
#      define GSL_COERCE_DBL(x) (x)
#    endif

/* Substitute gsl functions for missing system functions */

#    if !HAVE_DECL_HYPOT
#      define hypot gsl_hypot
#    endif

#    if !HAVE_DECL_LOG1P
#      define log1p gsl_log1p
#    endif

#    if !HAVE_DECL_EXPM1
#      define expm1 gsl_expm1
#    endif

#    if !HAVE_DECL_ACOSH
#      define acosh gsl_acosh
#    endif

#    if !HAVE_DECL_ASINH
#      define asinh gsl_asinh
#    endif

#    if !HAVE_DECL_ATANH
#      define atanh gsl_atanh
#    endif

#    if !HAVE_DECL_LDEXP
#      define ldexp gsl_ldexp
#    endif

#    if !HAVE_DECL_FREXP
#      define frexp gsl_frexp
#    endif

#    if !HAVE_DECL_ISINF
#      define isinf gsl_isinf
#    endif

#    if !HAVE_DECL_FINITE
#      if HAVE_DECL_ISFINITE
#        define finite isfinite
#      else
#        define finite gsl_finite
#      endif
#    endif

#    if !HAVE_DECL_ISNAN
#      define isnan gsl_isnan
#    endif

#    if defined(GSL_RANGE_CHECK_OFF) || !defined(GSL_RANGE_CHECK)
#      define GSL_RANGE_CHECK 0                   /* turn off range checking by default internally */
#    endif

#    define RETURN_IF_NULL(x) if (!x) { return ; }
#  else
#    ifdef __sparc__

/* define if bool is a built-in type */
#      define HAVE_BOOL

/* define if the compiler supports const_cast<> */
#      define HAVE_CONST_CAST

/* Define to 1 if you have the declaration of `acosh', and to 0 if you don't.
 */
#      define HAVE_DECL_ACOSH 1

/* Define to 1 if you have the declaration of `atanh', and to 0 if you don't.
 */
#      define HAVE_DECL_ATANH 1

/* Define to 1 if you have the declaration of `expm1', and to 0 if you don't.
 */
#      define HAVE_DECL_EXPM1 1

/* Define to 1 if you have the declaration of `feenableexcept', and to 0 if
   you don't. */
#      define HAVE_DECL_FEENABLEEXCEPT 0

/* Define to 1 if you have the declaration of `fesettrapenable', and to 0 if
   you don't. */
#      define HAVE_DECL_FESETTRAPENABLE 0

/* Define to 1 if you have the declaration of `finite', and to 0 if you don't.
 */
#      define HAVE_DECL_FINITE 0

/* Define to 1 if you have the declaration of `frexp', and to 0 if you don't.
 */
#      define HAVE_DECL_FREXP 1

/* Define to 1 if you have the declaration of `hypot', and to 0 if you don't.
 */
#      define HAVE_DECL_HYPOT 1

/* Define to 1 if you have the declaration of `isfinite', and to 0 if you
   don't. */
#      define HAVE_DECL_ISFINITE 0

/* Define to 1 if you have the declaration of `isinf', and to 0 if you don't.
 */
#      define HAVE_DECL_ISINF 0

/* Define to 1 if you have the declaration of `isnan', and to 0 if you don't.
 */
#      define HAVE_DECL_ISNAN 1

/* Define to 1 if you have the declaration of `ldexp', and to 0 if you don't.
 */
#      define HAVE_DECL_LDEXP 1

/* Define to 1 if you have the declaration of `log1p', and to 0 if you don't.
 */
#      define HAVE_DECL_LOG1P 1

/* define if the compiler supports default template parameters */
#      define HAVE_DEFAULT_TEMPLATE_PARAMETERS

/* Define if /dev/null exists */
#      define HAVE_DEV_NULL 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#      define HAVE_DLFCN_H 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
#      define HAVE_DOPRNT 1

/* define if the compiler supports dynamic_cast<> */
#      define HAVE_DYNAMIC_CAST

/* define if the compiler supports exceptions */
#      define HAVE_EXCEPTIONS

/* have exi success and failure */
#      define HAVE_EXIT_SUCCESS_AND_FAILURE

/* "HAVE_EXTENDED_PRECISION_REGISTERS" */
/* #undef HAVE_EXTENDED_PRECISION_REGISTERS */

/* Define to 1 if you have the <float.h> header file. */
#      define HAVE_FLOAT_H 1

/* Define to 1 if you have the `floor' function. */
/* #undef HAVE_FLOOR */

/* "Define this is IEEE comparisons work correctly (e.g. NaN != NaN)" */
#      define HAVE_IEEE_COMPARISONS 1

/* "Define this is IEEE denormalized numbers are available" */
#      define HAVE_IEEE_DENORMALS 1

/* "Define if you have inline" */
#      define HAVE_INLINE 1

/* Define to 1 if you have the <inttypes.h> header file. */
#      define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `iswprint' function. */
#      define HAVE_ISWPRINT 1

/* Define to 1 if you have the `m' library (-lm). */
#      define HAVE_LIBM 1

/* Define to 1 if you have the <limits.h> header file. */
#      define HAVE_LIMITS_H 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#      define HAVE_MALLOC 1

/* Define to 1 if you have the `memcpy' function. */
#      define HAVE_MEMCPY 1

/* Define to 1 if you have the `memmove' function. */
#      define HAVE_MEMMOVE 1

/* Define to 1 if you have the <memory.h> header file. */
#      define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#      define HAVE_MEMSET 1

/* define if the compiler implements namespaces */
#      define HAVE_NAMESPACES

/* define if the compiler accepts the new for scoping rules */
#      define HAVE_NEW_FOR_SCOPING

/* Define to 1 if you have the `pow' function. */
/* #undef HAVE_POW */

/* "Define this if printf can handle %Lf for long double" */
#      define HAVE_PRINTF_LONGDOUBLE 1

/* Define to 1 if the system has the type `ptrdiff_t'. */
#      define HAVE_PTRDIFF_T 1

/* Define to 1 if you have the `snprintf' function. */
#      define HAVE_SNPRINTF 1

/* Define to 1 if you have the `sqrt' function. */
/* #undef HAVE_SQRT */

/* Define to 1 if stdbool.h conforms to C99. */
#      define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stddef.h> header file. */
#      define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
/* #undef HAVE_STDINT_H */

/* Define to 1 if you have the <stdlib.h> header file. */
#      define HAVE_STDLIB_H 1

/* define if the compiler supports Standard Template Library */
#      define HAVE_STL

/* Define to 1 if you have the `strdup' function. */
#      define HAVE_STRDUP 1

/* Define to 1 if you have the <strings.h> header file. */
#      define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#      define HAVE_STRING_H 1

/* Define to 1 if you have the `strtol' function. */
#      define HAVE_STRTOL 1

/* Define to 1 if you have the `strtoul' function. */
#      define HAVE_STRTOUL 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#      define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#      define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#      define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vprintf' function. */
#      define HAVE_VPRINTF 1

/* Define to 1 if the system has the type `_Bool'. */
#      define HAVE__BOOL 1

/* Define to the address where bug reports for this package should be sent. */
#      define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#      define PACKAGE_NAME ""

/* Define to the full name and version of this package. */
#      define PACKAGE_STRING ""

/* Define to the one symbol short name of this package. */
#      define PACKAGE_TARNAME ""

/* Define to the version of this package. */
#      define PACKAGE_VERSION ""

/* "Defined if this is an official release" */
#      define RELEASED 1

/* Define to 1 if you have the ANSI C header files. */
#      define STDC_HEADERS 1

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#      ifndef __cplusplus
/* #undef inline */
#      endif

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to `unsigned' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */

#      if HAVE_EXTENDED_PRECISION_REGISTERS
#        define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#      else
#        define GSL_COERCE_DBL(x) (x)
#      endif

/* Substitute gsl functions for missing system functions */

#      if !HAVE_DECL_HYPOT
#        define hypot gsl_hypot
#      endif

#      if !HAVE_DECL_LOG1P
#        define log1p gsl_log1p
#      endif

#      if !HAVE_DECL_EXPM1
#        define expm1 gsl_expm1
#      endif

#      if !HAVE_DECL_ACOSH
#        define acosh gsl_acosh
#      endif

#      if !HAVE_DECL_ASINH
#        define asinh gsl_asinh
#      endif

#      if !HAVE_DECL_ATANH
#        define atanh gsl_atanh
#      endif

#      if !HAVE_DECL_LDEXP
#        define ldexp gsl_ldexp
#      endif

#      if !HAVE_DECL_FREXP
#        define frexp gsl_frexp
#      endif

#      if !HAVE_DECL_ISINF
#        define isinf gsl_isinf
#      endif

#      if !HAVE_DECL_FINITE
#        if HAVE_DECL_ISFINITE
#          define finite isfinite
#        else
#          define finite gsl_finite
#        endif
#      endif

#      if !HAVE_DECL_ISNAN
#        define isnan gsl_isnan
#      endif

#      if defined(GSL_RANGE_CHECK_OFF) || !defined(GSL_RANGE_CHECK)
#        define GSL_RANGE_CHECK 0                 /* turn off range checking by default internally */
#      endif

#      define RETURN_IF_NULL(x) if (!x) { return ; }
#    else

/* define if bool is a built-in type */
#      define HAVE_BOOL

/* define if the compiler supports const_cast<> */
#      define HAVE_CONST_CAST

#      define HAVE_DECL_LOG1P 1
#      define HAVE_DECL_ASINH 1

/* Define to 1 if you have the declaration of `acosh', and to 0 if you don't.
 */
#      define HAVE_DECL_ACOSH 1

/* Define to 1 if you have the declaration of `atanh', and to 0 if you don't.
 */
#      define HAVE_DECL_ATANH 1

/* Define to 1 if you have the declaration of `expm1', and to 0 if you don't.
 */
#      define HAVE_DECL_EXPM1 1

/* Define to 1 if you have the declaration of `feenableexcept', and to 0 if
   you don't. */
/* #undef HAVE_DECL_FEENABLEEXCEPT */

/* Define to 1 if you have the declaration of `fesettrapenable', and to 0 if
   you don't. */
/* #undef HAVE_DECL_FESETTRAPENABLE */

/* Define to 1 if you have the declaration of `finite', and to 0 if you don't.
 */
#      define HAVE_DECL_FINITE 1

/* Define to 1 if you have the declaration of `frexp', and to 0 if you don't.
 */
#      define HAVE_DECL_FREXP 1

/* Define to 1 if you have the declaration of `hypot', and to 0 if you don't.
 */
#      define HAVE_DECL_HYPOT 1

/* Define to 1 if you have the declaration of `isfinite', and to 0 if you
   don't. */
#      ifndef HAVE_DECL_ISFINITE
#        define HAVE_DECL_ISFINITE 0
#      endif

/* Define to 1 if you have the declaration of `isinf', and to 0 if you don't.
 */
#      define HAVE_DECL_ISINF 1

/* Define to 1 if you have the declaration of `isnan', and to 0 if you don't.
 */
#      define HAVE_DECL_ISNAN 1

/* Define to 1 if you have the declaration of `ldexp', and to 0 if you don't.
 */
#      define HAVE_DECL_LDEXP 1

/* Define to 1 if you have the declaration of `log1p', and to 0 if you don't.
 */
/* #undef HAVE_DECL_LOG1P */

/* define if the compiler supports default template parameters */
#      define HAVE_DEFAULT_TEMPLATE_PARAMETERS

/* Define if /dev/null exists */
#      define HAVE_DEV_NULL 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#      define HAVE_DLFCN_H 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* define if the compiler supports dynamic_cast<> */
#      define HAVE_DYNAMIC_CAST

/* define if the compiler supports exceptions */
#      define HAVE_EXCEPTIONS

/* have exi success and failure */
/* #undef HAVE_EXIT_SUCCESS_AND_FAILURE */

/* "HAVE_EXTENDED_PRECISION_REGISTERS" */
/* #undef HAVE_EXTENDED_PRECISION_REGISTERS */

/* Define to 1 if you have the <float.h> header file. */
#      define HAVE_FLOAT_H 1

/* Define to 1 if you have the `floor' function. */
/* #undef HAVE_FLOOR */

/* "Define this is IEEE comparisons work correctly (e.g. NaN != NaN)" */
/* #undef HAVE_IEEE_COMPARISONS */

/* "Define this is IEEE denormalized numbers are available" */
/* #undef HAVE_IEEE_DENORMALS */

/* "Define if you have inline" */
#      define HAVE_INLINE 1

/* Define to 1 if you have the <inttypes.h> header file. */
#      define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `iswprint' function. */
#      define HAVE_ISWPRINT 1

/* Define to 1 if you have the `m' library (-lm). */
/* #undef HAVE_LIBM */

/* Define to 1 if you have the <limits.h> header file. */
#      define HAVE_LIMITS_H 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#      define HAVE_MALLOC 1

/* Define to 1 if you have the `memcpy' function. */
/* #undef HAVE_MEMCPY */

/* Define to 1 if you have the `memmove' function. */
/* #undef HAVE_MEMMOVE */

/* Define to 1 if you have the <memory.h> header file. */
#      define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#      define HAVE_MEMSET 1

/* define if the compiler implements namespaces */
#      define HAVE_NAMESPACES

/* define if the compiler accepts the new for scoping rules */
#      define HAVE_NEW_FOR_SCOPING

/* Define to 1 if you have the `pow' function. */
/* #undef HAVE_POW */

/* "Define this if printf can handle %Lf for long double" */
/* #undef HAVE_PRINTF_LONGDOUBLE */

/* Define to 1 if the system has the type `ptrdiff_t'. */
#      define HAVE_PTRDIFF_T 1

/* Define to 1 if you have the `snprintf' function. */
#      define HAVE_SNPRINTF 1

/* Define to 1 if you have the `sqrt' function. */
/* #undef HAVE_SQRT */

/* Define to 1 if stdbool.h conforms to C99. */
#      define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stddef.h> header file. */
#      define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#      define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#      define HAVE_STDLIB_H 1

/* define if the compiler supports Standard Template Library */
#      define HAVE_STL

/* Define to 1 if you have the `strdup' function. */
/* #undef HAVE_STRDUP */

/* Define to 1 if you have the <strings.h> header file. */
#      define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#      define HAVE_STRING_H 1

/* Define to 1 if you have the `strtol' function. */
/* #undef HAVE_STRTOL */

/* Define to 1 if you have the `strtoul' function. */
#      define HAVE_STRTOUL 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#      define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#      define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#      define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vprintf' function. */
/* #undef HAVE_VPRINTF */

/* Define to 1 if the system has the type `_Bool'. */
#      define HAVE__BOOL 1

/* Define to the address where bug reports for this package should be sent. */
#      define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#      define PACKAGE_NAME ""

/* Define to the full name and version of this package. */
#      define PACKAGE_STRING ""

/* Define to the one symbol short name of this package. */
#      define PACKAGE_TARNAME ""

/* Define to the version of this package. */
#      define PACKAGE_VERSION ""

/* "Defined if this is an official release" */
/* #undef RELEASED */

/* Define to 1 if you have the ANSI C header files. */
#      define STDC_HEADERS 1

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#      ifndef __cplusplus
/* #undef inline */
#      endif

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to `unsigned' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */

#      if HAVE_EXTENDED_PRECISION_REGISTERS
#        define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#      else
#        define GSL_COERCE_DBL(x) (x)
#      endif

/* Substitute gsl functions for missing system functions */

#      if !HAVE_DECL_HYPOT
#        define hypot gsl_hypot
#      endif

#      if !HAVE_DECL_LOG1P
#        define log1p gsl_log1p
#      endif

#      if !HAVE_DECL_EXPM1
#        define expm1 gsl_expm1
#      endif

#      if !HAVE_DECL_ACOSH
#        define acosh gsl_acosh
#      endif

#      if !HAVE_DECL_ASINH
#        define asinh gsl_asinh
#      endif

#      if !HAVE_DECL_ATANH
#        define atanh gsl_atanh
#      endif

#      if !HAVE_DECL_LDEXP
#        define ldexp gsl_ldexp
#      endif

#      if !HAVE_DECL_FREXP
#        define frexp gsl_frexp
#      endif

#      if !HAVE_DECL_ISINF
#        define isinf gsl_isinf
#      endif

#      if !HAVE_DECL_FINITE
#        if HAVE_DECL_ISFINITE
#          define finite isfinite
#        else
#          define finite gsl_finite
#        endif
#      endif

#      if !HAVE_DECL_ISNAN
#        define isnan gsl_isnan
#      endif

#      if defined(GSL_RANGE_CHECK_OFF) || !defined(GSL_RANGE_CHECK)
#        define GSL_RANGE_CHECK 0                 /* turn off range checking by default internally */
#      endif

#      define RETURN_IF_NULL(x) if (!x) { return ; }
#    endif
#  endif
#endif
