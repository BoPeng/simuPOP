/* CDF_ERROR: call the error handler, and return a NAN. */

#define CDF_ERROR(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, GSL_NAN)
// this is added for simuPOP compiling usage because scons can not
// tell between gsl/error.h and gsl/specfunc/error.h
#include <specfunc/error.h>

