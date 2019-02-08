Module :mod:`simuPOP.gsl`
=========================

simuPOP makes use of many functions from the GUN Scientific Library. These
functions are used to generate random number and perform statistical tests
within simuPOP. Although these functions are not part of simuPOP, they can be
useful to users of simuPOP from time to time and it makes sense to expose these
functions directly to users.

Module :mod:`simuPOP.gsl` contains a number of GSL functions. Because only a
small proportion of GSL functions are used in simuPOP, this module is by no
means a comprehensive wrapper of GSL. Please refer to the simuPOP reference
manual for a list of functions included in this module, and the GSL manual for
more details. Because random number generation functions such as
``gsl_ran_gamma`` are already provided in the ``simuPOP.RNG`` class (e.g.
:func:`getRNG`\ ``.randGamma``), they are not provided in this module.


