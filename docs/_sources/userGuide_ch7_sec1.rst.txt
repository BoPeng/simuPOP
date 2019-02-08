Module :mod:`simuOpt` (function :func:`simuOpt.setOptions`)
===========================================================

Module :mod:`simuOpt` handles options to specify which simuPOP module to load
and how this module should be loaded, using function ``simuOpt.setOptions``\
with parameters *alleleType* (``short``, ``long``, or ``binary`` ), *optimized*
(``standard`` or ``optimized``), *gui* (whether or not use a graphical user
interface and which graphical toolkit to use), *revision* ``(``\ minimal
required version/revision), *quiet* (with or without banner message, and *debug*
(which debug code to turn on). These options have been discussed in Example
:ref:`lst_Use_of_standard_module <lst_Use_of_standard_module>` and
:ref:`lst_Use_of_optimized_module <lst_Use_of_optimized_module>` and other
related sections. Note that **most options can be set by environmental variables
and command line options** which are sometimes more versatile to use.


