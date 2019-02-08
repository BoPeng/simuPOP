Module :mod:`simuOpt`
=====================


.. module:: simuOpt

Module ``simuOpt`` provides a function ``simuOpt.setOptions`` to control which
simuPOP module to load, and how it is loaded, and a class ``simuOpt.Params``
that helps users manage simulation parameters.

When simuPOP is loaded, it checkes a few environmental variables
(``SIMUOPTIMIZED``, ``SIMUALLELETYPE``, and ``SIMUDEBUG``) to determine which
simuPOP module to load, and how to load it. More options can be set using the
``simuOpt.setOptions`` function. For example, you can suppress the banner
message when simuPOP is loaded and require a minimal version of simuPOP for
your script. simuPOP recognize the following commandline arguments

``--optimized``
    Load the optimized version of a simuPOP module.

``--gui=None|batch|interactive|True|wxPython|Tkinter``
    Whether or not use a graphical toolkit and which one to use.
    ``--gui=batch`` is usually used to run a script in batch mode (do not start
    a parameter input dialog and use all default values unless a parameter is
    specified from command line or a configuraiton file. If
    ``--gui=interactive``, an interactive shell will be used to solicit input
    from users. Otherwise, simuPOP will try to use a graphical parameter input
    dialog, and falls to an interactive mode when no graphical Toolkit is
    available. Please refer to parameter ``gui`` for ``simuOpt.setOptions``
    for details.

class ``params.Params`` provides a powerful way to handle commandline
arguments. Briefly speaking, a ``Params`` object can be created from a list
of parameter specification dictionaries. The parameters are then become
attributes of this object. A number of functions are provided to determine
values of these parameters using commandline arguments, a configuration
file, or a parameter input dialog (using ``Tkinter`` or ``wxPython``).
Values of these parameters can be accessed as attributes, or extracted
as a list or a dictionary. Note that the ``Params.getParam`` function
automatically handles the following commandline arguments.

``-h`` or ``--help``
    Print usage message.

``--config=configFile``
    Read parameters from a configuration file *configFile*.


Function setOptions
-------------------


.. function:: setOptions(alleleType=None, optimized=None, gui=None, quiet=None, debug=None, version=None, revision=None, numThreads=None, plotter=None)

   Set options before simuPOP is loaded to control which simuPOP module to
   load, and how the module should be loaded.
   
   alleleType
       Use the standard, binary,long or mutant allele version of the simuPOP
       module if ``alleleType`` is set to 'short', 'binary', 'long', 'mutant',
       or 'lineage' respectively. If this parameter is not set, this function
       will try to get its value from environmental variable ``SIMUALLELETYPE``.
       The standard (short) module will be used if the environmental variable
       is not defined.
   
   optimized
       Load the optimized version of a module if this parameter is set to
       ``True`` and the standard version if it is set to ``False``. If this
       parameter is not set (``None``), the optimized version will be used
       if environmental variable ``SIMUOPTIMIZED`` is defined. The standard
       version will be used otherwise.
   
   gui
       Whether or not use graphical user interfaces, which graphical toolkit
       to use and how to process parameters in non-GUI mode. If this parameter
       is ``None`` (default), this function will check environmental variable
       ``SIMUGUI`` or commandline option ``--gui`` for a value, and assume
       ``True`` if such an option is unavailable. If ``gui=True``, simuPOP
       will use ``wxPython``-based dialogs if ``wxPython`` is available, and
       use ``Tkinter``-based dialogs if ``Tkinter`` is available and use an
       interactive shell if no graphical toolkit is available.
       ``gui='Tkinter'`` or ``'wxPython'`` can be used to specify the
       graphical toolkit to use. If ``gui='interactive'``, a simuPOP script
       prompt users to input values of parameters. If ``gui='batch'``,
       default values of unspecified parameters will be used. In any case,
       commandline arguments and a configuration file specified by parameter
       --config will be processed. This option is usually left to ``None`` so
       that the same script can be run in both GUI and batch mode using
       commandline option ``--gui``.
   
   plotter
       (Deprecated)
   
   quiet
       If set to ``True``, suppress the banner message when a simuPOP module
       is loaded.
   
   debug
       A list of debug code (as string) that will be turned on when simuPOP
       is loaded. If this parameter is not set, a list of comma separated
       debug code specified in environmental variable ``SIMUDEBUG``, if
       available, will be used. Note that setting ``debug=[]`` will remove
       any debug code that might have been by variable ``SIMUDEBUG``.
   
   version
       A version string (e.g. 1.0.0) indicating the required version number
       for the simuPOP module to be loaded. simuPOP will fail to load if the
       installed version is older than the required version.
   
   revision
       Obsolete with the introduction of parameter version.
       
   numThreads
       Number of Threads that will be used to execute a simuPOP script. The
       values can be a positive number (number of threads) or 0 (all available
       cores of the computer, or whatever number set by environmental variable
       ``OMP_NUM_THREADS``). If this parameter is not set, the number of
       threads will be set to 1, or a value set by environmental variable
       ``OMP_NUM_THREADS``.


