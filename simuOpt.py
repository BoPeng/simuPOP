#!/usr/bin/env python

#
# $File: simuOpt.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit http://simupop.sourceforge.net
# for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

'''
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

'''


__all__ = [
    'simuOptions',
    'setOptions'
]

import os, sys, re, time, textwrap
#
# simuOptions that will be checked when simuPOP is loaded. This structure
# can be changed by function setOptions
#

simuOptions = {
    'Optimized': False,
    'AlleleType': 'short',
    'Debug': [],
    'Quiet': False,
    'Version': None,
    'Revision': None,
    'GUI': True,
    'Plotter': None,
    'NumThreads': 1,
}

# Optimized: command line option --optimized or environmental variable SIMUOPTIMIZED
if '--optimized' in sys.argv or os.getenv('SIMUOPTIMIZED') is not None:
    simuOptions['Optimized'] = True

# AlleleType: from environmental variable SIMUALLELETYPE
if os.getenv('SIMUALLELETYPE') in ['short', 'long', 'binary', 'mutant', 'lineage']:
   simuOptions['AlleleType'] = os.getenv('SIMUALLELETYPE')
elif os.getenv('SIMUALLELETYPE') is not None:
    print('Environmental variable SIMUALLELETYPE can only be short, long, binary, mutant, or lineage.')

# Debug: from environmental variable SIMUDEBUG
if os.getenv('SIMUDEBUG') is not None:
    simuOptions['Debug'].extend(os.getenv('SIMUDEBUG').split(','))

# openMP number of threads
if os.getenv('OMP_NUM_THREADS') is not None:
    try:
        simuOptions['NumThreads'] = int(os.getenv('OMP_NUM_THREADS'))
    except:
        print('Ignoring invalid value for environmental variable OMP_NUM_THREADS')

# GUI: from environmental variable SIMUGUI
if os.getenv('SIMUGUI') is not None:
    _gui = os.getenv('SIMUGUI')
elif '--gui' in sys.argv:
    if sys.argv[-1] == '--gui':
        raise ValueError('An value is expected for command line option --gui')
    _gui = sys.argv[sys.argv.index('--gui') + 1]
elif True in [x.startswith('--gui=') for x in sys.argv]:
    _gui = sys.argv[[x.startswith('--gui=') for x in sys.argv].index(True)][len('--gui='):]
else:
    _gui = None

if _gui in ['True', 'true', '1']:
    simuOptions['GUI'] = True
elif _gui in ['False', 'false', '0']:
    simuOptions['GUI'] = False
elif _gui in ['wxPython', 'Tkinter', 'batch', 'interactive']:
    simuOptions['GUI'] = _gui
elif _gui is not None:
    print("Invalid value '%s' for environmental variable SIMUGUI or commandline option --gui." % _gui)

def setOptions(alleleType=None, optimized=None, gui=None, quiet=None,
        debug=None, version=None, revision=None, numThreads=None, plotter=None):
    '''Set options before simuPOP is loaded to control which simuPOP module to
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
    '''
    # if the module has already been imported, check which module
    # was imported
    try:
        _imported = sys.modules['simuPOP'].moduleInfo()
    except Exception as e:
        _imported = {}
    # Allele type
    if alleleType in ['long', 'binary', 'short', 'mutant', 'lineage']:
        # if simuPOP has been imported and re-imported with a different module name
        # the existing module will be used so moduleInfo() will return a different
        # module type from what is specified in simuOptions.
        if _imported and _imported['alleleType'] != alleleType:
            raise ImportError(('simuPOP has already been imported with allele type %s (%s) and cannot be '
                're-imported with allele type %s. Please make sure you import module simuOpt before '
                'any simuPOP module is imported.') % (
                    _imported['alleleType'], ('optimized' if _imported['optimized'] else 'standard'),
                    alleleType))
        simuOptions['AlleleType'] = alleleType
    elif alleleType is not None:
        raise TypeError('Parameter alleleType can be either short, long, binary, mutant or lineage.')
    # Optimized
    if optimized in [True, False]:
        # if simuPOP has been imported and re-imported with a different module name
        # the existing module will be used so moduleInfo() will return a different
        # module type from what is specified in simuOptions.
        if _imported and _imported['optimized'] != optimized:
            raise ImportError(('simuPOP has already been imported with allele type %s (%s) and cannot be '
                're-imported in %s mode. Please make sure you import module simuOpt before '
                'any simuPOP module is imported.') % (
                    _imported['alleleType'], ('optimized' if _imported['optimized'] else 'standard'),
                    'optimized' if optimized else 'standard'))
        simuOptions['Optimized'] = optimized
    elif optimized is not None:
        raise TypeError('Parameter optimized can be either True or False.')        
    # Graphical toolkit
    if gui in [True, False, 'wxPython', 'Tkinter', 'batch']:
        simuOptions['GUI'] = gui
    elif gui is not None:
        raise TypeError('Parameter gui can be True/False, wxPython or Tkinter.')
    # Quiet
    if quiet in [True, False]:
        simuOptions['Quiet'] = quiet
    elif quiet is not None:
        raise TypeError('Parameter quiet can be either True or False.')
    # Debug
    if debug is not None:
        if type(debug) == str:
            simuOptions['Debug'] = [debug]
        else:
            simuOptions['Debug'] = debug
    # Version
    if type(version) == str:
        try:
            major, minor, release = [int(x) for x in re.sub('\D', ' ', version).split()]
        except:
            print('Invalid version string %s' % simuOptions['Version'])
        simuOptions['Version'] = version
    elif version is not None:
        raise TypeError('A version string is expected for parameter version.')
    # Revision
    if type(revision) == int:
        simuOptions['Revision'] = revision
    elif revision is not None:
        raise TypeError('A revision number is expected for parameter revision.')
    # NumThreads
    if type(numThreads) == int:
        simuOptions['NumThreads'] = numThreads
    elif numThreads is not None:
        raise TypeError('An integer number is expected for parameter numThreads.')
    if plotter is not None:
        sys.stderr.write('WARNING: plotter option is deprecated because of the removal of rpy/rpy2 support\n')

