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
# Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
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
simuPOP module to load, and how it is loaded, and a class ``simuOpt.simuOpt``
that helps users manage script options.

When simuPOP is loaded, it checkes a few environmental variables
(``SIMUOPTIMIZED``, ``SIMUALLELETYPE``, and ``SIMUDEBUG``) to determine which
simuPOP module to load, and how to load it. More options can be set using the
``simuOpt.setOptions`` function. For example, you can suppress the banner
message when simuPOP is loaded and require a minimal revision of simuPOP for
your script. simuPOP recognize the following commandline arguments

``--optimized``
    Load the optimized version of a simuPOP module.

``--gui=True|False|wxPython|Tkinter``
    Whether or not use a graphical toolkit and which one to use.
    ``--gui=False`` is usually used to run a script in batch mode (do not start
    a parameter input dialog and use interactive user input if a parameter can
    not be determined from command line or a configuraiton file, and it does not
    use its default value (``useDefault`` not set). Please refer to parameter
    ``gui`` for ``simuOpt.setOptions`` for details.

'''

import os, sys

#
# simuOptions that will be checked when simuPOP is loaded. This structure
# can be changed by function setOptions
#

simuOptions = {
    'Optimized': False,
    'AlleleType': 'short',
    'Debug': [],
    'Quiet': False,
    'Revision': None,
    'GUI': True,
}

# Optimized: command line option --optimized or environmental variable SIMUOPTIMIZED
if '--optimized' in sys.argv or os.getenv('SIMUOPTIMIZED') is not None:
    simuOptions['Optimized'] = True

# AlleleType: from environmental variable SIMUALLELETYPE
if os.getenv('SIMUALLELETYPE') in ['short', 'long', 'binary']:
   simuOptions['AlleleType'] = os.getenv('SIMUALLELETYPE')
elif os.getenv('SIMUALLELETYPE') is not None:
    print 'Environmental variable SIMUALLELETYPE can only be short, long, or binary'

# Debug: from environmental variable SIMUDEBUG
if os.getenv('SIMUDEBUG') is not None:
    simuOptions['Debug'].extend(os.getenv('SIMUDEBUG').split(','))

# GUI: from environmental variable SIMUGUI
if os.getenv('SIMUGUI') is not None:
    _gui = os.getenv('SIMUGUI')
elif '--gui' in sys.argv:
    if sys.argv[-1] == '-gui':
        raise exceptions.ValueError('An value is expected for command line option --gui')
    _gui = sys.argv[sys.argv.index('--gui') + 1]
elif True in [x.startswith('--gui=') for x in sys.argv]:
    _gui = sys.argv[[x.startswith('--gui=') for x in sys.argv].index(True)][len('--gui='):]
else:
    _gui = None

if _gui in ['True', 'true', '1']:
    simuOptions['GUI'] = True
elif _gui in ['False', 'false', '0']:
    simuOptions['GUI'] = False
elif _gui == 'wxPython':
    simuOptions['GUI'] = 'wxPython'
elif _gui == 'Tkinter':
    simuOptions['GUI'] = 'Tkinter'
elif _gui is not None:
    print "Invalid value '%s' for environmental variable SIMUGUI or commandline option --gui." % _gui

def setOptions(alleleType=None, optimized=None, gui=None, quiet=None, debug=None, revision=None):
    '''Set options before simuPOP is loaded to control which simuPOP module to
    load, and how the module should be loaded.

    alleleType
        Use the standard, binary or long allele version of the simuPOP module
        if ``alleleType`` is set to 'short', 'binary', or 'long' respectively.
        If this parameter is not set, this function will try to get its value
        from environmental variable ``SIMUALLELETYPE``. The standard (short)
        module will be used if the environmental variable is not defined.

    optimized
        Load the optimized version of a module if this parameter is set to
        ``True`` and the standard version if it is set to ``False``. If this
        parameter is not set (``None``), the optimized version will be used
        if environmental variable ``SIMUOPTIMIZED`` is defined. The standard
        version will be used otherwise.
    
    gui
        Whether or not use graphical user interfaces, and which graphical
        toolkit to use. If this parameter is ``None`` (default), this function
        will check environmental variable ``SIMUGUI`` for a value, and assume
        ``True`` if such an option is unavailable.
        
        gui=True
            allows simuPOP to use ``wxPython``-based dialogs if ``wxPython``
            is available, and use ``Tkinter``-based dialogs if ``Tkinter``
            is available.
        
        gui='Tkinter'
            Use ``Tkinter`` based dialogs even if ``wxPython`` is available.

        gui='wxPython'
            Use ``wxPython``  based dialogs. Usually not needed.
        
        gui=False
            Do not use any graphical toolkit. Run the script in batch mode.

        This option is usually left to ``None`` so that the same script can
        be ran in both GUI and batch mode using commandline option ``--gui``.

    quiet
        If set to ``True``, suppress the banner message when a simuPOP module
        is loaded.

    debug
        A list of debug code (as string) that will be turned on when simuPOP
        is loaded. If this parameter is not set, a list of comma separated
        debug code specified in environmental variable ``SIMUDEBUG``, if
        available, will be used. Note that setting ``debug=[]`` will remove
        any debug code that might have been by variable ``SIMUDEBUG``.

    revision
        A number indicating the required revision number for the simuPOP module
        to be loaded. simuPOP will fail to load if the installed simuPOP is
        older than the required revision. Please check simuPOP ChangeLog for
        the revision number of distributed versions.
    '''
    # Optimized
    if optimized in [True, False]:
        simuOptions['Optimized'] = optimized
    elif optimized is not None:
        raise exceptions.TypeError('Parameter optimized can be either True or False.')
    # Allele type
    if alleleType in ['long', 'binary', 'short']:
        simuOptions['AlleleType'] = alleleType
    elif alleleType is not None:
        raise exceptions.TypeError('Parameter alleleType can be either short, long, or binary.')
    # Graphical toolkit
    if gui in [True, False, 'wxPython', 'Tkinter']:
        simuOptions['GUI'] = gui
    elif gui is not None:
        raise exceptions.TypeError('Parameter gui can be True/False, wxPython or Tkinter.')
    # Quiet
    if quiet in [True, False]:
        simuOptions['Quiet'] = quiet
    elif quiet is not None:
        raise exceptions.TypeError('Parameter quiet can be either True or False.')
    # Debug
    if debug is not None:
        if type(debug) == type(''):
            simuOptions['Debug'] = [debug]
        else:
            simuOptions['Debug'] = debug
    # Revision
    if type(revision) == type(1):
        simuOptions['Revision'] = revision
    elif revision is not None:
        raise exceptions.TypeError('A revision number is expected for parameter revision.')


