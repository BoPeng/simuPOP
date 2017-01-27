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
    'setOptions',
    'valueNot',
    'valueOr',
    'valueAnd',
    'valueOneOf',
    'valueTrueFalse',
    'valueBetween',
    'valueGT',
    'valueGE',
    'valueLT',
    'valueLE',
    'valueEqual',
    'valueNotEqual',
    'valueIsInteger',
    'valueIsNum',
    'valueIsList',
    'valueValidDir',
    'valueValidFile',
    'valueListOf',
    'valueSumTo',
    'Params',
    'param',
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

#
# define some validataion functions
#
def valueNot(t):
    '''Return a function that returns true if passed option does not equal t,
    or does not passes validator t'''
    def func(val):
        if callable(t):
            return not t(val)
        else:
            return val != t
    return func


def valueOr(t1, t2):
    '''Return a function that returns true if passed option passes validator
    t1 or t2'''
    def func(val):
        if callable(t1) and callable(t2):
            return t1(val) or t2(val)
        else:
            raise ValueError("We expect a function valueXXX")
    return func


def valueAnd(t1, t2):
    '''Return a function that returns true if passed option passes validator
    t1 and t2'''
    def func(val):
        if callable(t1) and callable(t2):
            return t1(val) and t2(val)
        else:
            raise ValueError("We expect a function valueXXX")
    return func


def valueOneOf(*args):
    '''Return a function that returns true if passed option is one of the
    parameters, or one of the values in the only parameter'''
    if len(args) == 1:
        t = args[0]
    else:
        t = args
    if not type(t) in [list, tuple]:
        raise ValueError('argument of valueOneOf should be a list')
    def func(val):
        yes = False
        for item in t:
            if item == val:    # equal value
                return True
            if callable(item): # a test function
                if item(val):
                    return True
        return False
    return func


def valueTrueFalse():
    '''Return a function that returns true if passed option is True or False'''
    return valueOneOf([True, False])


def valueBetween(a,b):
    '''Return a function that returns true if passed option is between value a and b
    (a and b included)
    '''
    def func(val):
        return val >= a and val <=b
    return func


def valueGT(a):
    '''Return a function that returns true if passed option is greater than a'''
    def func(val):
        return val > a
    return func


def valueGE(a):
    '''Return a function that returns true if passed option is greater than or
    equal to a'''
    def func(val):
        return val >= a
    return func


def valueLT(a):
    '''Return a function that returns true if passed option is less than a'''
    def func(val):
        return val < a
    return func


def valueLE(a):
    '''Return a function that returns true if passed option is less than or
    equal to a'''
    def func(val):
        return val <= a
    return func


def valueEqual(a):
    'Return a function that returns true if passed option equals a'
    def func(val):
        return val == a
    return func


def valueNotEqual(a):
    'Return a function that returns true if passed option does not equal a'
    def func(val):
        return val != a
    return func

def valueIsInteger():
    'Return a function that returns true if passed option is an integer (int, long)'
    def func(val):
        return isinstance(val, (int, long))
    return func

def valueIsNum():
    'Return a function that returns true if passed option is a number (int, long or float)'
    def func(val):
        return isinstance(val, (int, long, float))
    return func


def valueIsList(size=None):
    '''Return a function that returns true if passed option is a sequence.
    If a ``size`` is given, the sequence must have the specified size (e.g.
    ``size=3``), or within the range of sizes (e.g. ``size=[1, 5]``). A ``None``
    can be used as unspecified lower or upper bound.'''
    def func(val):
        if not hasattr(val, '__iter__'):
            return False
        if size is not None:
            if isinstance(size, (int, long)):
                return len(val) == size
            elif hasattr(size, '___iter__') and len(size) == 2:
                if size[0] is not None and size[1] is not None:
                    return len(val) >= size[0] and len(val) <= size[1]
                elif size[0] is None:
                    return len(val) <= size[1]
                else:
                    return len(val) >= size[0]
            else:
                raise ValueError('Invalid size specification for simuOpt.valueIsList')
        return True
    return func

def valueSumTo(a, eps=1e-7):
    '''Return a function that returns true if passed option sum up to a.
    An eps value can be specified to allowed for numerical error.'''
    def func(val):
        return abs(sum(val) - a) < eps
    return func

def valueValidDir():
    '''Return a function that returns true if passed option val if a valid
    directory'''
    def func(val):
        '''Params.valueValidDir''' 
        return os.path.isdir(val)
    return func


def valueValidFile():
    '''Return a function that returns true if passed option val if a valid
    file'''
    def func(val):
        '''Params.valueValidFile''' 
        return os.path.isfile(val)
    return func


def valueListOf(t, size=None):
    '''Return a function that returns true if passed option val is a list of
    type ``t`` if ``t`` is a type, if ``v`` is one of ``t`` if ``t`` is a list,
    or if ``v`` passes test ``t`` if ``t`` is a validator (a function). If a
    ``size`` is given, the sequence must have the specified size (e.g.
    ``size=3``), or within the range of sizes (e.g. ``size=[1, 5]``). A ``None``
    can be used as unspecified lower or upper bound.'''
    def func(val):
        if not hasattr(val, '__iter__'):
            return False
        if type(t) in [list, tuple]:
            for i in val:
                if not i in t:
                    return False
        elif callable(t):
            for i in val:
                if not t(i):
                    return False
        else:
            for i in val:
                if type(i) != t:
                    return False
        if size is not None:
            if isinstance(size, (int, long)):
                return len(val) == size
            elif hasattr(size, '___iter__') and len(size) == 2:
                if size[0] is not None and size[1] is not None:
                    return len(val) >= size[0] and len(val) <= size[1]
                elif size[0] is None:
                    return len(val) <= size[1]
                else:
                    return len(val) >= size[0]
            else:
                raise ValueError('Invalid size specification for simuOpt.valueIsList')
        return True
    return func


def _prettyString(value, quoted=False, outer=True):
    '''Return a value in good format, the main purpose is to
      avoid [0.90000001, 0.2].
    '''
    def quote(txt):
        if not quoted:
            return txt
        if not "'" in txt:
            return "'%s'" % txt
        elif not '"' in txt:
            return '"%s"' % txt
        elif not "'''" in txt:
            return "'''%s'''" % txt
        elif not '"""' in txt:
            return '"""%s"""' % txt
        else:
            return "'%s'" % txt.replace("'", "\\'")
    #
    if type(value) in [list, tuple]:
        txt = '[' + ', '.join([_prettyString(x, True, False) for x in value]) + ']'
        if outer:
            return quote(txt)
        else:
            return txt
    elif type(value) == str:
        return quote(value)
    elif outer:
        return quote(str(value))
    else:
        return str(value)


def _prettyDesc(text, indent='', width=80):
    '''Reformat description of options. All lines are joined and all extra
    spaces are removed before the text is wrapped with specified *indent*.
    An exception is that lines with '|' as the first non-space/tab character
    will be outputed as is without the leading '|' symbol**.
    '''
    blocks = []
    for line in text.split('\n'):
        txt = line.strip()
        if txt.startswith('|'):
            blocks.append(txt)
            continue
        txt = re.sub('\s+', ' ', txt)
        if len(blocks) == 0 or blocks[-1].startswith('!'):
            blocks.append(txt)
        else:
            blocks[-1] += ' ' + txt
    txt = []
    for blk in blocks:
        if blk == '|':
            txt.append('')
        elif blk.startswith('|'):
            blk_indent = ''
            for c in blk[1:]:
                if c not in [' ', '\t']:
                    break
                blk_indent += c
            txt.extend(textwrap.wrap(blk[1:], width=width,
                initial_indent=indent + blk_indent, # keep spaces after | character
                subsequent_indent=indent))
        else:
            txt.extend(textwrap.wrap(blk, width=width, initial_indent=indent,
                subsequent_indent=indent))
    return '\n'.join(txt)


def _paramType(opt):
    '''Determine the type of an option'''
    

def _validate(value, opt, options=[]):
    '''validate an option against other options'''
    # if no validator is specified
    if not opt.has_key('validator'):
        return True
    # if a function is given, easy
    if callable(opt['validator']):
        return opt['validator'](value)
    # we need a dictionary
    env = {}
    for o in options:
        if o.has_key('separator'):
            continue
        name = o['name']
        env[name] = o['value']
    env[opt['name']] = value
    return eval(opt['validator'], globals(), env) is True


def _usage(options, msg='', details='', usage='usage: %prog [--opt[=arg]] ...'):
    'Return a usage message.'
    if msg != '':
        message = _prettyDesc(msg, width=80) + '\n\n' + _prettyDesc(details, width=80) + '\n'
    else:
        message = _prettyDesc(details, width=80) + '\n'
    #
    message += '''\n%s

options:
  -h, --help
        Display this help message and exit.

  --config=ARG (default: None)
        Load parameters from a configuration file ARG.

  --optimized
        Run the script using an optimized simuPOP module.

  --gui=[batch|interactive|True|Tkinter|wxPython] (default: None)
        Run the script in batch, interactive or GUI mode.

''' % usage.replace('%prog', os.path.basename(sys.argv[0]))
    for opt in options:
        if opt.has_key('separator'):
            continue
        if opt.has_key('separator'):
            continue
        name = ''
        # this is obsolete
        if opt.has_key('arg'):
            if opt['gui_type'] != 'boolean':
                name += '-%s=ARG, ' % opt['arg']
            else:
                name += '-%s, ' % opt['arg']
        #
        if opt['gui_type'] != 'boolean':
            name += '--%s=ARG' % opt['name']
        else:
            name += '--%s' % opt['name']
        #
        if isinstance(opt['default'], str):
            defaultVal = "'%s'" % opt['default']
        elif opt['default'] is None:
            defaultVal = 'None'
        else:
            defaultVal = _prettyString(opt['default'])
        #
        message += '  %s  (default: %s)\n' % (name, defaultVal)
        if opt.has_key('description'):
            message += _prettyDesc(opt['description'], indent=' '*8) + '\n'
        elif opt.has_key('label'):
            message += _prettyDesc(opt['label'], indent=' '*8) + '\n'
        message += '\n'
    return message.strip()

def _getParamValue(p, val, options):
    ''' try to get a value from value, raise exception if error happens. '''
    if p['gui_type'] == 'separator':
        raise ValueError('Cannot get a value for separator')
    # if we are giving a unicode string, convert!
    if type(val) == unicode:
        val = str(val)
    # remove quotes from string?
    if p.has_key('allowedTypes') and str in p['allowedTypes'] and str == type(val):
        for quote in ['"', "'", '"""', "'''"]:
            if val.startswith(quote) and val.endswith(quote):
                val = val[len(quote):-len(quote)]
                break
    if (not p.has_key('allowedTypes')) or type(val) in p['allowedTypes']:
        if not _validate(val, p, options):
            raise ValueError("Value '%s' is not allowed for parameter %s" % \
                (str(val), p['name']))
        return val
    # handle another 'auto-boolean' case
    elif p['gui_type'] == 'boolean':
        if val in ['1', 'true', 'True']:
            return True
        elif val in ['0', 'false', 'False']:
            return False
        else:
            raise ValueError('Expect 0/1, true/false for boolean values for parameter %s ' % p['name'])
    # other wise, need conversion
    if type(val) in [str, unicode]:
        try:
            val = eval(val)
        except:
            # may be we have a list of string?
            items = val.split(',')
            if len(items) > 1: # is actually a list
                val = []
                for i in items:
                    val.append(i.strip())
    # evaluated type is OK now.
    if type(val) in p['allowedTypes']:
        if not _validate(val, p, options):
            raise ValueError("Default value '" + str(val) + "' for option '" + p['name'] + "' does not pass validation.")
        return val
    elif list in p['allowedTypes'] or tuple in p['allowedTypes']:
        if not _validate([val], p, options):
            raise ValueError("Value "+str([val])+' does not pass validation')
        return [val]
    elif type(val) == bool and int in p['allowedTypes']: # compatibility problem
        return val
    elif type(val) == unicode and str in p['allowedTypes']:
        return str(val)
    else:
        raise ValueError('Type of input parameter "' + str(val) + '" is disallowed for option ' +
            p['name'])
    #print(p, val)

class _paramDialog:
    def __init__(self, options, title = '', description='', details='', nCol=1):
        #
        # now, initialize variables
        self.options = options
        self.title = title
        if nCol is None:
            self.nCol = len(self.options)//20 + 1
        else:
            self.nCol = nCol
        self.description = _prettyDesc(description, indent='', width=55*self.nCol)
        self.details = details

    def setLayout(self, useTk=False):
        '''Design a layout for the parameter dialog. Currently only used by wxPython'''
        row = 0
        for opt in self.options:
            if not (opt.has_key('label') or opt.has_key('separator')):
                continue
            if opt.has_key('separator'):
                row += 2
                continue
            elif opt.has_key('label'):
                row += 1
            if opt.has_key('chooseFrom'):
                if len(opt['chooseFrom']) <= 3 or useTk:
                    rspan = len(opt['chooseFrom'])
                else:
                    rspan = len(opt['chooseFrom']) * 4 / 5
                row += rspan - 1
            elif opt.has_key('chooseOneOf') and useTk:
                rspan = len(opt['chooseOneOf'])
                row += rspan - 1
        nRow = row
        if nRow // self.nCol * self.nCol == nRow:
            nRow = nRow // self.nCol
        else:
            nRow = nRow // self.nCol + 1
        # set row and col for each item.
        headerRow = 0
        r = 0
        c = 0
        for opt in self.options:
            if not (opt.has_key('label') or opt.has_key('separator')):
                continue
            if opt.has_key('separator'):
                # start of a column
                if r == 0:
                    opt['layout'] = r, c, 1
                    r += 1
                else:
                    r += 1
                    # middle, the last one?
                    if r == nRow - 1:
                        if c + 1 == self.nCol:
                            opt['layout'] = r - 1, c, 1
                        else:
                            c += 1
                            opt['layout'] = 0, c, 1
                            r = 1
                    else:
                        opt['layout'] = r, c, 1
                        r += 1
                continue
            if opt.has_key('chooseFrom'):
                if len(opt['chooseFrom']) <= 3 or useTk:
                    rspan = len(opt['chooseFrom'])
                else:
                    rspan = len(opt['chooseFrom']) * 4 / 5
                opt['layout'] = r, c, rspan
                r += rspan
                if r >= nRow and c + 1 < self.nCol:
                    nRow = r
                    r = 0
                    c += 1
            elif opt.has_key('chooseOneOf') and useTk:
                rspan = len(opt['chooseOneOf'])
                opt['layout'] = r, c, rspan
                r += rspan
                if r >= nRow and c + 1 < self.nCol:
                    nRow = r
                    r = 0
                    c += 1
            else:
                opt['layout'] = r, c, 1
                r += 1
                if r >= nRow and c + 1 < self.nCol:
                    nRow = r
                    r = 0
                    c += 1
        return nRow, self.nCol

    def createDialog(self):
        raise SystemError('Please define createDialog')

    def runDialog(self):
        raise SystemError('Please define runDialog')

    def getParam(self):
        '''Create, run a dialog, return result'''
        self.createDialog()
        self.runDialog()
        #
        # after the run has completed
        return not self.cancelled


class _tkParamDialog(_paramDialog):
    def __init__(self, options, title = '', description='', details='', nCol=1):
        ''' get options from a given options structure '''
        _paramDialog.__init__(self, options, title, description, details, nCol)
        # style wise, this seems to be very bad
        import Tkinter
        import tkFont
        if sys.version_info[0] < 3:
            globals()['Tkinter'] = Tkinter
            globals()['tkFont'] = tkFont
        else:
            globals()['tkinter'] = tkinter

    def denyWindowManagerClose(self):
        '''Don't allow WindowManager close'''
        x = Tkinter.Tk()
        x.withdraw()
        x.bell()
        x.destroy()

    def onHelpOK(self, event):
        self.helpDlg.destroy()

    def onOpen(self, event):
        widget = event.widget
        opt = self.options[self.entryWidgets.index(widget)]
        import tkFileDialog as fileDlg
        if opt['gui_type'] == 'browseFile':
            filename = fileDlg.askopenfilename(title=opt['label'])
            if filename is not None:
                # only available in Python 2.6
                widget.delete(0, Tkinter.END)
                if 'relpath' in dir(os.path):
                    widget.insert(0, os.path.relpath(filename))
                else:
                    widget.insert(0, filename)
        else:
            dirname = fileDlg.askdirectory(title=opt['label'])
            if dirname is not None:
                widget.delete(0, Tkinter.END)
                if 'relpath' in dir(os.path):
                    widget.insert(0, os.path.relpath(dirname))
                else:
                    widget.insert(0, dirname)

    def createHelpDialog(self):
        self.helpDlg = Tkinter.Toplevel(self.app)
        self.helpDlg.title('Help for ' + self.title)
        #
        msg = Tkinter.Text(self.helpDlg, wrap=Tkinter.WORD)
        msg.insert(Tkinter.END, _usage(self.options, '', self.details))
        msg.grid(column=0, row=0, pady=10, padx=10,
            sticky = Tkinter.E + Tkinter.W + Tkinter.N + Tkinter.S)
        # scrollbar
        sb = Tkinter.Scrollbar(self.helpDlg)
        sb.config(command=msg.yview)
        sb.grid(column=1, row=0, sticky=Tkinter.N + Tkinter.S,
            padx=0, pady=10)
        msg["yscrollcommand"] = sb.set
        self.helpDlg.columnconfigure(0, weight=1)
        self.helpDlg.rowconfigure(0, weight=1)
        self.helpDlg.rowconfigure(1, weight=0)
        okButton = Tkinter.Button(self.helpDlg, takefocus=1, text="OK")
        okButton.grid(column=0, row=1, columnspan=2, pady=10, padx=10)
        # bind the keyboard events to the widget
        okButton.bind("<Return>", self.onHelpOK)
        okButton.bind("<Button-1>", self.onHelpOK)
        self.helpDlg.bind("<Escape>", self.onHelpOK)

    def onHelp(self, event):
        self.createHelpDialog()
        self.app.wait_window(self.helpDlg)

    def onCancel(self, event):
        '''When ESC is pressed cancel'''
        self.cancelled = True
        self.app.quit()

    def onOK(self, event):
        ''' get result and convert values '''
        for g in range(len(self.entryWidgets)):
            if self.entryWidgets[g] == None:
                # for the case of non-GUI default parameter that do not pass validation (e.g. 5 for [5])
                if self.options[g].has_key('value'):
                    self.options[g]['value'] = _getParamValue(self.options[g], self.options[g]['value'], self.options)
                continue
            try:
                if self.options[g]['gui_type'] == 'boolean':
                    # gets 0/1 for false/true
                    var = self.options[g]['value'].get()
                    val = _getParamValue(self.options[g], var == 1, self.options)
                elif self.options[g]['gui_type'] == 'chooseOneOf':
                    sel = self.entryWidgets[g].curselection()
                    items = self.options[g]['chooseOneOf'][int(sel[0])]
                    val = _getParamValue(self.options[g], items, self.options)
                elif self.options[g]['gui_type'] == 'chooseFrom':
                    sel = self.entryWidgets[g].curselection()
                    items = [self.options[g]['chooseFrom'][int(x)] for x in sel]
                    val = _getParamValue(self.options[g], items, self.options)
                else:
                    val = _getParamValue(self.options[g], self.entryWidgets[g].get(), self.options)
            except Exception as e:
                for lab in self.labelWidgets:
                    if lab is not None:
                        lab.configure(fg='black')
                # set this one to red
                print('Error handling paramter %s: %s' % (self.options[g]['name'], e))
                self.labelWidgets[g].configure(fg='red')
                self.entryWidgets[g].focus_force()
                return
            else:
                # convert to values
                self.options[g]['value'] = val
        # get all results and return
        self.cancelled = False
        self.app.quit()

    def createDialog(self):
        self.app = Tkinter.Tk()
        self.app.protocol('WM_DELETE_WINDOW', self.denyWindowManagerClose)
        self.app.title(self.title)
        #
        # the main window
        self.entryWidgets = [None]*len(self.options)
        self.labelWidgets = [None]*len(self.options)
        # all use grid management
        # top message
        topMsg = Tkinter.Label(self.app, text=self.description.strip(), justify=Tkinter.LEFT)
        topMsg.grid(row=0, column=0, columnspan = 2 * self.nCol, sticky=Tkinter.W,
            padx=10, pady=10)
        # find out number of items etc
        numRows, numCols = self.setLayout(True)
        # all entries
        for g,opt in enumerate(self.options):
            if not (opt.has_key('label') or opt.has_key('separator')):
                continue
            r, c, rspan = opt['layout']
            # skip the top label...
            r += 1
            # use different entry method for different types
            if opt['gui_type'] == 'separator':
                self.labelWidgets[g] = Tkinter.Label(self.app, text=opt['separator'])
                f = tkFont.Font(font=self.labelWidgets[g]["font"]).copy()
                f.config(weight='bold')
                self.labelWidgets[g].config(font=f)
                self.labelWidgets[g].grid(column=c*2, row= r,
                    columnspan=2, ipadx=0, padx=10, sticky=Tkinter.W + Tkinter.N + Tkinter.S)
                self.entryWidgets[g] = None
                continue
            value = self.options[g]['value']
            if value is None:
                value = self.options[g]['default']
            if opt['gui_type'] == 'chooseOneOf':    # single choice
                self.labelWidgets[g] = Tkinter.Label(self.app, text=opt['label'])
                self.labelWidgets[g].grid(column=c*2, row= r,
                    padx=10, rowspan = 1, sticky=Tkinter.W, pady=2)
                self.entryWidgets[g] = Tkinter.Listbox(self.app, selectmode=Tkinter.SINGLE,
                    exportselection=0, height = rspan)
                self.entryWidgets[g].grid(column=c*2+1, row=r,
                    padx=10, rowspan = rspan, pady=2)
                for entry in opt['chooseOneOf']:
                    self.entryWidgets[g].insert(Tkinter.END, str(entry))
                if value is not None:
                    self.entryWidgets[g].select_set(opt['chooseOneOf'].index(value))
            elif opt['gui_type'] == 'chooseFrom':    # multiple choice
                self.labelWidgets[g] = Tkinter.Label(self.app, text=opt['label'])
                self.labelWidgets[g].grid(column=c*2, row=r,
                    padx=10, sticky=Tkinter.W, pady=2)
                self.entryWidgets[g] = Tkinter.Listbox(self.app, selectmode=Tkinter.EXTENDED,
                    exportselection=0, height = rspan)
                self.entryWidgets[g].grid(column=c*2+1, row=r,
                    padx=10, rowspan = rspan, pady=2)
                for entry in opt['chooseFrom']:
                    self.entryWidgets[g].insert(Tkinter.END, str(entry))
                if value is not None:
                    if type(value) in [tuple, list]:
                        for val in value:
                            self.entryWidgets[g].select_set( opt['chooseFrom'].index(val))
                    else:
                        self.entryWidgets[g].select_set( opt['chooseFrom'].index( value ))
            elif opt['gui_type'] == 'boolean':
                iv = Tkinter.IntVar()
                iv.set(self.options[g]['value'] == True) # value can be None, True or False
                self.options[g]['value'] = iv
                self.entryWidgets[g] = Tkinter.Checkbutton(self.app, height=1,
                    text = opt['label'], variable=self.options[g]['value'])
                self.entryWidgets[g].grid(column=c*2, row=r, padx=10, columnspan=2, 
                    sticky=Tkinter.W)
            else:
                self.labelWidgets[g] = Tkinter.Label(self.app, text=opt['label'])
                self.labelWidgets[g].grid(column=c*2, row=r,
                    padx=10, sticky=Tkinter.W, pady=2)
                self.entryWidgets[g] = Tkinter.Entry(self.app)
                self.entryWidgets[g].grid(column=c*2+1, row=r,
                    padx=10, ipadx=0, pady=2)
                if opt['gui_type'] in ['browseFile', 'browseDir']:
                    self.entryWidgets[g].bind('<Double-Button-1>', self.onOpen)
                # put default value into the entryWidget
                self.entryWidgets[g].insert(0, _prettyString(value))
            self.entryWidgets[g].bind("<Return>", self.onOK)
            self.entryWidgets[g].bind("<Escape>", self.onCancel)
        # help button: left
        helpButton = Tkinter.Button(self.app, takefocus=1, text="Help")
        helpButton.bind("<Return>", self.onHelp)
        helpButton.bind("<Button-1>", self.onHelp)
        helpButton.grid(column=0, columnspan=self.nCol, row = numRows+1, pady=20, sticky='w', padx=20)
        # ok button: right
        okButton = Tkinter.Button(self.app, takefocus=1, text="Run!")
        okButton.bind("<Return>", self.onOK)
        okButton.bind("<Button-1>", self.onOK)
        okButton.grid(column=self.nCol, columnspan=self.nCol, row = numRows+1, sticky='e', pady=20, padx=10)
        # cancel button: middle
        cancelButton = Tkinter.Button(self.app, takefocus=1, text="Cancel")
        cancelButton.bind("<Return>", self.onCancel)
        cancelButton.bind("<Button-1>", self.onCancel)
        cancelButton.grid(column=0, columnspan=2*self.nCol, row = numRows+1, pady=20)
        #
        self.app.bind("<Escape>", self.onCancel)
        # first un-none
        for g in range(len(self.options)):
            if self.entryWidgets[g] is not None:
                self.entryWidgets[g].focus_force()
                break

    def runDialog(self):
        self.app.mainloop()    # run it!
        self.app.destroy()    # button_click didn't destroy self.app, so we do it now


# get options from a given options structure
class _wxParamDialog(_paramDialog):
    def __init__(self, options, title = '', description='', details='', nCol=1):
        _paramDialog.__init__(self, options, title, description, details, nCol)
        import wx
        import wx.lib.filebrowsebutton
        globals()['wx'] = wx

    def onHelp(self, event):
        # open another dialog
        helpDlg = wx.Dialog(parent=self.dlg, id=-1,
                    title='Help for ' + self.title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(wx.TextCtrl(parent=helpDlg, id=-1, size=[450, 400],
            style=wx.TE_MULTILINE | wx.TE_READONLY,
            value=_usage(self.options, '', self.details)), 0, wx.ALL, 10)
        self.addButton(wx.ID_OK, "", lambda event:helpDlg.EndModal(wx.ID_OK), 
                    helpDlg, box)
        helpDlg.SetSizerAndFit(box)
        helpDlg.Layout()
        helpDlg.ShowModal()
        helpDlg.Destroy()

    def onOK(self, event):
        ''' get result and convert values '''
        for g in range(len(self.entryWidgets)):
            if self.entryWidgets[g] == None:
                # for the case of non-GUI default parameter that do not pass validation (e.g. 5 for [5])
                if self.options[g].has_key('value'):
                    self.options[g]['value'] = _getParamValue(self.options[g], self.options[g]['value'], self.options)
                continue
            try:
                if self.options[g]['gui_type'] == 'chooseOneOf':
                    val = _getParamValue(self.options[g],
                            self.options[g]['chooseOneOf'][int(self.entryWidgets[g].GetSelection())],
                            self.options)
                elif self.options[g]['gui_type'] == 'chooseFrom':
                    items = []
                    for s in range(len(self.options[g]['chooseFrom'])):
                        if self.entryWidgets[g].IsChecked(s):
                            items.append(self.options[g]['chooseFrom'][s])
                    val = _getParamValue(self.options[g], items, self.options)
                else:
                    val = _getParamValue(self.options[g], self.entryWidgets[g].GetValue(), self.options)
            except Exception as e:
                # incorrect value
                # set to red
                # clear other red colors
                for lab in self.labelWidgets:
                    if lab is not None:
                        lab.SetForegroundColour('black')
                        lab.Refresh()
                print('Error handling paramter %s' % self.options[g]['name'])
                print(e)
                if self.labelWidgets[g] is not None:
                    # happen to a boolean set...
                    # set this one to red
                    self.labelWidgets[g].SetForegroundColour('red')
                    self.entryWidgets[g].SetFocus()
                    self.labelWidgets[g].Refresh()
                return
            else:
                # convert to values
                self.options[g]['value'] = val
        #
        self.cancelled = False
        self.dlg.EndModal(wx.ID_OK)

    def onCancel(self, event):
        '''When ESC is pressed cancel, clear values and return'''
        self.cancelled = True
        self.dlg.EndModal(wx.ID_CANCEL)

    def addButton(self, ID, text, func, parent, box):
        button = wx.Button(parent, ID, text)
        box.Add(button, 0, wx.ALIGN_CENTER | wx.ALL, border=5)
        parent.Bind(wx.EVT_BUTTON, func, button)

    def createDialog(self):
        self.app = wx.PySimpleApp(0)
        self.dlg = wx.Dialog(parent=None, id=-1, title = self.title)
        self.entryWidgets = [None]*len(self.options)
        self.labelWidgets = [None]*len(self.options)
        #
        # the main window
        box = wx.BoxSizer(wx.VERTICAL)
        # do not use a very long description please
        topLabel = wx.StaticText(parent=self.dlg, id=-1, label=self.description.strip())
        box.Add(topLabel, 0, wx.EXPAND | wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, 15)
        # add a box for all ...
        paraBox = wx.FlexGridSizer(cols = self.nCol)
        for col in range(self.nCol):
            paraBox.AddGrowableCol(col)
        #
        # add several GridBagSizer
        gridBox = []
        for col in range(self.nCol):
            gridBox.append(wx.GridBagSizer(vgap=2, hgap=5))
            #gridBox[-1].AddGrowableCol(0)
            gridBox[-1].AddGrowableCol(1)
            paraBox.Add(gridBox[-1], 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 10)
        box.Add(paraBox, 1, wx.EXPAND | wx.ALL, 5)
        # count numbers of valid parameters..
        # chooseFrom count as many
        numRows, numCols = self.setLayout(False)
        # all entries
        for g,opt in enumerate(self.options):
            if opt['gui_type'] == 'hidden':
                continue
            r, c, rspan = opt['layout']
            if opt['gui_type'] == 'separator':
                self.labelWidgets[g] = wx.StaticText(parent=self.dlg, id=-1, label=opt['separator'])
                f = self.labelWidgets[g].GetFont()
                f.SetWeight(wx.BOLD)
                self.labelWidgets[g].SetFont(f)
                gridBox[c].Add(self.labelWidgets[g], (r, 0), span=(1, 2),
                        border=2, flag=wx.ALIGN_LEFT | wx.ALIGN_BOTTOM | wx.BOTTOM)
                # no entry widget
                self.entryWidgets[g] = None
                continue
            value = self.options[g]['value']
            if value is None:
                value = self.options[g]['default']
            # label
            # use different entry method for different types
            if opt.has_key('description'):
                tooltip = _prettyDesc(opt['description'])
            else:
                tooltip = 'arg: ' + opt['name']
            if opt['gui_type'] == 'chooseOneOf':    # single choice
                self.labelWidgets[g] = wx.StaticText(parent=self.dlg, id=-1, label=opt['label'])
                gridBox[c].Add(self.labelWidgets[g], (r, 0), flag=wx.ALIGN_LEFT
                    | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, border=2)
                self.entryWidgets[g] = wx.Choice(parent=self.dlg, id=g,
                    choices = [str(x) for x in opt['chooseOneOf']])
                gridBox[c].Add(self.entryWidgets[g], (r, 1), flag=wx.EXPAND
                    | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, border=2)
                # if an value is given through command line argument or configuration file
                if value is not None:
                    try:
                        self.entryWidgets[g].SetSelection(opt['chooseOneOf'].index(value))
                    except:
                        raise ValueError('Value: %s is not one of %s.' % (str(value), str(opt['chooseOneOf'])))
            elif opt['gui_type'] == 'chooseFrom':    # multiple choice
                self.labelWidgets[g] = wx.StaticText(parent=self.dlg, id=-1, label=opt['label'])
                gridBox[c].Add(self.labelWidgets[g], (r, 0), flag=wx.ALIGN_LEFT
                    | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, border=2)
                # the height is a little bit too much...
                self.entryWidgets[g] = wx.CheckListBox(parent=self.dlg, id=g,
                    choices = [str(x) for x in opt['chooseFrom']])
                if value is not None:
                    if type(value) in [list, tuple]:
                        for val in value:
                            self.entryWidgets[g].Check(opt['chooseFrom'].index(val))
                    else:
                        self.entryWidgets[g].Check(opt['chooseFrom'].index(value))
                gridBox[c].Add(self.entryWidgets[g], (r, 1), span=(rspan, 1), flag=wx.EXPAND
                    | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, border=2)
            elif opt['gui_type'] == 'boolean':
                self.entryWidgets[g] = wx.CheckBox(parent=self.dlg, id=g, label = opt['label'])
                if value is not None:
                    self.entryWidgets[g].SetValue(value)
                gridBox[c].Add(self.entryWidgets[g], (r, 0), span=(1, 2), flag=wx.EXPAND
                    | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, border=2)
            elif opt['gui_type'] == 'browseFile':
                self.entryWidgets[g] = wx.lib.filebrowsebutton.FileBrowseButton(parent=self.dlg, id=g,
                    labelText=opt['label'], initialValue=value)
                gridBox[c].Add(self.entryWidgets[g], (r, 0), span=(1,2), flag=wx.EXPAND | wx.ALIGN_LEFT
                    | wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, border=2)
            elif opt['gui_type'] == 'browseDir':
                self.entryWidgets[g] = wx.lib.filebrowsebutton.DirBrowseButton(parent=self.dlg, id=g,
                    labelText=opt['label'], startDirectory=value)
                self.entryWidgets[g].SetValue(value)
                gridBox[c].Add(self.entryWidgets[g], (r, 0), span=(1, 2), flag=wx.EXPAND | wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL
                    | wx.BOTTOM | wx.TOP, border=2)
            else: # an edit box
                self.labelWidgets[g] = wx.StaticText(parent=self.dlg, id=-1, label=opt['label'])
                gridBox[c].Add(self.labelWidgets[g], (r, 0), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL
                    | wx.BOTTOM | wx.TOP, border=2)
                txt = _prettyString(value)
                self.entryWidgets[g] = wx.TextCtrl(parent=self.dlg, id=g, value=txt)
                gridBox[c].Add(self.entryWidgets[g], (r, 1), flag=wx.EXPAND | wx.BOTTOM | wx.TOP | wx.ALIGN_CENTER_VERTICAL, border=2)
            self.entryWidgets[g].SetToolTipString(tooltip)
        # help button
        buttonBox = wx.GridSizer(cols=3)
        self.addButton(wx.ID_HELP, '', self.onHelp, self.dlg, buttonBox)
        self.addButton(wx.ID_CANCEL, '', self.onCancel, self.dlg, buttonBox)
        self.addButton(wx.ID_OK, '', self.onOK, self.dlg, buttonBox)
        box.Add(buttonBox, 0, wx.ALL | wx.EXPAND, 5)
        self.dlg.SetSizerAndFit(box)
        self.dlg.Layout()
        # set focus to the first un-none entry
        for g in range(len(self.options)):
            if self.entryWidgets[g] is not None:
                self.entryWidgets[g].SetFocus()
                break

    def runDialog(self):
        self.dlg.ShowModal()
        self.dlg.Destroy()
        # This app may get in the way when getParam() or another wxPython
        # related function is called.
        del self.app

def param(name='', default=None, **kwargs):
    '''A simple wrapper that allows the specification of a parameter using a
    function instead of a dictionary. Please refer to class
    ``simuOpt.Params`` for allowed keyword arguments and their meanings.'''
    return kwargs.update({'name':name, 'default':default})

class Params:
    '''
    class Params provides a uniform interface for simuPOP scripts to handle
    parameters. It allows users to get parameters from command line options,
    a configuration file, a parameter input dialog (*tkInter* or *wxPython*) or
    from interative input. This class provides parameter validation, conversion
    and and some utility functions to print, save and restore parameters.

    A Params object accepts a parameter specification list that consists of
    dictionaries with pre-defined keys. Each item defines an option in terms of
    command line option, entry name in a configuration file, label in a
    parameter input dialog, acceptable types, validation rules and a default
    value. The following keys are currently supported:

    name
        Long command line option name.  For example ``'version'``  checks the
        presence of argument ``--version``. For example, ``'mu'`` matches
        command line option ``--mu=0.001`` or ``--mu 0.001``. **This item
        defines the name of an option and cannot be ignored**. An options that
        does not expect a value is identified as a single BooleanType in
        allowedTypes or a default value ``False`` when no allowedTypes is
        defined. Such a value should have default value ``False`` and the
        presence of this argument in the command line (e.g. ``--verbose``)
        change it to ``True``.

    label
        The label of the input field in a parameter input dialog. It will also
        be used as the prompt for this option during interactive parameter
        input. **Options without a label will not be displayed in the parameter
        input dialog and will not be saved to a configuration file**. A typical
        example of such an option is ``--version``.

    default
        Default value for this parameter. It is used as the default value in
        the parameter input dialog, and as the option value when a user presses
        ``Enter`` directly during interactive parameter input. **A default
        value is required for all options**.

    description
        A long description of this parameter. This description will be put into
        the usage information, and as parameter tooltip in the parameter input
        dialog. This string will be reformatted when it is written to a usage
        string (remove newlines and extra spaces and re-indent), with the
        exception that **lines with '|' as the first non-space/tab character
        will be outputed as is without the leading '|' symbol**.

    allowedTypes
        A list of acceptable types of this option. class ``Params`` will try
        to convert user input to these types. For example, if ``allowedTypes``
        is ``list`` or ``tuple`` and the user's input is a
        scalar, the input will be converted to a list automatically. An option
        will not be accepted if such conversion fails. If this item is not
        specified, the type of the default value will be used. If only one type
        is acceptable, a single value can be used as input (ignore []).

    validator
        An expression or a function to validate the parameter. If an expression
        (a string) is used, it will be evaluated using current values of
        parameters as inputs. If a function is specified, it will be called
        with the value of the parameter. The option will not be accepted if
        the expression is evalulated as ''False'' or if the function returns
        ``False``. This module defines a large number of such validation
        functions but user defined functions are also acceptable.

    type
        Type of the input parameter. Which can be a datatype (e.g. ``bool``),
        a list of acceptable types (e.g. ``(int, long)``), or one of 
        ``'chooseOneOf', values``, ``'chooseFrom', values``, ``'filename'``
        (a valid filename), ``'dirname'`` (a valid directory name),
        ``'integer'`` (``int`` or ``long``), ``'integers'`` (list of integers),
        ``'number'`` (``int``, ``long`` or ``float``), ``'numbers'`` (list of
        numbers), ``'string'`` (the same as ``str``), ``'strings'`` (list of
        strings). These types will advise class ``simuOpt.Params`` how to
        accept parameters, and how to convert user input to appropriate
        types. For example, a file browser will be used to browse for a valid
        filename in a parameter input dialog (if ``--gui=True``) for a
        parameter of type ``'filename'``, and input ``'52'`` will be
        automatically converted to ``(52,)`` for a parameter of type
        ``'integers'``.

    separator
        This item specifies a separator (group header) in the parameter input
        dialog. All other fields are ignored.
   
    arg, longarg, useDefault, chooseFrom, chooseOneOf, allowedTypes
        These parameters are deprecated because of the introduction of the
        'name', 'gui_type' keys and the 'batch' mode.

    Not all keys need to be specified in each option description. Missing
    values are handled using some internal rules. For example, items without
    a ``label`` will not be displayed on the parameter dialog. This will
    effectively *hide* a parameter although users who know this parameter
    can set it using command line options.

    The ``Params.Params`` class defines a number of functions to collect,
    validate, and manipulate parameters using this parameter specification
    list.

    As a shortcut to create a Params object with a number of attributes,
    a Params object can be created with additional ``key=value`` pairs that
    could be assessed as attributes. This is used to create a ``Params``
    object in which *parameters* are assigned directly.
    '''
    def __init__(self, options=[], doc='', details='', **kwargs):
        '''Create a ``Params`` oject using a list of parameter specification
        dictionaries *options*. Additional *doc* and *details* can be
        specified which will be displayed as script summary (on the top of
        a parameter input dialog) and script introduction (the first part of
        a help message), respectively. Additional attributes could be assigned
        to a ``Params`` object as keyword arguments. Note that it is customary
        to use module document (the first string object in a Python script) as
        *details*, using parameter ``details=__doc__``.
        '''
        #
        # validator
        if type(options) != list:
            raise ValueError('An option specification list is expected')
        #
        self.options = []
        self.dict = {}
        for opt in options:
            if type(opt) != dict:
                raise ValueError('An option specification dictionary is expected')
            self.addOption(**opt)
        #
        if type(doc) != str:
            raise ValueError('Parameter doc must be a string.')
        if type(details) != str:
            raise ValueError('Parameter details must be a string.')
        self.doc = doc.strip()
        self.details = details.strip()
        self.processedArgs = []
        # allow the change of default parameter or addition of parameters
        # using additional key=value pairs.
        for key in kwargs.keys():
            if self.dict.has_key(key):
                self.dict[key]['value'] = kwargs[key]
            else:
                self.__dict__[key] = kwargs[key]

    def __getattr__(self, name):
        'Return the value of a parameter as an attribute.'
        if self.dict.has_key(name):
            return self.dict[name]['value']
        raise AttributeError('Can not locate attribute %s.' % name)

    def __setattr__(self, name, value):
        'Set the value of a parameter as an attribute.'
        if self.__dict__.has_key('dict') and self.dict.has_key(name):
            self.dict[name]['value'] = value
        else:
            self.__dict__[name] = value

    def addOption(self, name='', default=None, **kwargs):
        '''
        Append an entry to the parameter specification list. Dictionary
        entries should be specified as keyword arguments such as
        ``name='option'``. More specifically, you can specify parameters
        ``name`` (required), ``label``, ``default`` (required),
        ``description``, ``validator``, ``type``, and ``separator``.
        This option will have a name specified by ``name`` and an initial
        default value specified by ``default``.
        '''
        allowed_keys = ['name', 'default', 'label', 'description', 'validator',
            'separator', 'type', 'arg', 'longarg', 'allowedTypes', 'useDefault',
            'validate', 'chooseOneOf', 'chooseFrom']
        #
        methods = ['asDict', 'asList', 'getParam', 'loadConfig', 'saveConfig',
            'usage', 'processArgs', 'guiGetParam', 'termGetParam', 'addOption']
        #
        reserved_options = ['optimized', 'gui', 'config', 'help']
        #
        opt = {}
        # allow the input of single value for allowed types.
        for key in kwargs:
            if key not in allowed_keys:
                raise ValueError('Invalid option specification key %s' % key)
            # I used not hasattr(kwargs[key], '__iter__') to test single element but
            # hasattr(tuple, '__iter__') returns True. Using isinstance
            # solves this problem.
            if key == 'allowedTypes':
                if 'DBG_COMPATIBILITY' in simuOptions['Debug']:
                    sys.stderr.write('WARNING: allowedTypes is obsolete and might be removed from a future version of simuPOP.\n')
                if isinstance(kwargs[key], type(int)):
                    opt[key] = [kwargs[key]]
                else:
                    # must be a list.
                    opt[key] = kwargs[key]
            elif key == 'validate':
                if 'DBG_COMPATIBILITY' in simuOptions['Debug']:
                    sys.stderr.write('WARNING: key validate has been renamed to validator.\n')
                opt['validator'] = kwargs[key]
            elif key == 'longarg':
                if 'DBG_COMPATIBILITY' in simuOptions['Debug']:
                    sys.stderr.write('WARNING: longarg is obsolete and might be removed from a future version of simuPOP.\n')
                opt['name'] = kwargs[key].rstrip('=')
                if not kwargs[key].endswith('='):
                    opt['gui_type'] = 'boolean'
            elif key == 'arg':
                if 'DBG_COMPATIBILITY' in simuOptions['Debug']:
                    sys.stderr.write('WARNING: arg is obsolete and might be removed from a future version of simuPOP.\n')
                opt[key] = kwargs[key].rstrip(':')
                if not kwargs[key].endswith(':'):
                    opt['gui_type'] = 'boolean'
                if len(opt['arg']) != 1 or not opt['arg'][0].isalpha():
                    raise ValueError('Short arg should have one and only one alphabetic character.')
            else:
                if key == 'chooseOneOf':
                    if 'DBG_COMPATIBILITY' in simuOptions['Debug']:
                        sys.stderr.write('WARNING: chooseOneOf is obsolete and might be removed from a future version of simuPOP. Use type instead.\n')
                elif key == 'chooseFrom':
                    if 'DBG_COMPATIBILITY' in simuOptions['Debug']:
                        sys.stderr.write('WARNING: chooseFrom is obsolete and might be removed from a future version of simuPOP. Use type instead.\n')
                opt[key] = kwargs[key]
            if key == 'useDefault' and 'DBG_COMPATIBILITY' in simuOptions['Debug']:
                sys.stderr.write('WARNING: useDefault is obsolete and might be removed from a future version of simuPOP.\n')
        #
        if not opt.has_key('name') and name != '':
            opt['name'] = name
        if not opt.has_key('default'):
            opt['default'] = default
        self.options.append(opt)
        if opt.has_key('separator'):
            opt['gui_type'] = 'separator'
            return
        if 'name' not in opt.keys():
            raise ValueError('Item name cannot be ignored in an option specification dictionary')
        # allow alphabet, number and underscore (_).
        if not opt['name'].replace('_', '').isalnum() or not opt['name'][0].isalpha():
            raise ValueError('Invalid option name %s' % opt['name'])
        if kwargs.has_key('arg') and sum([x.has_key('arg') and x['arg'][0] == opt['arg'][0] for x in self.options]) > 1:
            raise ValueError("Duplicated short argument '%s'" % opt['arg'])
        if 'default' not in opt.keys() and 'separator' not in opt.keys():
            raise ValueError('A default value is not provided for option "%s"' % opt['name'])
        if opt.has_key('arg') and sum([x.has_key('arg') and x['arg'] == opt['arg'] for x in self.options]) > 1:
            raise ValueError("Duplicated short argument '%s'" % opt['arg'])
        if opt['name'] in reserved_options:
            raise ValueError("Option '--%s' is reserved. Please use another name." % opt['name'])
        if opt['name'] in methods:
            raise ValueError("Option '%s' conflicts with the '%s' member function of the Params class." % \
                (opt['name'], (opt['name'])))
        if opt['name'] in self.__dict__.keys():
            raise ValueError("Option '%s' conflicts with attribute '%s' of this Params object." % \
                (opt['name'], (opt['name'])))
        #
        # 
        # simuPOP 1.0.3 disallow invalid default value in parameter specification
        # dictionary. This seemed to be a logic change but it turned out that invalid
        # default value cannot always be avoided (e.g. a valid filename that cannot have
        # a valid default value). This version allows invalid default value again.
        #
        #if not _validate(opt['value'], opt, self.options):
        #    raise ValueError("Default value '%s' for option '%s' does not pass validation." % (str(opt['default']), opt['name']))
        opt['value'] = opt['default']
        opt['processed'] = False
        #
        if self.dict.has_key(opt['name']):
            raise ValueError('Option %s already exists.' % opt['name'])
        # process raw_type
        if opt.has_key('type'):
            if type(opt['type']) in (tuple, list) and \
                len(opt['type']) == 2 and opt['type'][0] == 'chooseOneOf':
                opt['gui_type'] = 'chooseOneOf'
                opt['chooseOneOf'] = opt['type'][1]
                if len(opt['chooseOneOf']) == 0:
                    raise ValueError('Empty list to choose from')
                if not opt.has_key('validator'):
                    opt['validator'] = valueOneOf(opt['chooseOneOf'])
                if not opt.has_key('allowedTypes'):
                    opt['allowedTypes'] = [type(opt['chooseOneOf'][0])]
            elif type(opt['type']) in (tuple, list) and \
                len(opt['type']) == 2 and opt['type'][0] == 'chooseFrom':
                opt['gui_type'] = 'chooseFrom'
                opt['chooseFrom'] = opt['type'][1]
                if len(opt['chooseFrom']) == 0:
                    raise ValueError('Empty list to choose from')
                if not opt.has_key('validator'):
                    opt['validator'] = valueOneOf(opt['chooseFrom'])
                if not opt.has_key('allowedTypes'):
                    opt['allowedTypes'] = [type(opt['chooseFrom'][0])]
            elif opt['type'] == 'integer':
                opt['gui_type'] = 'others'
                opt['allowedTypes'] = [int, long]
                if not opt.has_key('validator'):
                    opt['validator'] = valueIsInteger()
            elif opt['type'] == 'integers':
                opt['gui_type'] = 'others'
                opt['allowedTypes'] = [list, tuple]
                if not opt.has_key('validator'):
                    opt['validator'] = valueListOf(valueIsInteger())
            elif opt['type'] == 'number':
                opt['gui_type'] = 'others'
                opt['allowedTypes'] = [int, long, float]
                if not opt.has_key('validator'):
                    opt['validator'] = valueIsNum()
            elif opt['type'] == 'numbers':
                opt['gui_type'] = 'others'
                opt['allowedTypes'] = [list, tuple]
                if not opt.has_key('validator'):
                    opt['validator'] = valueListOf(valueIsNum())
            elif opt['type'] == 'string':
                opt['gui_type'] = 'others'
                opt['allowedTypes'] = [str]
            elif opt['type'] == 'strings':
                opt['gui_type'] = 'others'
                opt['allowedTypes'] = [list, tuple]
                if not opt.has_key('validator'):
                    opt['validator'] = valueListOf(str)
            elif opt['type'] == 'filename':
                opt['gui_type'] = 'browseFile'
                opt['allowedTypes'] = [str]
                #if not opt.has_key('validator'):
                #    opt['validator'] = valueValidFile()
            elif opt['type'] == 'dirname':
                opt['gui_type'] = 'browseDir'
                opt['allowedTypes'] = [str]
                #if not opt.has_key('validator'):
                #    opt['validator'] = valueValidDir()
            else:
                if isinstance(opt['type'], type(int)):
                    opt['allowedTypes'] = [opt['type']]
                elif type(opt['type']) in [list, tuple]:
                    # must be a list.
                    for v in opt['type']:
                        if not isinstance(v, type(int)):
                            raise ValueError('Key type is not one of the specified types, or a list of types.')
                    opt['allowedTypes'] = opt['type']            
                else:
                    raise ValueError('Key type is not one of the specified types, or a list of types.')
        # if raw_type is not specified.
        # determine type of option
        if opt.has_key('gui_type'):
            if opt['gui_type'] == 'boolean' and (opt.has_key('chooseOneOf') or opt.has_key('chooseFrom')):
	             raise ValueError('Directive chooseOneOf or chooseFrom can only be used for option that expects a value.');
        elif not (opt.has_key('label') or opt.has_key('separator')):
            opt['gui_type'] = 'hidden'
        elif opt.has_key('separator'):
            opt['gui_type'] = 'separator'
        elif opt.has_key('allowedTypes') and len(opt['allowedTypes']) == 1 and opt['allowedTypes'][0] == bool:
            opt['gui_type'] = 'boolean'
        elif opt.has_key('chooseFrom'):
            opt['gui_type'] = 'chooseFrom'
            if not opt.has_key('validator'):
                opt['validator'] = valueListOf(valueOneOf(opt['chooseFrom']))
        elif opt.has_key('chooseOneOf'):
            opt['gui_type'] = 'chooseOneOf'
            if not opt.has_key('validator'):
                opt['validator'] = valueOneOf(opt['chooseOneOf'])
        elif opt.has_key('validator') and callable(opt['validator']) and \
            opt['validator'].__doc__ == valueValidFile().__doc__:
            opt['gui_type'] = 'browseFile'
        elif opt.has_key('validator') and callable(opt['validator']) and \
            opt['validator'].__doc__ == valueValidDir().__doc__:
            opt['gui_type'] = 'browseDir'
        else:
            opt['gui_type'] = 'others'
        #
        if opt['gui_type'] == 'boolean' and opt['default'] is True and 'DBG_COMPATIBILITY' in simuOptions['Debug']:
            sys.stderr.write('WARNING: the default value for a boolean parameter should be False.\n')
        # is default value allowed?
        if not opt.has_key('allowedTypes'):
            if opt.has_key('chooseFrom'):
                opt['allowedTypes'] = [list, tuple]
            else:
                opt['allowedTypes'] = [type(opt['default'])]
        if opt.has_key('allowedTypes') and type(opt['default']) not in opt['allowedTypes']:
            if list in opt['allowedTypes']:
                opt['default'] = [opt['default']]
            elif tuple in opt['allowedTypes']:
                opt['default'] = (opt['default'],)
            else:
                raise ValueError('Default value "%s" is not of one of the allowed types.' % str(opt['default']))
        self.dict[opt['name']] = opt


    def saveConfig(self, file, params=[]):
        '''Write a configuration file to *file*. This file can be later read
        with command line option ``-c`` or ``--config``. All parameters with a
        ``label`` entry are saved unless a list of parameters are specified
        in *params*. In addition to parameter definitions, command lines
        options to specify the same set of parameters are saved to the
        configuration file.
        '''
        cfg = open(file,'w')
        try:
            # sys.argv[0] might not exist
            cfg.write("# configuration file for program " + sys.argv[0] + "\n")
        except:
            pass
        cfg.write("# Configuration saved at at " + time.asctime() + "\n")
        for opt in self.options:
            if opt.has_key('separator'):
                continue
            if len(params) > 0 and opt['name'] not in params:
                continue
            # no label, and is not specified in params
            if not opt.has_key('label') and len(params) == 0:
                continue
            cfg.write('\n')
            # write arg and long arg
            if opt.has_key('label'):
                cfg.write("# label:\t%s\n" % opt['label'])
            # write description
            if opt.has_key('description'):
                desc = opt['description'].splitlines()
                cfg.write("# description:\n")
                for d in desc:
                    cfg.write("#\t" + d.strip() + "\n")
            # write out option value, try to make it python readable
            cfg.write("%s = %s\n" % (opt['name'], _prettyString(opt['value'], quoted=True)))
        cfg.write("\n\n#The same options can be given by command line options (subject to minor changes)\n")
        cmd = "#    --gui=False "
        # shorter version
        scmd = "#    --gui=False "
        for opt in self.options:
            if opt.has_key('separator') or not opt.has_key('label'):
                continue
            defaultVal = opt.has_key('useDefault') and opt['useDefault'] and str(opt['value']) == str(opt['default'])
            if opt['gui_type'] != 'boolean':
                if ',' in str(opt['value']):    # has single quote
                    arg = " --" + opt['name'] + '=' + _prettyString(opt['value'], quoted=True)
                    cmd += arg
                    if not defaultVal:
                        scmd += arg
                else:
                    arg = " --" + opt['name'] + "=" + _prettyString(opt['value'], quoted=True)
                    cmd += arg
                    if not defaultVal:
                        scmd += arg
            elif opt['value']: # this option is True
                cmd += " --" + opt['name']
                if not defaultVal:
                    scmd += " --" + opt['name']
        cfg.write(' \\\n#    '.join(textwrap.wrap(cmd, break_long_words=False)) + "\n")
        # print out shorter version
        cfg.write("\n\n#Or a shorter version if default arguments are ignored\n")
        cfg.write(' \\\n#    '.join(textwrap.wrap(scmd, break_long_words=False)) + "\n")
        cfg.close()

    def loadConfig(self, file, params=[]):
        '''Load configuration from a file. If a list of parameters are
        specified in *params*, only these parameters will be processed.
        '''
        cfg = open(file)
        for line in cfg.readlines():
            if line.strip() == '' or line.startswith('#'):
                continue
            for opt in self.options:
                if opt.has_key('separator'):
                    continue
                if len(params) > 0 and opt['name'] not in params:
                    continue
                name = opt['name']
                scan = re.compile(name + r'\s*=\s*(.*)')
                if scan.match(line):
                    value = scan.match(line).groups()[0]
                    opt['value'] = _getParamValue(opt, value.strip('''"'\n'''), self.options)
                    opt['processed'] = True
        cfg.close()

    def processArgs(self, args=None, params=[]):
        '''try to get parameters from a list of arguments *args* (default to
        ``sys.argv``). If ``-h`` or ``--help`` is in *args*, this function
        prints out a usage message and returns ``False``. If a list of
        parameters are specified in *params*, only these parameters will be
        processed.
        '''
        if args is None:
            cmdArgs = sys.argv
        else:
            cmdArgs = args
        if '-h' in cmdArgs or '--help' in cmdArgs:
            print(self.usage())
            return False
        for opt in self.options:
            if opt['gui_type'] == 'separator':
                continue
            if len(params) > 0 and opt['name'] not in params:
                continue
            if opt['gui_type'] == 'boolean': # do not expect an argument, simple
                # this part is complex in order to be backward compatible.
                # In the new version, only --name is allowed.
                value = None
                indexes = []
                for idx,arg in enumerate(cmdArgs):
                    if arg == '--' + opt['name']:
                        value = True
                        indexes.append(idx)
                        if idx < len(cmdArgs) - 1 and cmdArgs[idx+1] in ['True', 'true', '1']:
                            indexes.append(idx+1)
                        elif idx < len(cmdArgs) - 1 and cmdArgs[idx+1] in ['False', 'false', '0']:
                            value = False
                            indexes.append(idx+1)
                    elif arg in ['--%s=%s' % (opt['name'], x) for x in  ['True', 'true', '1']]:
                        value = True
                        indexes.append(idx)
                    elif arg in ['--%s=%s' % (opt['name'], x) for x in  ['False', 'false', '0']]:
                        value = False
                        indexes.append(idx)
                    elif opt.has_key('arg') and arg == '-' + opt['arg']:
                        value = True
                        indexes.append(idx)
                        if idx < len(cmdArgs) - 1 and cmdArgs[idx+1] in ['True', 'true', '1']:
                            indexes.append(idx+1)
                        elif idx < len(cmdArgs) - 1 and cmdArgs[idx+1] in ['False', 'false', '0']:
                            value = False
                            indexes.append(idx+1)
                    elif opt.has_key('arg') and arg in ['-%s=%s' % (opt['arg'], x) for x in  ['True', 'true', '1']]:
                        value = True
                        indexes.append(idx)
                    elif opt.has_key('arg') and arg in ['-%s=%s' % (opt['arg'], x) for x in  ['False', 'false', '0']]:
                        value = False
                        indexes.append(idx)
                if value is not None:
                    for idx in indexes:
                        if idx in self.processedArgs:
                            raise ValueError("Parameter " + cmdArgs[idx] + " has been processed before.")
                        self.processedArgs.append(idx)
                    opt['value'] = value
                    opt['processed'] = True
                continue
            # this is a more complicated case
            name = opt['name']
            hasArg = [x.startswith('--%s=' % name) or x == '--' + name for x in cmdArgs]
            if True in hasArg:
                idx = hasArg.index(True)
                # case 1: --arg something
                if cmdArgs[idx] == '--' + name:
                    if idx in self.processedArgs or idx+1 in self.processedArgs:
                        raise ValueError("Parameter " + cmdArgs[idx] + " has been processed before.")
                    try:
                        val = _getParamValue(opt, cmdArgs[idx+1], self.options)
                        self.processedArgs.extend([idx, idx+1])
                        opt['value'] = val
                        opt['processed'] = True
                    except Exception as e:
                        print("ERROR: Failed to assign parameter %s with value '%s'" % (opt['name'],
                            cmdArgs[idx+1]))
                        print(e)
                        continue
                # case 2 --arg=something
                else:
                    try:
                        val = _getParamValue(opt, cmdArgs[idx][(len(name)+3):], self.options)
                        self.processedArgs.append(idx)
                        opt['value'] = val
                        opt['processed'] = True
                    except Exception as e:
                        print("ERROR: Failed to assign parameter %s with value '%s'" % (opt['name'],
                            cmdArgs[idx][(len(name)+3):]))
                        print(e)
                        continue
            if not opt.has_key('arg') or len(opt['arg']) == 0:
                continue
            hasArg = [x[:2] == '-' + opt['arg'][0] for x in cmdArgs]
            if True in hasArg:
                idx = hasArg.index(True)
                # has something like -a
                # case 1: -a file
                if cmdArgs[idx] == '-' + opt['arg'][0]:
                    if idx in self.processedArgs or idx+1 in self.processedArgs:
                        raise ValueError("Parameter " + cmdArgs[idx] + " has been processed before.")
                    try:
                        val = _getParamValue(opt, cmdArgs[idx+1], self.options)
                        self.processedArgs.extend([idx, idx+1])
                        opt['value'] = val
                        opt['processed'] = True
                    except Exception as e:
                        print("ERROR: Failed to assign parameter %s with value '%s'" % (opt['name'],
                            cmdArgs[idx+1]))
                        print(e)
                        continue
                # case 2: -aopt or -a=opt
                else:
                    if idx in self.processedArgs:
                        raise ValueError("Parameter " + cmdArgs[idx] + " has been processed before.")
                    try:
                        arg = cmdArgs[idx]
                        if len(arg) > 3 and arg[2] == '=':
                            val = _getParamValue(opt, cmdArgs[idx][3:], self.options)
                        else:
                            val = _getParamValue(opt, cmdArgs[idx][2:], self.options)
                        self.processedArgs.append(idx)
                        opt['value'] = val
                        opt['processed'] = True
                    except Exception as e:
                        print("ERROR: Failed to assign parameter %s with value '%s'" % (opt['name'],
                            cmdArgs[idx]))
                        print(e)
                        continue
        return True

    def usage(self, usage='usage: %prog [-opt [arg] | --opt [=arg]] ...'):
        '''Reutn the usage message from the option description list.
        ``'%prog'`` in parameter *usage* will be replaced by
        ``os.path.basename(sys.argv[0])``.
        '''
        return _usage(self.options, self.doc, self.details)

    def termGetParam(self, params=[]):
        '''Get parameters from interactive user input. By default, all
        parameters are processed unless one of the following conditions is
        met:

        1. Parameter without a label
        2. Parameter with ``useDefault`` set to ``True`` (deprecated)
        3. Parameter that have been determined from command line options or a
           configuration file
        4. Parameter that have been determined by a previous call of this
           function.

        If a list of parameters are given in *params*, these parameters are
        processed regardless the mentioned conditions.
        '''
        #
        for opt in self.options:
            if opt.has_key('separator'):
                continue
            if len(params) == 0 and (
                opt['processed'] or \
                (not opt.has_key('label')) or \
                ('DBG_BATCHTESTING' in simuOptions['Debug']) or \
                (opt.has_key('useDefault') and opt['useDefault'])):
                continue
            if len(params) > 0 and opt['name'] not in params:
                continue
            # prompt
            if opt.has_key('label'):
                prompt = '%s (%s): ' % (opt['label'], str(opt['value']))
            else:
                prompt = '%s (%s): ' % (opt['name'], str(opt['value']))
            while True:
                value = raw_input('\n' + prompt)
                if value == '':
                    # use existing value...
                    opt['processed'] = True
                    break
                try:
                    opt['value'] = _getParamValue(opt, value, self.options)
                    opt['processed'] = True
                    break
                except:
                    print("Invalid input.\n")
            if value is None and not (opt.has_key('allowedTypes') and type(None) in opt['allowedTypes']):
                # should have a valid value now.
                return False
        # successfully handled all parameters
        return True

    def guiGetParam(self, nCol = None, gui=None):
        '''Get parameter from a ``Tkinter`` or ``wxPython`` dialog. The
        parameter will try to arrange parameters optimaly but you can also
        set the number of columns using parameter *nCol*. If both GUI toolkits
        are available, ``wxPython`` will be used unless *gui* is set to
        ``Tkinter``. If none of the toolkits are available, this function will
        raise an ``ImportError``.

        If ``Params.valueValidFile`` or ``Params.valueValidDir`` is used
        to validate a parameter, double click the text input box of this
        parameter will open a file or directory browse dialog.
        '''
        title = os.path.split(sys.argv[0])[-1]
        if gui != 'Tkinter':
            try:
                return _wxParamDialog(self.options, title, self.doc, self.details, nCol).getParam()
            except ImportError:
                # continue to try tk
                pass
        # try tkinter
        return _tkParamDialog(self.options, title, self.doc, self.details, nCol).getParam()

    # get parameter
    def getParam(self, gui=None, nCol=None, configFile=None, args=None,
        checkArgs=True):
        '''Get parameters from commandline option, configuration file, a
        parameter input dialog and from interactive user input.

        gui
            Whether or not use a dialog and which graphical toolkit to use.
            Global gui setting is used by default but you can also set this
            parameter to ``True``, ``False``, ``Tkinter`` or ``wxPython`` to
            override the global setting.

        nCol
            Number of columns in the parameter input dialog. This is usual
            determine automatically depending on the number of options.

        configFile
            Configuration file from which to load values of parameters. If
            unspecified, it will be determined from command line option
            ``--config``.

        args
            Command line arguments are obtained from ``sys.argv`` unless a
            list of options are provided in this argument.

        checkArgs
            This function by default checks if all commandline arguments have
            been processed, you can set ``chekArgs`` to ``False`` if some of
            the arguments are intended to be processed separately.

        '''
        if args is None:
            cmdArgs = sys.argv
        else:
            cmdArgs = args
        #
        if gui is None:
            # command line option --gui should have been processed...
            self.gui = simuOptions['GUI']
        else:
            self.gui = gui
        #
        # Start processing
        #
        self.processedArgs = []
        # first assign values from non-GUI sources
        if configFile is not None:
            self.loadConfig(configFile)
        elif '--config' in cmdArgs:
            idx = cmdArgs.index('--config')
            if idx == len(cmdArgs) - 1:
                raise exception.ValueError('Expect a filename after --config')
            self.loadConfig(cmdArgs[idx+1])
        elif True in [x.startswith('--config=') for x in cmdArgs]:
            idx = [x.startswith('--config=') for x in cmdArgs].index(True)
            file = cmdArgs[idx][len('--config='):]
            for quote in ['"', "'"]:
                if file.startswith(quote) and file.endswith(quote):
                    file = file[1:-1]
            self.loadConfig(file)
        #
        if not self.processArgs(cmdArgs):
            # encounter -h or --help
            return False
        #
        if checkArgs:
            # look if any argument was not processed
            for i in range(1, len(cmdArgs)):
                if i in self.processedArgs:
                    continue
                elif cmdArgs[i] in ['-h', '--help', '--optimized', '--config', '--gui']:
                    continue
                elif cmdArgs[i].startswith('--config=') or cmdArgs[i].startswith('--gui='):
                    continue
                elif i > 0 and cmdArgs[i-1] in ['--config', '--gui']:
                    continue
                raise ValueError('Command line argument %s is not process.' % cmdArgs[i] +
                    'You may have misspelled the argument name or passed it an invalid value.')
        #
        if self.gui in [False, 'interactive']:
            return self.termGetParam()
        elif self.gui == 'batch':
            # valid values because some default values can be false
            for opt in self.options:
                if opt.has_key('separator'):
                    continue
                # default value might not have the required type (e.g. 9 for [9])
                opt['value'] = _getParamValue(opt, opt['value'], self.options)
                if opt.has_key('allowedTypes') and type(opt['value']) not in opt['allowedTypes']:
                    raise ValueError("Value '%s' is not of allowed type for parameter '%s'." % \
                        (str(opt['value']), opt['name']))
                if not _validate(opt['value'], opt, self.options):
                    raise ValueError("Value '%s' is not allowed for parameter '%s'." % \
                        (str(opt['value']), opt['name']))
            return True
        # GUI
        try:
            return self.guiGetParam(nCol, gui=self.gui)
        except ImportError:
            return self.termGetParam()
        return False

    def asDict(self):
        '''Return parameters as a dictionary.'''
        res = {}
        for opt in self.options:
            if opt.has_key('separator'):
                continue
            name = opt['name']
            res[name] = opt['value']
        return res

    def asList(self):
        '''Return parameters as a list.'''
        res = []
        for opt in self.options:
            if opt.has_key('separator'):
                continue
            res.append(opt['value'])
        return res
