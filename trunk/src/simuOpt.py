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


class ``simuOpt.simuOpt`` provides a powerful way to handle commandline
arguments. Briefly speaking, a ``simuOpt`` object can be created from a list
of parameter specification dictionaries. The parameters are then become
attributes of this object. A number of functions are provided to determine
values of these parameters using commandline arguments, a configuration
file, or a parameter input dialog (using ``Tkinter`` or ``wxPython``).
Values of these parameters can be accessed as attributes, or extracted
as a list or a dictionary. Note that the ``simuOpt.getParam`` function
automatically handles the following commandline arguments.

``-h`` or ``--help``
    Print usage message.

``--config=configFile``
    Read parameters from a configuration file *configFile*.

'''

# First try to get environmental variable

import os, sys, exceptions, types, re, time, textwrap


#
# define some validataion functions
#
def valueNot(t):
    '''Return a function that returns true if passed option does not passes
    validator t'''
    def func(val):
        if type(t) == types.FunctionType:
            return not t(val)
        else:
            raise exceptions.ValueError("We expect a function valueXXX")
    return func


def valueOr(t1, t2):
    '''Return a function that returns true if passed option passes validator
    t1 or t2'''
    def func(val):
        if type(t1) == types.FunctionType and type(t2) == types.FunctionType:
            return t1(val) or t2(val)
        else:
            raise exceptions.ValueError("We expect a function valueXXX")
    return func


def valueAnd(t1, t2):
    '''Return a function that returns true if passed option passes validator
    t1 and t2'''
    def func(val):
        if type(t1) == types.FunctionType and type(t2) == types.FunctionType:
            return t1(val) and t2(val)
        else:
            raise exceptions.ValueError("We expect a function valueXXX")
    return func


def valueOneOf(t):
    '''Return a function that returns true if passed option is one of the values
    list in t'''
    if not type(t) in [types.ListType, types.TupleType]:
        raise exceptions.ValueError('argument of valueOneOf should be a list')
    def func(val):
        yes = False
        for item in t:
            if item == val:    # equal value
                return True
            if type(item) == types.FunctionType: # a test function
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


def valueIsNum():
    'Return a function that returns true if passed option is a number (int, long or float)'
    def func(val):
        return type(val) in [types.IntType, types.LongType, types.FloatType]
    return func


def valueIsList():
    'Return a function that returns true if passed option is a list (or tuple)'
    def func(val):
        return type(val) in [types.ListType, types.TupleType]
    return func


def valueValidDir():
    '''Return a function that returns true if passed option val if a valid
    directory'''
    def func(val):
        return os.path.isdir(val)
    return func


def valueValidFile():
    '''Return a function that returns true if passed option val if a valid
    file'''
    def func(val):
        return os.path.isfile(val)
    return func


def valueListOf(t):
    '''Return a function that returns true if passed option val is a list of
    type t. If t is a function (validator), check if all v in val pass t(v)
    '''
    def func(val):
        if not type(val) in [types.ListType, types.TupleType]:
            return False
        if type(t) in [types.ListType, types.TupleType]:
            for i in val:
                if not type(i) in t:
                    return False
        elif type(t) == types.FunctionType:
            for i in val:
                if not t(i):
                    return False
        else:
            for i in val:
                if type(i) != t:
                    return False
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
    if type(value) in [types.ListType, types.TupleType] and len(value)>1:
        txt = '[' + ', '.join([_prettyString(x, True, False) for x in value]) + ']'
        if outer:
            return quote(txt)
        else:
            return txt
    elif type(value) == types.StringType:
        return quote(value)
    elif outer:
        return quote(str(value))
    else:
        return str(value)


def _prettyDesc(text, indent=''):
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
            txt.extend(textwrap.wrap(blk[1:], width=80,
                initial_indent=indent + blk_indent, # keep spaces after | character
                subsequent_indent=indent))
        else:
            txt.extend(textwrap.wrap(blk, width=80, initial_indent=indent,
                subsequent_indent=indent))
    return '\n'.join(txt)


def _usage(options, msg='', usage='usage: %prog [-opt [arg] | --opt [=arg]] ...'):
    'Return a usage message.'
    if msg != '':
        message = msg + '\n'
    else:
        message = ''
    message += '''%s

options:
  -h, --help            show this help message and exit
  --config ARG          load parameters from ARG
  --optimized           run the script using an optimized simuPOP module
  --gui ARG             which graphical toolkit to use
''' % usage.replace('%prog', os.path.basename(sys.argv[0]))
    for opt in options:
        if opt.has_key('separator'):
            continue
        name = ''
        if opt.has_key('arg'):
            if opt['arg'].endswith(':'):
                name += '-%s ARG, ' % opt['arg'].rstrip(':')
            else:
                name += '-%s, ' % opt['arg']
        #
        if opt['longarg'].endswith('='):
            name += '--%s ARG' % opt['longarg'].rstrip('=')
        else:
            name += '--%s' % opt['longarg']
        #
        if opt.has_key('label'):
            label = opt['label'] + ' '
        else:
            label = ''
        if len(name) >= 22:
            name += '\n' + ' '*24
        else:
            name = '%-21s ' % name
        if opt['default'] == '':
            defaultVal = "''"
        elif opt['default'] is None:
            defaultVal = 'None'
        else:
            defaultVal = _prettyString(opt['default'])
        message += '  %s%s[default: %s ]\n' % (name, label, defaultVal)
        if opt.has_key('description'):
            message += _prettyDesc(opt['description'], indent=' '*24) + '\n'
    return message

def _getParamValue(p, val):
    ''' try to get a value from value, raise exception if error happens. '''
    if p.has_key('separator'):
        raise exceptions.ValueError('Cannot get a value for separator')
    # if we are giving a unicode string, convert!
    if type(val) == types.UnicodeType:
        val = str(val)
    # remove quotes from string?
    if p.has_key('allowedTypes') and type('') in p['allowedTypes']:
        for quote in ['"', "'", '"""', "'''"]:
            if val.startswith(quote) and val.endswith(quote):
                val = val[len(quote):-len(quote)]
                break
    if (not p.has_key('allowedTypes')) or type(val) in p['allowedTypes']:
        if p.has_key('validate') and not p['validate'](val):
                raise exceptions.ValueError("Value '%s' is not allowed for parameter %s" % \
                    (str(val), p['longarg'].rstrip('=')))
        return val
    # handle another 'auto-boolean' case
    elif not (p['longarg'].endswith('=')):
        if val in ['1', 'true', 'True']:
            return True
        elif val in ['0', 'false', 'False']:
            return False
        else:
            raise exceptions.ValueError('Expect 0/1, true/false for boolean values for parameter %s ' % p['longarg'].rstrip('='))
    # other wise, need conversion
    if type(val) in [types.StringType, types.UnicodeType]:
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
        if p.has_key('validate') and not p['validate'](val):
                raise exceptions.ValueError("Value "+str(val)+' does not pass validation')
        return val
    elif types.ListType in p['allowedTypes'] or types.TupleType in p['allowedTypes']:
        if p.has_key('validate') and not p['validate']([val]):
                raise exceptions.ValueError("Value "+str([val])+' does not pass validation')
        return [val]
    elif type(val) == type(True) and types.IntType in p['allowedTypes']: # compatibility problem
        return val
    elif type[val] == types.UnicodeType and types.StringType in p['allowedTypes']:
        return str(val)
    else:
        raise exceptions.ValueError('Type of input parameter ' + str(val) + " is incorrect. (param " \
            + p.setdefault('longarg','none') +")")


class _paramDialog:
    def __init__(self, options, title = '', description='', details='', nCol=1):
        if len(options) == 0:
            raise exceptions.ValueError("Empty field names...")    # some behaviors
        #
        # now, initialize variables
        self.options = options
        self.title = title
        self.description = description
        self.details = details
        if nCol is None:
            self.nCol = len(self.options)/20 + 1
        else:
            self.nCol = nCol


    def getNumOfRows(self, multiLineChooseOneOf=False):
        '''Count the number of rows that is needed for all parameters'''
        row = 0
        for opt in self.options:
            if opt.has_key('label') or opt.has_key('separator'):
                row += 1
            if opt.has_key('chooseFrom'):
                row += len(opt['chooseFrom']) - 1
            elif opt.has_key('chooseOneOf') and multiLineChooseOneOf:
                row += len(opt['chooseOneOf']) - 1
        if row / self.nCol * self.nCol == row:
            row /= self.nCol
        else:
            row = row/self.nCol + 1
        # it is possible but the row'th row sits between chooseOneOf or chooseFrom ...
        r = 0
        for opt in self.options:
            if opt.has_key('label') or opt.has_key('separator'):
                r += 1
            if opt.has_key('chooseFrom'):
                r += len(opt['chooseFrom']) - 1
            elif opt.has_key('chooseOneOf') and multiLineChooseOneOf:
                r += len(opt['chooseOneOf']) - 1
            if r >= row:
                row = r
                # starts new
                r = 0
        return row

    def createDialog(self):
        raise exceptions.SystemError('Please define createDialog')

    def runDialog(self):
        raise exceptions.SystemError('Please define runDialog')

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
        import Tkinter as tk
        import tkFont
        globals()['tk'] = tk
        globals()['tkFont'] = tkFont

    def denyWindowManagerClose(self):
        '''Don't allow WindowManager close'''
        x = tk.Tk()
        x.withdraw()
        x.bell()
        x.destroy()

    def onHelpOK(self, event):
        self.helpDlg.destroy()

    def createHelpDialog(self):
        self.helpDlg = tk.Toplevel(self.app)
        self.helpDlg.title('Help for ' + self.title)
        #
        msg = tk.Text(self.helpDlg, wrap=tk.WORD)
        msg.insert(tk.END, _usage(self.options, self.details))
        msg.grid(column=0, row=0, pady=10, padx=10,
            sticky = tk.E + tk.W + tk.N + tk.S)
        # scrollbar
        sb = tk.Scrollbar(self.helpDlg)
        sb.config(command=msg.yview)
        sb.grid(column=1, row=0, sticky=tk.N + tk.S,
            padx=0, pady=10)
        msg["yscrollcommand"] = sb.set
        self.helpDlg.columnconfigure(0, weight=1)
        self.helpDlg.rowconfigure(0, weight=1)
        self.helpDlg.rowconfigure(1, weight=0)
        okButton = tk.Button(self.helpDlg, takefocus=1, text="OK")
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
                continue
            try:
                # get text from different type of entries
                if self.entryWidgets[g].winfo_class() == "Entry":    # an entry box?
                    val = _getParamValue(self.options[g], self.entryWidgets[g].get())
                elif self.entryWidgets[g].winfo_class() == "Listbox":    # a listbox
                    sel = self.entryWidgets[g].curselection()
                    if len(sel) == 1:
                        items = self.entryWidgets[g].get(sel)
                    else:
                        items = []
                        for s in sel:
                            items.append(self.entryWidgets[g].get(s))
                    val = _getParamValue(self.options[g], items)
                elif self.entryWidgets[g].winfo_class() == "Checkbutton":    # a checkbutton (true or false)
                    # gets 0/1 for false/true
                    var = self.options[g]['value'].get()
                    val = _getParamValue(self.options[g], var == 1)
            except Exception,e:
                for lab in self.labelWidgets:
                    if lab is not None:
                        lab.configure(fg='black')
                # set this one to red
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
        self.app = tk.Tk()
        self.app.protocol('WM_DELETE_WINDOW', self.denyWindowManagerClose)
        self.app.title(self.title)
        #
        # the main window
        self.entryWidgets = [None]*len(self.options)
        self.labelWidgets = [None]*len(self.options)
        # all use grid management
        # top message
        topMsg = tk.Label(self.app, text=self.description)
        topMsg.grid(row=0, column=0, columnspan = 2 * self.nCol, sticky=tk.E + tk.W,
            padx=10, pady=5)
        # find out number of items etc
        numRows = self.getNumOfRows(True)
        rowIndex = 0
        # all entries
        for g,opt in enumerate(self.options):
            if not (opt.has_key('label') or opt.has_key('separator')):
                continue
            colIndex = rowIndex / numRows
            # use different entry method for different types
            if opt.has_key('separator'):
                self.labelWidgets[g] = tk.Label(self.app, text=opt['separator'],
                    font=tkFont.Font(size=12, weight='bold'))
                self.labelWidgets[g].grid(column=colIndex*2, row= rowIndex % numRows + 1,
                    ipadx=0, padx=0, sticky=tk.E)
                emptyInput = tk.Label(self.app, text='')
                emptyInput.grid(column = colIndex*2 + 1, row= rowIndex % numRows + 1,
                    ipadx=0, padx=0)
                self.entryWidgets[g] = None
                rowIndex += 1
                continue
            value = self.options[g]['value']
            if value is None:
                value = self.options[g]['default']
            if opt.has_key('chooseOneOf'):    # single choice
                height = len(opt['chooseOneOf'])
                self.labelWidgets[g] = tk.Label(self.app, text=opt['label'])
                self.labelWidgets[g].grid(column=colIndex*2, row=rowIndex%numRows+1,
                    padx=5, rowspan = height, sticky=tk.E)
                self.entryWidgets[g] = tk.Listbox(self.app, selectmode=tk.SINGLE,
                    exportselection=0, height = height)
                self.entryWidgets[g].grid(column=colIndex*2+1, row=rowIndex%numRows+1,
                    padx=5, rowspan = height)
                rowIndex += height
                for entry in opt['chooseOneOf']:
                    self.entryWidgets[g].insert(tk.END, str(entry))
                if value is not None:
                    self.entryWidgets[g].select_set(opt['chooseOneOf'].index(value))
            elif opt.has_key('chooseFrom'):    # multiple choice
                height = len(opt['chooseFrom'])
                self.labelWidgets[g] = tk.Label(self.app, text=opt['label'])
                self.labelWidgets[g].grid(column=colIndex*2, row=rowIndex%numRows+1,
                    padx=5, rowspan = height, sticky=tk.E)
                self.entryWidgets[g] = tk.Listbox(self.app, selectmode=tk.EXTENDED,
                    exportselection=0, height = height)
                self.entryWidgets[g].grid(column=colIndex*2+1, row=rowIndex%numRows+1,
                    padx=5, rowspan = height)
                rowIndex += height
                for entry in opt['chooseFrom']:
                    self.entryWidgets[g].insert(tk.END, str(entry))
                if value is not None:
                    if type(value) in [types.TupleType, types.ListType]:
                        for val in value:
                            self.entryWidgets[g].select_set( opt['chooseFrom'].index(val))
                    else:
                        self.entryWidgets[g].select_set( opt['chooseFrom'].index( value ))
            elif (opt.has_key('arg') and opt['arg'][-1] != ':') or \
                 (opt.has_key('longarg') and opt['longarg'][-1] != '='):  # true or false
                self.labelWidgets[g] = tk.Label(self.app, text=opt['label'])
                self.labelWidgets[g].grid(column=colIndex*2, row=rowIndex%numRows+1, padx=10,
                    rowspan = 1, sticky=tk.E)
                # replace value by a tk IntVar() because tk.Checkbutton has to store
                # its value in such a variable. value.get() will be used to return the
                # state of this Checkbutton.
                # c.f. http://infohost.nmt.edu/tcc/help/pubs/tkinter/control-variables.html
                iv = tk.IntVar()
                iv.set(self.options[g]['value'] == True) # value can be None, True or False
                self.options[g]['value'] = iv
                self.entryWidgets[g] = tk.Checkbutton(self.app, height=1,
                    text = "Yes / No", variable=self.options[g]['value'])
                self.entryWidgets[g].grid(column=colIndex*2+1, row=rowIndex%numRows+1, padx=5,
                    sticky=tk.W)
                rowIndex += 1
                self.entryWidgets[g].deselect()
            else:
                self.labelWidgets[g] = tk.Label(self.app, text=opt['label'])
                self.labelWidgets[g].grid(column=colIndex*2, row=rowIndex%numRows+1,
                    padx=5, sticky=tk.E)
                self.entryWidgets[g] = tk.Entry(self.app)
                self.entryWidgets[g].grid(column=colIndex*2+1, row=rowIndex%numRows+1,
                    padx=5, ipadx=0)
                rowIndex += 1
                 # put default value into the entryWidget
                if value is not None:
                    self.entryWidgets[g].insert(0, _prettyString(value))
            self.entryWidgets[g].bind("<Return>", self.onOK)
            self.entryWidgets[g].bind("<Escape>", self.onCancel)
        # help button
        helpButton = tk.Button(self.app, takefocus=1, text="Help")
        helpButton.bind("<Return>", self.onHelp)
        helpButton.bind("<Button-1>", self.onHelp)
        helpButton.grid(column=0, columnspan=self.nCol, row = numRows+1, pady=20)
        # ok button
        okButton = tk.Button(self.app, takefocus=1, text="Run!")
        okButton.bind("<Return>", self.onOK)
        okButton.bind("<Button-1>", self.onOK)
        okButton.grid( column=self.nCol, columnspan=self.nCol, row = numRows+1, pady=20)
        # cancel button
        cancelButton = tk.Button(self.app, takefocus=1, text="Cancel")
        cancelButton.bind("<Return>", self.onCancel)
        cancelButton.bind("<Button-1>", self.onCancel)
        cancelButton.grid( column=0, columnspan=2*self.nCol, row = numRows+1, pady=20)
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
        globals()['wx'] = wx

    def onHelp(self, event):
        # open another dialog
        helpDlg = wx.Dialog(parent=self.dlg, id=-1, title='Help for ' + self.title)
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(wx.TextCtrl(parent=helpDlg, id=-1, size=[600,400],
            style=wx.TE_MULTILINE | wx.TE_READONLY,
            value=_usage(self.options, self.details)), 0, wx.ALL, 20)
        self.addButton(wx.ID_OK, "OK", lambda event:helpDlg.EndModal(wx.ID_OK), helpDlg, box)
        helpDlg.SetSizerAndFit(box)
        helpDlg.Layout()
        helpDlg.ShowModal()
        helpDlg.Destroy()

    def onOK(self, event):
        ''' get result and convert values '''
        for g in range(len(self.entryWidgets)):
            if self.entryWidgets[g] == None:
                continue
            try:
                # get text from different type of entries
                try:    # an entry box or check box
                    val = _getParamValue(self.options[g], self.entryWidgets[g].GetValue())
                except:
                    try:    # a list box?
                        val = _getParamValue(self.options[g],
                            self.options[g]['chooseOneOf'][int(self.entryWidgets[g].GetSelection())])
                    except: # a checklist box?
                        items = []
                        for s in range(len(self.options[g]['chooseFrom'])):
                            if self.entryWidgets[g].IsChecked(s):
                                items.append(self.options[g]['chooseFrom'][s])
                        val = _getParamValue(self.options[g], items)
            except exceptions.Exception, e:
                # incorrect value
                # set to red
                # clear other red colors
                for lab in self.labelWidgets:
                    if lab is not None:
                        lab.SetForegroundColour('black')
                # set this one to red
                self.labelWidgets[g].SetForegroundColour('red')
                self.entryWidgets[g].SetFocus()
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
        box.Add(button, 0, wx.ALIGN_CENTER)
        self.dlg.Bind(wx.EVT_BUTTON, func, button)

    def createDialog(self):
        self.app = wx.PySimpleApp(0)
        self.dlg = wx.Dialog(parent=None, id=-1, title = self.title)
        self.entryWidgets = [None]*len(self.options)
        self.labelWidgets = [None]*len(self.options)
        #
        # the main window
        box = wx.BoxSizer(wx.VERTICAL)
        # do not use a very long description please
        topLabel = wx.StaticText(parent=self.dlg, id=-1, label='\n' + self.description)
        box.Add(topLabel, 0, wx.EXPAND | wx.LEFT | wx.ALIGN_CENTER, 50)
        topLabel.SetFont(wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL))
        # add a box for all ...
        paraBox = wx.FlexGridSizer(cols = self.nCol)
        for col in range(self.nCol):
            paraBox.AddGrowableCol(col)
        #
        # add several FlexGridSizer
        gridBox = []
        for col in range(self.nCol):
            gridBox.append(wx.FlexGridSizer(cols=2, vgap=2, hgap=5))
            gridBox[-1].AddGrowableCol(0)
            gridBox[-1].AddGrowableCol(1)
            paraBox.Add(gridBox[-1], 1, wx.EXPAND | wx.ALL, 10)
        box.Add(paraBox, 1, wx.EXPAND | wx.ALL, 5)
        # count numbers of valid parameters..
        # chooseFrom count as many
        numRows = self.getNumOfRows()
        rowIndex = 0
        # all entries
        for g,opt in enumerate(self.options):
            if not (opt.has_key('label') or opt.has_key('separator')):
                continue
            colIndex = rowIndex / numRows
            if opt.has_key('separator'):
                self.labelWidgets[g] = wx.StaticText(parent=self.dlg, id=-1, label=opt['separator'])
                self.labelWidgets[g].SetFont(wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
                gridBox[colIndex].Add(self.labelWidgets[g], 0, wx.ALIGN_LEFT )
                # no entry widget
                self.entryWidgets[g] = None
                # add an empty string to the right... it would be good to extend.
                gridBox[colIndex].Add(wx.StaticText(parent=self.dlg, id=-1, label=''), 1, wx.ALIGN_LEFT )
                rowIndex += 1
                continue
            value = self.options[g]['value']
            if value is None:
                value = self.options[g]['default']
            # label
            self.labelWidgets[g] = wx.StaticText(parent=self.dlg, id=-1, label=opt['label'])
            gridBox[colIndex].Add(self.labelWidgets[g], 0, wx.ALIGN_LEFT )
            # use different entry method for different types
            if opt.has_key('chooseOneOf'):    # single choice
                self.entryWidgets[g] = wx.Choice(parent=self.dlg, id=g, choices = opt['chooseOneOf'])
                if opt.has_key('description'):
                    self.entryWidgets[g].SetToolTipString(_prettyDesc(opt['description']))
                gridBox[colIndex].Add(self.entryWidgets[g], 1, wx.EXPAND )
                # if an value is given through command line argument or configuration file
                if value is not None:
                    try:
                        self.entryWidgets[g].SetSelection(opt['chooseOneOf'].index(value))
                    except:
                        raise ValueError('Value: %s is not one of %s.' % (str(value), str(opt['chooseOneOf'])))
                rowIndex += 1
            elif opt.has_key('chooseFrom'):    # multiple choice
                # the height is a little bit too much...
                self.entryWidgets[g] = wx.CheckListBox(parent=self.dlg, id=g,
                    choices = opt['chooseFrom'])
                if opt.has_key('description'):
                    self.entryWidgets[g].SetToolTipString(_prettyDesc(opt['description']))
                if value is not None:
                    if type(value) in [types.ListType, types.TupleType]:
                        for val in value:
                            self.entryWidgets[g].Check(opt['chooseFrom'].index(val))
                    else:
                        self.entryWidgets[g].Check(opt['chooseFrom'].index(value))
                gridBox[colIndex].Add(self.entryWidgets[g], 1, wx.EXPAND)
                rowIndex += len(opt['chooseFrom'])
            elif (opt.has_key('arg') and opt['arg'][-1] != ':') or \
                 (opt.has_key('longarg') and opt['longarg'][-1] != '='):  # true or false
                self.entryWidgets[g] = wx.CheckBox(parent=self.dlg, id=g, label = 'Yes / No')
                if opt.has_key('description'):
                    self.entryWidgets[g].SetToolTipString(_prettyDesc(opt['description']))
                if value is not None:
                    self.entryWidgets[g].SetValue(value)
                gridBox[colIndex].Add(self.entryWidgets[g], 1, wx.EXPAND)
                rowIndex += 1
            else: # an edit box
                # put default value into the entryWidget
                txt = ''
                if value is not None:
                    txt = _prettyString(value)
                self.entryWidgets[g] = wx.TextCtrl(parent=self.dlg, id=g, value=txt)
                if opt.has_key('description'):
                    self.entryWidgets[g].SetToolTipString(_prettyDesc(opt['description']))
                gridBox[colIndex].Add(self.entryWidgets[g], 1, wx.EXPAND )
                rowIndex += 1
        # help button
        buttonBox = wx.GridSizer(cols=3)
        self.addButton(wx.ID_HELP, 'Help', self.onHelp, self.dlg, buttonBox)
        self.addButton(wx.ID_CANCEL, 'Cancel', self.onCancel, self.dlg, buttonBox)
        self.addButton(wx.ID_OK, 'OK', self.onOK, self.dlg, buttonBox)
        box.Add(buttonBox, 0, wx.ALL | wx.EXPAND, 20)
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


class simuOpt:
    '''
    class simuOpt provides a uniform interface for simuPOP scripts to handle
    parameters. It allows users to get parameters from command line options,
    a configuration file, a parameter input dialog (*tkInter* or *wxPython*) or
    from interative input. This class provides parameter validation, conversion
    and and some utility functions to print, save and restore parameters.

    A simuOpt object accepts a parameter specification list that consists of
    dictionaries with pre-defined keys. Each item defines an option in terms of
    command line option, entry name in a configuration file, label in a
    parameter input dialog, acceptable types, validation rules and a default
    value. The following keys are currently supported:

    arg
        Short command line option name. For example ``'c'`` checks the presence
        of argument ``-c``. If a value is expected, a comma should be appened
        to the option name. For example, ``'p:'`` matches command line option
        ``-p=100`` or ``-p 100``. An options that does not expect a value is
        displayed in the parameter input dialog as an on/off switch.

    longarg
        Long command line option name.  For example ``'version'``  checks the
        presence of argument ``--version``. A equal character should be
        appended to the option name if a value is expected. For example,
        ``'mu='`` matches command line option ``--mu=0.001`` or ``--mu 0.001``.
        **This item defines the name of an option and cannot be ignored**.
        An options that does not expect a value is displayed in the parameter
        input dialog as an on/off switch.

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

    useDefault
        Use default value without asking, if the value can not be determined
        from GUI, command line option or config file. This is usually used for
        options that rarely need to be changed. Setting ``useDefault`` to such
        options simplifies user input.

    description
        A long description of this parameter. This description will be put into
        the usage information, and as parameter tooltip in the parameter input
        dialog. This string will be reformatted when it is written to a usage
        string (remove newlines and extra spaces and re-indent), with the
        exception that **lines with '|' as the first non-space/tab character
        will be outputed as is without the leading '|' symbol**.

    allowedTypes
        A list of acceptable types of this option. class ``simuOpt`` will try
        to convert user input to these types. For example, if ``allowedTypes``
        is ``types.ListType`` or  ``types.TupleType`` and the user's input is a
        scalar, the input will be converted to a list automatically. An option
        will not be accepted if such conversion fails.

    validate
        A function to validate the parameter. The function will be applied to
        user input. The option will not be accepted if this function returns
        ``False``. This module defines a large number of such validation
        functions but user defined functions are also acceptable.

    chooseOneOf
        If specified, a list of specified values will be displayed in the
        parameter input dialog and users are allowed to choose one of them.

    chooseFrom
        If specified, a list of specified values will be displayed in the
        parameter input dialog and users are allowed to choose one or more
        of them.

    separator
        This item specifies a separator (group header) in the parameter input
        dialog. All other fields are ignored.
   
    Not all keys need to be specified in each option description. Missing
    values are handled using some internal rules. For example, items without
    a ``label`` will not be displayed on the parameter dialog. This will
    effectively *hide* a parameter although users who know this parameter
    can set it using command line options.

    The ``simuOpt.simuOpt`` class defines a number of functions to collect,
    validate, and manipulate parameters using this parameter specification
    list.

    As a shortcut to create a simuOpt object with a number of attributes,
    a simuOpt object can be created with additional ``key=value`` pairs that
    could be assessed as attributes. This is used to create a ``simuOpt``
    object in which *parameters* are assigned directly.
    '''
    def __init__(self, options=[], doc='', details='', **kwargs):
        '''Create a ``simuOpt`` oject using a list of parameter specification
        dictionaries *options*. Additional *doc* and *details* can be
        specified which will be displayed as script summary (on the top of
        a parameter input dialog) and script introduction (the first part of
        a help message), respectively. Additional attributes could be assigned
        to a ``simuOpt`` object as keyword arguments. Note that it is customary
        to use module document (the first string object in a Python script) as
        *details*, using parameter ``details=__doc__``.
        '''
        #
        # validate
        if type(options) != type([]):
            raise exceptions.ValueError('An option specification list is expected')
        #
        self.options = []
        self.dict = {}
        for opt in options:
            if type(opt) != type({}):
                raise exceptions.ValueError('An option specification dictionary is expected')
            self.addOption(**opt)
        #
        if type(doc) != type(''):
            raise exceptions.ValueError('Parameter doc must be a string.')
        if type(details) != type(''):
            raise exceptions.ValueError('Parameter details must be a string.')
        self.doc = doc
        self.details = details
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
        raise exceptions.AttributeError('Can not locate attribute %s.' % name)

    def __setattr__(self, name, value):
        'Set the value of a parameter as an attribute.'
        if self.__dict__.has_key('dict') and self.dict.has_key(name):
            self.dict[name]['value'] = value
        else:
            self.__dict__[name] = value

    def addOption(self, pos=-1, **kwargs):
        '''
        Append an entry to the parameter specification list. Dictionary
        entries should be specified as keyword arguments such as
        ``longarg='option='``. More specifically, you can specify parameters
        ``arg``, ``longarg`` (required), ``label``, ``allowedTypes``,
        ``useDefault``, ``default`` (required), ``description``, ``validate``,
        ``chooseOneOf``, ``chooseFrom`` and ``separator``. This option will
        have a name specified by ``longarg`` (without optional trailing ``=``)
        and an initial default value specified by ``default``.

        An optional parameter *pos* can be given to specify an index before
        which this option will be inserted.
        '''
        allowed_keys = ['arg', 'longarg', 'label', 'allowedTypes',
            'useDefault', 'default', 'description',
            'validate', 'chooseOneOf', 'chooseFrom', 'separator']
        #
        methods = ['asDict', 'asList', 'getParam', 'loadConfig', 'saveConfig',
            'usage', 'processArgs', 'guiGetParam', 'termGetParam', 'addOption']
        #
        reserved_options = ['optimized', 'gui', 'config', 'help']
        #
        opt = {}
        for key in kwargs:
            if key in allowed_keys:
                opt[key] = kwargs[key]
            else:
                raise exceptions.ValueError('Invalid option specification key %s' % key)
        #
        if pos >= 0 and pos < len(self.options):
            self.options.insert(pos, opt)
        else:
            self.options.append(opt)
        if opt.has_key('separator'):
            return
        if 'longarg' not in opt.keys():
            raise exceptions.ValueError('Item longarg cannot be ignored in an option specification dictionary')
        # allow alphabet, number and underscore (_).
        if not opt['longarg'].strip('=').replace('_', '').isalnum() or not opt['longarg'][0].isalpha():
            raise exceptions.ValueError('Invalid option name %s' % opt['longarg'].strip('='))
        if 'default' not in opt.keys() and 'separator' not in opt.keys():
            raise exceptions.ValueError('A default value must be provided for all options')
        if opt.has_key('arg') and \
            opt['arg'].endswith(':') != opt['longarg'].endswith('='):
            raise exceptions.ValueError('Error: arg and longarg should both accept or not accept an argument')
        if opt.has_key('arg') and (len(opt['arg'].rstrip(':')) != 1 or not opt['arg'][0].isalpha()):
            raise exceptions.ValueError('Short arg should have one and only one alphabetic character.')
        if opt.has_key('arg') and sum([x.has_key('arg') and x['arg'][0] == opt['arg'][0] for x in self.options]) > 1:
            raise exceptions.ValueError("Duplicated short argument '%s'" % opt['arg'].rstrip(':'))
        if opt['longarg'].rstrip('=') in reserved_options:
            raise exceptions.ValueError("Option '--%s' is reserved. Please use another name." % opt['longarg'].rstrip('='))
        if (not opt['longarg'].endswith('=')) and (opt.has_key('chooseOneOf') or \
            opt.has_key('chooseFrom')):
            raise exceptions.ValueError('Directive chooseOneOf or chooseFrom can only be used for option that expects a value.');
        if opt['longarg'].rstrip('=') in methods:
            raise exceptions.ValueError("Option '%s' conflicts with the '%s' member function of the simuOpt class." % \
                (opt['longarg'].rstrip('='), (opt['longarg'].rstrip('='))))
        if opt['longarg'].rstrip('=') in self.__dict__.keys():
            raise exceptions.ValueError("Option '%s' conflicts with attribute '%s' of this simuOpt object." % \
                (opt['longarg'].rstrip('='), (opt['longarg'].rstrip('='))))
        if not opt['longarg'].endswith('=') and opt.has_key('allowedTypes') and type(True) not in opt['allowedTypes']:
            raise exceptions.ValueError("Boolean type (True/False) should be allowed in boolean option %s." % opt['longarg'])
        #
        opt['value'] = opt['default']
        opt['processed'] = False
        #
        name = opt['longarg'].rstrip('=')
        if self.dict.has_key(name):
            raise exceptions.ValueError('Option %s already exists.' % name)
        self.dict[name] = opt


    def saveConfig(self, file, params=[]):
        '''Write a configuration file to *file*. This file can be later read
        with command line option ``-c`` or ``--config``. All parameters with a
        ``label`` entry are saved unless a list of parameters are specified
        in *params*.
        '''
        cfg = open(file,'w')
        try:
            # sys.argv[0] might not exist
            print >> cfg, "# configuration file for program", sys.argv[0]
        except:
            pass
        print >> cfg, "# Configuration saved at at ", time.asctime()
        for opt in self.options:
            if opt.has_key('separator'):
                continue
            if len(params) > 0 and opt['longarg'].rstrip('=') not in params:
                continue
            # no label, and is not specified in params
            if not opt.has_key('label') and len(params) == 0:
                continue
            print >> cfg
            # write arg and long arg
            if opt.has_key('label'):
                print >> cfg, "# label:\t%s" % opt['label']
            if opt.has_key('arg'):
                if opt['arg'].endswith(':'):
                    print >> cfg, "# shortarg:\t-%s = %s" % (opt['arg'].rstrip(':'), _prettyString(opt['value']))
                else:
                    print >> cfg, "# shortarg:\t-%s" % opt['arg']
            # write description
            if opt.has_key('description'):
                desc = opt['description'].splitlines()
                print >> cfg, "# description:"
                for d in desc:
                    print >> cfg, "#\t", d.strip()
            # write out option value, try to make it python readable
            print >> cfg, "%s = %s" % (opt['longarg'].rstrip('='), _prettyString(opt['value'], quoted=True))
        print >> cfg, "\n\n#The same options can be given by command line options (subject to minor changes)"
        cmd = "#    --gui=False "
        # shorter version
        scmd = "#    --gui=False "
        for opt in self.options:
            if opt.has_key('separator') or not opt.has_key('label'):
                continue
            defaultVal = opt.has_key('useDefault') and opt['useDefault'] and str(opt['value']) == str(opt['default'])
            if opt['longarg'][-1] == '=':
                if ',' in str(opt['value']):    # has single quote
                    arg = " --" + opt['longarg'][0:-1] + '=' + _prettyString(opt['value'], quoted=True)
                    cmd += arg
                    if not defaultVal:
                        scmd += arg
                else:
                    arg = " --" + opt['longarg'][0:-1] + "=" + _prettyString(opt['value'], quoted=True)
                    cmd += arg
                    if not defaultVal:
                        scmd += arg
            elif opt['value']: # this option is True
                cmd += " --" + opt['longarg']
                if not defaultVal:
                    scmd += " --" + opt['longarg']
        print >> cfg, ' \\\n#    '.join(textwrap.wrap(cmd, break_long_words=False))
        # print out shorter version
        print >> cfg, "\n\n#Or a shorter version if default arguments are ignored"
        print >> cfg, ' \\\n#    '.join(textwrap.wrap(scmd, break_long_words=False))
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
                if len(params) > 0 and opt['longarg'].rstrip('=') not in params:
                    continue
                name = opt['longarg'].rstrip('=')
                scan = re.compile(name + r'\s*=\s*(.*)')
                if scan.match(line):
                    value = scan.match(line).groups()[0]
                    opt['value'] = _getParamValue(opt, value.strip('''"'\n'''))
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
            print self.usage()
            return False
        for opt in self.options:
            if opt.has_key('separator'):
                continue
            if len(params) > 0 and opt['longarg'].rstrip('=') not in params:
                continue
            if not opt['longarg'].endswith('='): # do not expect an argument, simple
                if '--' + opt['longarg'] in cmdArgs:
                    idx = cmdArgs.index('--'+opt['longarg'])
                elif opt.has_key('arg') and '-' + opt['arg'] in cmdArgs:
                    idx = cmdArgs.index('-' + opt['arg'])
                else:
                    continue
                if idx in self.processedArgs:
                    raise exceptions.ValueError("Parameter " + cmdArgs[idx] + " has been processed before.")
                self.processedArgs.append(idx)
                opt['value'] = True
                opt['processed'] = True
                continue
            # this is a more complicated case
            name = opt['longarg'].rstrip('=')
            hasArg = [x.startswith('--%s=' % name) or x == '--' + name for x in cmdArgs]
            if True in hasArg:
                idx = hasArg.index(True)
                # case 1: --arg something
                if cmdArgs[idx] == '--' + name:
                    if idx in self.processedArgs or idx+1 in self.processedArgs:
                        raise exceptions.ValueError("Parameter " + cmdArgs[idx] + " has been processed before.")
                    try:
                        val = _getParamValue(p, cmdArgs[idx+1])
                        self.processedArgs.extend([idx, idx+1])
                        opt['value'] = val
                        opt['processed'] = True
                    except Exception, e:
                        print "ERROR: Failed to assign parameter %s with value '%s'" % (opt['longarg'].rstrip('='),
                            cmdArgs[idx+1])
                        print e
                        continue
                # case 2 --arg=something
                else:
                    try:
                        val = _getParamValue(opt, cmdArgs[idx][(len(name)+3):])
                        self.processedArgs.append(idx)
                        opt['value'] = val
                        opt['processed'] = True
                    except Exception, e:
                        print "ERROR: Failed to assign parameter %s with value '%s'" % (opt['longarg'].rstrip('='),
                            cmdArgs[idx][(len(name)+3):])
                        print e
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
                        raise exceptions.ValueError("Parameter " + cmdArgs[idx] + " has been processed before.")
                    try:
                        val = _getParamValue(opt, cmdArgs[idx+1])
                        self.processedArgs.extend([idx, idx+1])
                        opt['value'] = val
                        opt['processed'] = True
                    except Exception, e:
                        print "ERROR: Failed to assign parameter %s with value '%s'" % (opt['longarg'].rstrip('='),
                            cmdArgs[idx+1])
                        print e
                        continue
                # case 2: -aopt or -a=opt
                else:
                    if idx in self.processedArgs:
                        raise exceptions.ValueError("Parameter " + cmdArgs[idx] + " has been processed before.")
                    try:
                        arg = cmdArgs[idx]
                        if len(arg) > 3 and arg[2] == '=':
                            val = _getParamValue(opt, cmdArgs[idx][3:])
                        else:
                            val = _getParamValue(opt, cmdArgs[idx][2:])
                        self.processedArgs.append(idx)
                        opt['value'] = val
                        opt['processed'] = True
                    except Exception, e:
                        print "ERROR: Failed to assign parameter %s with value '%s'" % (opt['longarg'].rstrip('='),
                            cmdArgs[idx])
                        print e
                        continue
        return True

    def usage(self, usage='usage: %prog [-opt [arg] | --opt [=arg]] ...'):
        '''Reutn the usage message from the option description list.
        ``'%prog'`` in parameter *usage* will be replaced by
        ``os.path.basename(sys.argv[0])``.
        '''
        return _usage(self.options, self.doc + '\n' + self.details)

    def termGetParam(self, params=[]):
        '''Get parameters from interactive user input. By default, all
        parameters are processed unless one of the following conditions is
        met:

        1. Parameter without a label
        2. Parameter with ``useDefault`` set to ``True``
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
                (opt.has_key('useDefault') and opt['useDefault'])):
                continue
            if len(params) > 0 and opt['longarg'].rstrip('=') not in params:
                continue
            # prompt
            if opt.has_key('label'):
                prompt = '%s (%s): ' % (opt['label'], str(opt['value']))
            else:
                prompt = '%s (%s): ' % (opt['longarg'].rstrip('='), str(opt['value']))
            while True:
                value = raw_input('\n' + prompt)
                if value == '':
                    # use existing value...
                    opt['processed'] = True
                    break
                try:
                    opt['value'] = _getParamValue(opt, value)
                    opt['processed'] = True
                    break
                except:
                    print "Invalid input.\n"
            if value is None and not (opt.has_key('allowedTypes') and types.NoneType in opt['allowedTypes']):
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
        '''
        title = os.path.split(sys.argv[0])[-1]
        if gui != 'Tkinter':
            try:
                return _wxParamDialog(self.options, title, self.doc, self.details, nCol).getParam()
            except exceptions.ImportError:
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
                raise exceptions.ValueError("Unprocessed command line argument: " + cmdArgs[i])
        #
        if self.gui == False:
            return self.termGetParam()
        # GUI
        try:
            return self.guiGetParam(nCol, gui=self.gui)
        except exceptions.ImportError:
            return self.termGetParam()
        return False

    def asDict(self):
        '''Return parameters as a dictionary.'''
        res = {}
        for opt in self.options:
            if opt.has_key('separator'):
                continue
            name = opt['longarg'].rstrip('=')
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

#
# simuOptions that will be checked by simuPOP.py when simuPOP is loaded.
# This structure can be changed by function setOptions
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


