#!/usr/bin/env python

############################################################################
#  Copyright (C) 2004 by Bo Peng
#  bpeng@rice.edu
#
#  $LastChangedDate$
#  $Rev$
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the
#  Free Software Foundation, Inc.,
#  59 Temple Place - Suite 330, Boston, MA    02111-1307, USA.
############################################################################

'''
Module  simuOpt  can be used to control which simuPOP module to load, and
how it is loaded using function  setOptions  . It also provides a simple
way to set simulation options, from user input, command line, configuration
file or a parameter dialog. All you need to do is to define an option
description list that lists all parameters in a given format, and call
the getParam function.

This module, if loaded, pre-process the command line options. More specifically,
it checks command line option:

-c configfile: read from a configuration file

--config configfile: the same as -c

--optimized: load optimized modules, unless setOption explicitly use non-optimized
  modules.

-q: Do not display banner information when simuPOP is loaded

--quiet: the same as -q

--useTkinter: force the use of Tcl/Tk dialog even when wxPython is available. By
  default, wxPython is used whenever possible.

--noDialog: do not use option dialog. If the options can not be obtained from
  command line or configuraiton file, users will be asked to input them interactively.


Because these options are reserved, you can not use them in your simuPOP script.
'''

# First try to get environmental variable

import os, sys, exceptions, types, re, time, imp, textwrap

allowed_keys = ['arg', 'longarg', 'label', 'allowedTypes', 'prompt', 'useDefault', 'jump', \
    'jumpIfFalse', 'default', 'description', 'validate', 'chooseOneOf', 'chooseFrom', 'separator']

allowed_commandline_options = ['-c', '--config', '--optimized', \
    '-q', '--useTkinter', '--quiet', '--noDialog']

def _getParamShortArg(p, processedArgs):
    ''' try to get a param from short arg '''
    if not p.has_key('arg'):
        return None
    if p['arg'] == 'c':
        raise exceptions.ValueError("-c option is reserved for config file.")
    if p['arg'][-1] == ':': # expecting an argument
        try:
            idx = map(lambda x:x[:2]=='-'+p['arg'][0], sys.argv[1:]).index(True)
            # has something like -a
            # case 1: -a file
            if sys.argv[idx+1] == '-'+p['arg'][0]:
                if idx+1 in processedArgs or idx+2 in processedArgs:
                    raise exceptions.ValueError("Parameter " + sys.argv[idx+1] + " has been processed before.")
                try:
                    val = _getParamValue(p, sys.argv[idx+2])
                    processedArgs.append(idx+1)
                    processedArgs.append(idx+2)
                    return val
                except:
                    return None
            # case 2: -aopt or -a=opt
            else:
                if idx+1 in processedArgs:
                    raise exceptions.ValueError("Parameter " + sys.argv[idx+1] + " has been processed before.")
                try:
                    opt = sys.argv[idx+1]
                    if len(opt) > 3 and opt[2] == '=':
                        val = _getParamValue(p, sys.argv[idx+1][3:])
                    else:
                        val = _getParamValue(p, sys.argv[idx+1][2:])
                    processedArgs.append(idx+1)
                    return val
                except:
                    return None
        except:
            # not available
            return None
    else:     # true or false
        # handle -h option, as a special case
        if '-'+p['arg'] in sys.argv[1:]:
            idx = sys.argv[1:].index('-'+p['arg'])
            if idx+1 in processedArgs:
                raise exceptions.ValueError("Parameter " + sys.argv[idx+1] + " has been processed before.")
            processedArgs.append(idx+1)
            return True
        else:
            return None


def _getParamLongArg(p, processedArgs):
    ''' get param from long arg '''
    if not p.has_key('longarg'):
        return None
    if p['longarg'] == 'config':
        raise exceptions.ValueError("--config option is reserved for config gile.")
    if p['longarg'][-1] == '=': # expecting an argument
        try:
            endChar = len(p['longarg'].split('=')[0])
            idx = map(lambda x:x[:(endChar+2)]=='--'+p['longarg'][0:endChar], sys.argv[1:]).index(True)
            # case 1: --arg something
            if sys.argv[idx+1] == '--'+p['longarg'][0:-1]:
                if idx+1 in processedArgs or idx+2 in processedArgs:
                    raise exceptions.ValueError("Parameter " + sys.argv[idx+1] + " has been processed before.")
                try:
                    val = _getParamValue(p, sys.argv[idx+2])
                    processedArgs.append(idx+1)
                    processedArgs.append(idx+2)
                    return val
                except:
                    return None
            # case 2 --arg=something
            else:
                if sys.argv[idx+1][endChar+2] != '=':
                    raise exceptions.ValueError("Parameter " + sys.argv[idx+1] + " is invalid. (--longarg=value)")
                try:
                    val = _getParamValue(p, sys.argv[idx+1][(endChar+3):])
                    processedArgs.append(idx+1)
                    return val
                except:
                    return None
        except:
            # not available
            return None
    else:     # true or false
        if '--'+p['longarg'] in sys.argv[1:]:
            idx = sys.argv[1:].index('--'+p['longarg'])
            if idx+1 in processedArgs:
                raise exceptions.ValueError("Parameter " + sys.argv[idx+1] + " has been processed before.")
            processedArgs.append(idx+1)
            return True


def _getParamConfigFile(p, processedArgs):
    ''' get param from configuration file    '''
    if not p.has_key('longarg'):
        return None
    try:         # check -c and --config
        idx = sys.argv[1:].index('-c')
        processedArgs.append(idx+1)
        processedArgs.append(idx+2)
        config = sys.argv[idx+2]
    except:
        try:
            idx = sys.argv[1:].index('--config')
            processedArgs.append(idx+1)
            processedArgs.append(idx+2)
            config = sys.argv[idx+2]
        except:
            return None
    # OK
    # read configuration file
    # deal with () in label.
    if p['longarg'][-1] == '=':
        name = p['longarg'][0:-1]
    else:
        name = p['longarg']
    scan = re.compile(name+r'\s*=\s*(.*)')
    try:
        file = open(config)
        for l in file.readlines():
            try:
                (value,) = scan.match(l).groups()
            except:
                # does not match
                continue
            else:
                file.close()
                try:
                    return _getParamValue(p, value.strip('''"'\n'''))
                except:
                    return None
        file.close()
        # get nothing
        return None
    except:    # can not open file
        print "Can not open configuration file ", config
        return None


def _getParamUserInput(p):
    ''' get param from user input '''
    # prompt
    if p.has_key('prompt'):
        prompt = p['prompt']
    elif p.has_key('label'):
        prompt = '%s (%s): ' % (p['label'], str(p['default']))
    elif p.has_key('longarg'):
        prompt = '--%s (%s): ' % (p['longarg'], str(p['default']))
    elif p.has_key('shortarg'):
        prompt = '-%s (%s): ' % (p['shortarg'], str(p['default']))
    else:
        raise exceptions.ValueError('Do not know how to prompt for user input (no label, longarg etc)')
    while True:
        value = raw_input('\n' + prompt)
        if value == '':
            value = None    # will use default value
            break
        else:
            try:
                return _getParamValue(p, value)
            except:
                print "Invalid input.\n"
                continue
    if value == None:
        if p.has_key('default'):
            return p['default']
        else:
            raise exceptions.ValueError("Can not get param for parameter (no default value): " + str(p['longarg']))


def _getParamValue(p, val):
    ''' try to get a value from value, raise exception if error happens. '''
    # if we are giving a unicode string, convert!
    if type(val) == types.UnicodeType:
        val = str(val)
    if (not p.has_key('allowedTypes')) or type(val) in p['allowedTypes']:
        if p.has_key('validate') and not p['validate'](val):
                raise exceptions.ValueError("Value "+str(val)+' does not pass validation')
        return val
    # handle another 'auto-boolean' case
    elif (p.has_key('arg') and p['arg'][-1] != ':') or \
        (p.has_key('longarg') and p['longarg'][-1] != '='):
        if val in ['1', 'true', 'True']:
            return True
        elif val in ['0', 'false', 'False']:
            return False
        else:
            raise exceptions.ValueError('Expect 0/1, true/false for boolean values for parameter %s ' % p['longarg'])
    # other wise, need conversion
    if type(val) in [types.StringType, types.UnicodeType] :
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
        raise ValueError('Type of input parameter ' + str(val) + " is incorrect. (param " \
            + p.setdefault('longarg','none') +")")


def _termGetParam(options, useDefault=False, checkUnprocessedArgs=False):
    ''' using user input to get param '''
    # get param from short arg
    processedArgs = []
    # process all options
    values = []
    goto = 0
    for opt in range(0, len(options)):
        p = options[opt]
        # validate p
        for k in p.keys():
            if not k in allowed_keys:
                raise exceptions.ValueError("Unrecognized option entry " + k )
        if p.has_key('separator'):
            continue
        val = _getParamShortArg(p, processedArgs)
        if val == None:
            val = _getParamLongArg(p, processedArgs)
        if val == None:
            val = _getParamConfigFile(p, processedArgs)
        if val == None:
            if (useDefault or (not p.has_key('label')) or (p.has_key('useDefault') and p['useDefault'])) and p.has_key('default'):
                val = p['default']
            elif opt >= goto:
                val = _getParamUserInput(p)
        # these parameters are skipped, but still processed to check unprocessed args
        if opt < goto:
            values.append(val)
        elif val == None:
            # should have a valid value now.
            raise exceptions.ValueError("Failed to get parameter " + p.setdefault("label",'') + " " + p.setdefault("longarg",''))
        else:
            values.append( _getParamValue(p, val))
        # now we really short have something not None, unless the default is None
        # if a string is fine
        # now, deal with jump option
        if (values[-1] == True and p.has_key('jump')) or \
            (values[-1] == False and p.has_key('jumpIfFalse')):
            if p.has_key('jump'):
                jumpTo = p['jump']
            else:
                jumpTo = p['jumpIfFalse']
            if jumpTo in [-1, None, '']:    # go to last
                goto = len(options)
            elif type(jumpTo) == type(''):  # go to another parameter
                goto = -1
                for s_idx in range(opt + 1, len(options)):
                    s_arg = options[s_idx]
                    if s_arg.has_key('longarg') and \
                        ((s_arg['longarg'][-1] == '=' and s_arg['longarg'][:-1] == jumpTo) or \
                         (s_arg['longarg'][-1] != '=' and s_arg['longarg'] == jumpTo)):
                        goto = s_idx;
                        break
                if goto == -1:
                    raise ValueError('Failed to jump to option %s.' % jumpTo)
            elif jumpTo <= opt:
                raise ValueError("Can not stay or jump backwards when processing options.")
            else:
                goto = jumpTo
    # look if any argument was not processed
    if checkUnprocessedArgs:
        for i in range(1, len(sys.argv)):
            if (not sys.argv[i] in allowed_commandline_options) and (not i in processedArgs):
                raise exceptions.ValueError("Unprocessed command line argument: " + sys.argv[i])
    return values


class _paramDialog:
    def __init__(self, options, title = '', description='', details='', nCol=1):
        if len(options) == 0:
            raise exceptions.ValueError("Empty field names...")    # some behaviors
        # values, not the final result
        # first set them with command line options etc
        self.values = []
        #
        processedArgs = []
        for opt in options:
            # validate opt
            for k in opt.keys():
                if not k in allowed_keys:
                    raise exceptions.ValueError("Unrecognized option entry " + k )
            #
            if opt.has_key('separator'):
                val = opt['separator']
            else:
                val = _getParamShortArg(opt, processedArgs)
                if val == None:
                    val = _getParamLongArg(opt, processedArgs)
                if val == None:
                    val = _getParamConfigFile(opt, processedArgs)
                if val == None:
                    if opt.has_key('default'):
                        val = opt['default']
            self.values.append(val)
        # look if any argument was not processed
        for i in range(1, len(sys.argv)):
            if (not sys.argv[i] in allowed_commandline_options) and (not i in processedArgs):
                raise exceptions.ValueError("Unprocessed command line argument: " + sys.argv[i])
        # now, initialize variables
        self.options = options
        self.title = title
        self.description = description
        self.details = details
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

    def formatDesc(self, text):
        # linux can auto wrap, windows can not but sometime wrap
        # at unexpected places... It is safer to wrap at original
        # place.
        return '\n'.join([x.strip() for x in text.splitlines()] )

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
        if len(self.values) == len(self.options):
            # remove values inserted by separators
            ret = []
            for p,opt in enumerate(self.options):
                if not opt.has_key('separator'):
                    ret.append(self.values[p])
            return ret
        else:
            return []


class _tkParamDialog(_paramDialog):
    def __init__(self, options, title = '', description='', details='', nCol=1):
        ''' get options from a given options structure '''
        _paramDialog.__init__(self, options, title, description, details, nCol)

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
        msg.insert(tk.END, usage(self.options, self.details))
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
        self.values = []
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
                            items.append(self.entryWidgets[g].get( s))
                    val = _getParamValue(self.options[g], items)
                elif self.entryWidgets[g].winfo_class() == "Checkbutton":    # a checkbutton (true or false)
                    # gets 0/1 for false/true
                    var = self.values[g].get()
                    val = _getParamValue(self.options[g], var == 1)
            except Exception,e:
                print e
                for lab in self.labelWidgets:
                    if lab is not None:
                        lab.configure(fg='black')
                # set this one to red
                self.labelWidgets[g].configure(fg='red')
                self.entryWidgets[g].focus_force()
                return
            else:
                # convert to values
                self.values[g] = val
        # get all results and return
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
            value = self.values[g]
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
                iv.set(self.values[g] == True) # value can be None, True or False
                self.values[g] = iv
                self.entryWidgets[g] = tk.Checkbutton(self.app, height=1,
                    text = "Yes / No", variable=self.values[g])
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
                    self.entryWidgets[g].insert(0, prettyOutput(value))
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

    def onHelp(self, event):
        # open another dialog
        helpDlg = wx.Dialog(parent=self.dlg, id=-1, title='Help for ' + self.title)
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(wx.TextCtrl(parent=helpDlg, id=-1, size=[600,400],
            style=wx.TE_MULTILINE | wx.TE_READONLY,
            value=usage(self.options, self.details)), 0, wx.ALL, 20)
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
                self.values[g] = val
        # get all results and return
        self.dlg.EndModal(wx.ID_OK)

    def onCancel(self, event):
        '''When ESC is pressed cancel, clear values and return'''
        self.values = []
        self.dlg.EndModal(wx.ID_CANCEL)

    def addButton(self, ID, text, func, parent, box):
        button = wx.Button(parent, ID, text)
        box.Add(button, 0, wx.ALIGN_CENTER)
        self.dlg.Bind(wx.EVT_BUTTON, func, button)

    def createDialog(self):
        self.app = wx.App(0)
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
            if not (opt.has_key('label') or opt.has_key('separator')) :
                continue
            value = self.values[g]
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
            # label
            self.labelWidgets[g] = wx.StaticText(parent=self.dlg, id=-1, label=opt['label'])
            gridBox[colIndex].Add(self.labelWidgets[g], 0, wx.ALIGN_LEFT )
            # use different entry method for different types
            if opt.has_key('chooseOneOf'):    # single choice
                self.entryWidgets[g] = wx.Choice(parent=self.dlg, id=g, choices = opt['chooseOneOf'])
                if opt.has_key('description'):
                    self.entryWidgets[g].SetToolTipString(self.formatDesc(opt['description']))
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
                    self.entryWidgets[g].SetToolTipString(self.formatDesc(opt['description']))
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
                    self.entryWidgets[g].SetToolTipString(self.formatDesc(opt['description']))
                if value is not None:
                    self.entryWidgets[g].SetValue(value)
                gridBox[colIndex].Add(self.entryWidgets[g], 1, wx.EXPAND)
                rowIndex += 1
            else: # an edit box
                # put default value into the entryWidget
                txt = ''
                if value is not None:
                    txt = prettyOutput(value)
                self.entryWidgets[g] = wx.TextCtrl(parent=self.dlg, id=g, value=txt)
                if opt.has_key('description'):
                    self.entryWidgets[g].SetToolTipString(self.formatDesc(opt['description']))
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


# get parameter
def getParam(options=[], doc="", details="", noDialog=False, UnprocessedArgs=True, verbose=False, nCol=1):
    """ Get parameters from either
            - a Tcl/Tk based, or wxPython based parameter dialog
              (wxPython is used if it is available)
            - command line argument
            - configuration file specified by  -c file   (  --config  file), or
            - prompt for user input

        The option description list consists of dictionaries with some
        predefined keys. Each dictionary defines an option. Each option
        description item can have the following keys:

        arg: short command line option name.  'h'  checks the presence of argument  -h  .
        If an argument is expected, add a comma to the option name. For example,  'p:'
        matches command line option  -p=100  or  -p 100  .

        longarg: long command line option name.  'help'  checks the presence of
            argument  '--help'  .  'mu='  matches command line
            option  --mu=0.001  or  -mu 0.001  .

        label: The label of the input field in a parameter dialog, and as the prompt for
          user input.

        default: default value for this parameter. It is used to as the default value
          in the parameter dialog, and as the option value when a user presses  'Enter'
          directly during interactive parameter input.

        useDefault: use default value without asking, if the value can not be determined
          from GUI, command line option or config file. This is useful for options that
          rarely need to be changed. Setting them to useDfault allows shorter command
          lines, and easy user input.

        description: a long description of this parameter, will be put into the usage
          information, which will be displayed with (  -h  ,  --help  command line option, or
          help button in parameter dialog).

        allowedTypes: acceptable types of this option. If  allowedTypes  is  types.ListType
          or  types.TupleType  and the user's input is a scalar, the input will be converted
          to a list automatically. If the conversion can not be done, this option will
          not be accepted.

        validate: a function to validate the parameter. You can define your own functions
          or use the ones defined in this module.

        chooseOneOf: if specified,  simuOpt  will choose one from a list of values using a
          listbox (Tk) or a combo box (wxPython) .

        chooseFrom: if specified,  simuOpt  will choose one or more items from a list of
          values using a listbox (tk) or a combo box (wxPython).

        separator: if specified, a blue label will be used to separate groups of
          parameters.

        jump: it is used to skip some parameters when doing the interactive user input.
          For example,  getParam  will skip the rest of the parameters if  -h  is specified
          if parameter  -h  has item  'jump':-1  which means jumping to the end.
          Another situation of using this value is when you have a hierarchical parameter
          set. For example, if mutation is on, specify mutation rate, otherwise proceed.
          The value of this option can be the absolute index or the  longarg  name of
          another option.

        jumpIfFalse: The same as jump but jump if current parameter is  False  .


        This function will first check command line argument. If the argument
        is available, use its value. Otherwise check if a config file is
        specified. If so, get the value from the config file. If both failed,
        prompt user to input a value. All input will be checked against types,
        if exists, an array of allowed types.

        Parameters of this function are:

        options: a list of option description dictionaries

        doc: short description put to the top of parameter dialog

        details: module help. Usually set to  __doc__  .

        noDialog: do not use a parameter dialog, used in batch mode. Default to False.

        checkUnprocessedArgs: obsolete because unused args are always checked.

        verbose: whether or not print detailed info

        nCol: number of columns in the parameter dialog.
    """
    # check if --noDialog, -h is present
    # or there is no 'label' in the options structure
    # for backward compatibility, change 'configName' to 'label'
    for opt in options:
        if opt.has_key('configName'):
            print 'Warning: configName is obsolete, please use "label" instead'
            opt['label'] = opt['configName']
        if not opt.has_key('default') and not opt.has_key('separator'):
            raise exceptions.ValueError('Error: a default value must be provided for all options')
        if opt.has_key('arg') and opt.has_key('longarg') and\
            opt['arg'].endswith(':') != opt['longarg'].endswith('='):
            raise exceptions.ValueError('Error: arg and longarg should both accept or not accept an argument')

    if noDialog or par_noDialog or '-h' in sys.argv[1:] or '--help' in sys.argv[1:] \
        or True not in map(lambda x:x.has_key('label'), options):
        return _termGetParam(options, False, True)
    else:
        title = os.path.split(sys.argv[0])[-1]
        if useTkinter:
            return _tkParamDialog(options, title, doc, details, nCol).getParam()
        elif useWxPython:
            return _wxParamDialog(options, title, doc, details, nCol).getParam()
        else:
            return _termGetParam(options, False, True)


def usage(options, before=''):
    """ Print usage information from the option description list. Used
    with  -h  (or  --help   ) option, and in the parameter input dialog.

    options: option description list.

    before: optional information
    """
    message = ''
    if before != '':
        message += '    ' + before + '\n'
    message += '\n' + sys.argv[0] + ' usage:\n'
    message += '    > ' + sys.argv[0] + ' options\n\n'
    message += '    Options: (-shortoption --longoption: description.)\n'
    message += '        -c xxx --config xxx :\n                Load parameters from file xxx\n'
    message += '        --noDialog :\n                Enter parameter from command line\n'
    message += '        --optimized :\n                Use optimized library (no error checking)\n'
    for p in options:
        message += "        "
        if p.has_key('arg'):
            if p['arg'][-1] == ':':
                message += '-'+ p['arg'][0:-1] + ' xxx '
            else:
                message += '-'+ p['arg'] + ' '
        if p.has_key('longarg'):
            if p['longarg'][-1] == '=':
                message += '--' + p['longarg'][0:-1] + ' xxx '
            else:
                message += '--' + p['longarg'] + ' '
        if p.has_key('label'):
            message += '(config file entry: ' + p['label'] + ')'
        message += ':\n                '
        if p.has_key('description'):
            message += p['description']
        message += '\n'
        if p.has_key('default') and p['default'] is not None:
            message +=    '                Default to ' + prettyOutput(p['default']) + '\n'
            message += '\n'
    return message


def prettyOutput(value, quoted=False, outer=True):
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
        txt = '[' + ', '.join([prettyOutput(x, True, False) for x in value]) + ']'
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


def saveConfig(opt, file, param):
    """ Write a configuration file. This file can be later read with
    command line option  -c  or --config  .

    opt: the option description list

    file: output file

    param: parameters returned from  getParam

    """
    try:
        f = open(file,'w')
    except:
        print 'Can not open ', file , ' to write.'
        return
    options = []
    for g in opt:
        # execlude separators from opt
        if not g.has_key('separator'):
            options.append(g)
    if len(options) != len(param):
        raise ValueError("Length of option specification and param should be the same.")
    print >> f, "# configuration file for program", sys.argv[0]
    print >> f, "# saved at ", time.asctime()
    for p in range(0, len(options)):
        if options[p].has_key('label'):
            print >> f
            # write arg and long arg
            if options[p].has_key('label'):
                print >> f, "# label:\t%s" % options[p]['label']
            if options[p].has_key('arg'):
                if options[p]['arg'][-1] == ':':
                    print >> f, "# shortarg:\t-%s = %s" % (options[p]['arg'][:-1], prettyOutput(param[p]))
                else:
                    print >> f, "# shortarg:\t-%s" % options[p]['arg']
            # write description
            if options[p].has_key('description'):
                desc = options[p]['description'].splitlines()
                print >> f, "# description:"
                for d in desc:
                    print >> f, "#\t", d.strip()
            if options[p].has_key('longarg'):
                if options[p]['longarg'][-1] == '=':
                    arg = options[p]['longarg'][:-1]
                else:
                    arg = options[p]['longarg']
                # write out option value, try to make it python readable
                print >> f, "%s = %s" % (arg, prettyOutput(param[p], quoted=True))
    print >> f, "\n\n#The same options can be given by command line options (subject to minor changes)"
    cmd = "#    --noDialog "
    # shorter version
    scmd = "#    --noDialog "
    for p in range(0, len(options)):
        if options[p].has_key('label') and options[p].has_key('longarg'):
            defaultVal = options[p].has_key('useDefault') and options[p]['useDefault'] \
                and str(param[p]) == str(options[p]['default'])
            if options[p]['longarg'][-1] == '=':
                if str(param[p]).find(",") >= 0:    # has single quote
                    arg = " --" + options[p]['longarg'][0:-1] \
                        + '=' + prettyOutput(param[p], quoted=True)
                    cmd += arg
                    if not defaultVal:
                        scmd += arg
                else:
                    arg = " --" + options[p]['longarg'][0:-1] \
                        + "=" + prettyOutput(param[p], quoted=True)
                    cmd += arg
                    if not defaultVal:
                        scmd += arg
            elif param[p]: # this option is True
                cmd += " --" + options[p]['longarg']
                if not defaultVal:
                    scmd += " --" + options[p]['longarg']
    print >> f, ' \\\n#    '.join(textwrap.wrap(cmd, break_long_words=False))
    # print out shorter version
    print >> f, "\n\n#Or a shorter version if default arguments are ignored"
    print >> f, ' \\\n#    '.join(textwrap.wrap(scmd, break_long_words=False))
    f.close()


def printConfig(opt, param, out=sys.stdout):
    """ Print configuration.

        opt: option description list

        param: parameters returned from  getParam()

        out: output
    """
    # remove separators from opt
    options = []
    for g in opt:
        if not g.has_key('separator'):
            options.append(g)
    if len(options) != len(param):
        raise ValueError("Length of option specification and param should be the same.")
    for p in range(0, len(options)):
        if options[p].has_key('label'):
            if type(param[p]) == types.StringType:
                print >> out, options[p]['label'], '\t"'+str(param[p])+'"'
            else:
                print >> out, options[p]['label'], '\t', str(param[p])


# define some validataion functions
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


env_optimized = os.getenv('SIMUOPTIMIZED')
env_longAllele = os.getenv('SIMUALLELETYPE')
env_debug = os.getenv('SIMUDEBUG')

[par_optimized] = _termGetParam([{'longarg':'optimized', \
    'default':''}], True, False)
[par_quiet] = _termGetParam([{'arg':'q','longarg':'quiet', \
    'default':False}], True, False)
[par_useTkinter] = _termGetParam([{'longarg':'useTkinter', \
    'default':False }], True, False)
[par_noDialog] = _termGetParam([{'longarg':'noDialog', \
    'default':False }], True, False)

# remove these parameters from sys.argv
for arg in ['--optimized', '--quiet', '-q', '--useTkinter']:
    try:
        sys.argv.remove(arg)
    except:
        pass

if par_optimized != '':
    _optimized = par_optimized
elif env_optimized is not None:
    _optimized = True
else:     # default to false
    _optimized = False

if env_longAllele in ['standard', 'short', 'long', 'binary']:
    _longAllele = env_longAllele
else:
    _longAllele = 'standard'

simuOptions = {'Optimized':_optimized,
    'AlleleType':_longAllele, 'Debug':[], 'Quiet':par_quiet}

if env_debug is not None:
    simuOptions['Debug'].extend( env_debug.split(',') )

def setOptions(optimized=None, mpi=None, chromMap=[], alleleType=None, quiet=None, debug=[]):
    '''set options before simuPOP is loaded to control which simuPOP module to load,
    and how the module should be loaded.

    optimized: whether or not load optimized version of a module. If not set,
        environmental variable SIMUOPTIMIZED, and commandline option --optimized
        will be used if available. If nothing is defined, standard version will
        be used.

    mpi: obsolete.

    chromMap: obsolete.

    alleleType: 'binary', 'short', or 'long'. 'standard' can be used as 'short'
        for backward compatibility. If not set, environmental variable
        SIMUALLELETYPE will be used if available. if it is not defined, the
        short allele version will be used.

    quiet: If True, supress banner information when simuPOP is loaded.

    debug: a list of debug code (or string). If not set, environmental variable
        SIMUDEBUG will be used if available.

    '''
    if optimized in [True, False]:
        simuOptions['Optimized'] = optimized
    if alleleType in ['standard', 'long', 'binary', 'short']:
        simuOptions['AlleleType'] = alleleType
    if quiet in [True, False]:
        simuOptions['Quiet'] = quiet
    if len(debug) > 0:
        simuOptions['Debug'].extend(debug)

# short = standard
if simuOptions['AlleleType'] == 'standard':
    simuOptions['AlleleType'] = 'short'
if simuOptions['Optimized'] not in [True, False]:
    simuOptions['Optimized'] = False

def requireRevision(rev):
    '''Compare the revision of this simuPOP module with given revision. Raise
    an exception if current module is out of date. '''
    if simuRev() <= rev:
        raise exceptions.SystemError('''This script requires simuPOP revision >= %d,
            please upgrade your installation. ''' % rev)

useTkinter = False
useWxPython = False

if not par_useTkinter:
    try:
        # wxPython might not exist
        imp.find_module('wx')
        if not par_noDialog:
            import wx
    except:
        useWxPython = False
    else:
        useWxPython = True


if par_useTkinter or not useWxPython:
    # Tkinter should almost always exists, but ...
    try:
        imp.find_module('Tkinter')
        if not par_noDialog:
            import Tkinter as tk
            import tkFont as tkFont
    except:
        print "Tkinter can not be loaded. Please check your Python installation."
        useTkinter = False
    else:
        useTkinter = True
        # this is not possible now because of the use of find_module
        #if TkVersion < 8.0 :
        #    useTkinter = False

