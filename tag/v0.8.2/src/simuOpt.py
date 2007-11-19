#!/usr/bin/env python

############################################################################
#        Copyright (C) 2004 by Bo Peng
#        bpeng@rice.edu
#
#        $LastChangedDate$
#        $Rev$
#
#        This program is free software; you can redistribute it and/or modify
#        it under the terms of the GNU General Public License as published by
#        the Free Software Foundation; either version 2 of the License, or
#        (at your option) any later version.
#
#        This program is distributed in the hope that it will be useful,
#        but WITHOUT ANY WARRANTY; without even the implied warranty of
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
#        GNU General Public License for more details.
#
#        You should have received a copy of the GNU General Public License
#        along with this program; if not, write to the
#        Free Software Foundation, Inc.,
#        59 Temple Place - Suite 330, Boston, MA    02111-1307, USA.
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

import os, sys, exceptions, types, re, time, imp

allowed_keys = ['arg', 'longarg', 'label', 'allowedTypes', 'prompt', 'useDefault', 'jump', \
    'jumpIfFalse', 'default', 'description', 'validate', 'chooseOneOf', 'chooseFrom', 'separator']

allowed_commandline_options = ['-c', '--config', '--optimized', '--mpi', '--chromMap', \
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


def _termGetParam(options, checkUnprocessedArgs=True, verbose=False, useDefault=False):
    ''' using user input to get param '''
    # get param from short arg
    processedArgs = []
    # process all options
    values = []
    goto = 0
    for opt in range(0, len(options)):
        if opt < goto:
            values.append(None)
            continue
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
            else:
                val = _getParamUserInput(p)
        # should have a valid value now.
        if val == None:
            raise exceptions.ValueError("Failed to get parameter " + p.setdefault("label",'') + " " + p.setdefault("longarg",''))
        values.append( _getParamValue(p, val))
        # now we really short have something not None, unless the default is None
        # if a string is fine
        # now, deal with jump option
        if values[-1] == True and p.has_key('jump'):
            if p['jump'] == -1:    # go to last
                goto = len(options)
            else:
                goto = p['jump']
            # can not jump backward
            if goto <= opt:
                raise ValueError("Can not stay or jump backwards when processing options.")
        if values[-1] == False and p.has_key('jumpIfFalse'):
            if p['jumpIfFalse'] == -1:    # go to last
                goto = len(options)
            else:
                goto = p['jumpIfFalse']
            # can not jump backward
            if goto <= opt:
                raise ValueError("Can not stay or jump backwards when processing options.")
    # look if any argument was not processed
    if checkUnprocessedArgs and (True not in map(lambda x:x.has_key('jump'), options) \
        and True not in map(lambda x:x.has_key('jumpIfFalse'), option)):
        for i in range(1, len(sys.argv)):
            if (not sys.argv[i] in allowed_commandline_options) and (not i in processedArgs):
                raise exceptions.ValueError("Unprocessed command line argument: " + sys.argv[i])
    return values


###
### currently, I do not want to specify font. Use defaults on each
### OS seems to be a better choice.
#DEFAULT_FONT_FAMILY     = ("MS", "Sans", "Serif")
#MONOSPACE_FONT_FAMILY = ("Courier")
#DEFAULT_FONT_SIZE         = 12
#BIG_FONT_SIZE = 12

#
# This function is adapted from easyGUI.
# It is easy dirty as it can be since I need a
# function to handle all the stuff.
# Maybe a class approach will be used later.
#
def _tkGetParam(opt, title = '', description='', details='', checkUnprocessedArgs=True, nCol=1):
    ''' get options from a given options structure '''
    import Tkinter as tk
    if len(opt) == 0:
        raise exceptions.ValueError("Empty field names...")    # some behaviors
    # remove separators, tk version does not do this
    options = []
    for g in opt:
        if not g.has_key('separator'):
            options.append(g)
    # values, not the final result
    # first set them with command line options etc
    values = []
    processedArgs = []
    # process all options
    goto = 0
    for opt in range(0, len(options)):
        if opt < goto:
            values.append(None)
            continue
        p = options[opt]
        # validate p
        for k in p.keys():
            if not k in allowed_keys:
                raise exceptions.ValueError("Unrecognized option entry " + k )
        val = _getParamShortArg(p, processedArgs)
        if val == None:
            val = _getParamLongArg(p, processedArgs)
        if val == None:
            val = _getParamConfigFile(p, processedArgs)
        if val == None:
            if p.has_key('default'):
                val = p['default']
        # ignore jump options
        values.append(val)
    # look if any argument was not processed
    if checkUnprocessedArgs:
        for i in range(1, len(sys.argv)):
            if (not sys.argv[i] in allowed_commandline_options) and (not i in processedArgs):
                raise exceptions.ValueError("Unprocessed command line argument: " + sys.argv[i])
    #
    def denyWindowManagerClose():
        """ don't allow WindowManager close	"""
        x = tk.Tk()
        x.withdraw()
        x.bell()
        x.destroy()
    root = tk.Tk()
    root.protocol('WM_DELETE_WINDOW', denyWindowManagerClose )
    entryWidgets = [None]*len(options)
    labelWidgets = [None]*len(options)
    root.title(title)
    root.iconname('Dialog')
    root.geometry("+300+200")
    # cancel
    def doCancel(event):
        " when ESC is pressed cancel "
        # scripple values by changing its length
        values.append(None)
        root.quit()
    # done
    def doGetText(event):
        ''' get result and convert values '''
        for g in range(len(entryWidgets)):
            if entryWidgets[g] == None:
                continue
            try:
                # get text from different type of entries
                if entryWidgets[g].winfo_class() == "Entry":    # an entry box?
                    val = _getParamValue( options[g], entryWidgets[g].get())
                elif entryWidgets[g].winfo_class() == "Listbox":    # a listbox
                    sel = entryWidgets[g].curselection()
                    if len(sel) == 1:
                        items = entryWidgets[g].get( sel)
                    else:
                        items = []
                        for s in sel:
                            items.append( entryWidgets[g].get( s))
                    val = _getParamValue( options[g], items)
                elif entryWidgets[g].winfo_class() == "Checkbutton":    # a checkbutton (true or false)
                    var = values[g].get()
                    val = _getParamValue( options[g], var)
            except:                
                #print "Invalid Value: ", entryWidgets[g].class()
                # incorrect value
                # set to red
                # clear other red colors
                for lab in labelWidgets:
                    if lab != None:
                        lab.configure(fg='black')
                # set this one to red
                labelWidgets[g].configure(fg='red')
                entryWidgets[g].focus_force()
                return
            else:
                # convert to values
                values[g] = val
        # get all results and return
        root.quit()
    # help
    def doHelp(event):
        # open another window
        root1 = tk.Tk()
        # OK for help
        def doOK(event):
            " OK buton is pressed "
            root1.quit()
        #global root1
        root1.title("Help for " + title)
        root1.iconname('Dialog')
        root1.geometry("+200+200")
        root1.minsize(400, 100)
        messageFrame = tk.Frame(root1)
        messageFrame.pack(side=tk.TOP, fill=tk.BOTH)
        scrollBar = tk.Scrollbar(messageFrame)
        scrollBar.pack(side=tk.RIGHT, fill=tk.Y)
        messageWidget = tk.Text(messageFrame, wrap=tk.WORD,
            yscrollcommand=scrollBar.set)
        messageWidget.insert(tk.END, usage(options, details))
        scrollBar.config(command=messageWidget.yview)
        #messageWidget.configure(font=(DEFAULT_FONT_FAMILY,DEFAULT_FONT_SIZE), state=DISABLED)
        messageWidget.pack(side=tk.TOP, expand=tk.YES, fill=tk.X, padx='3m', pady='3m')
        buttonFrame = tk.Frame(root1)
        buttonFrame.pack(side=tk.BOTTOM, fill=tk.BOTH)
        okButton = tk.Button(buttonFrame, takefocus=1, text="OK")
        okButton.pack(expand=tk.YES, padx='1m', pady='1m', ipadx='2m', ipady='1m')
        # bind the keyboard events to the widget
        okButton.bind("<Return>", doOK)
        okButton.bind("<Button-1>", doOK)
        root1.bind("<Escape>", doOK)
        # put the focus on the first button
        root1.mainloop()
        root1.destroy()
    #
    # the main window
    #
    root.bind("<Escape>", doCancel)
    # all use grid management
    # top message
    # do not use a very long description please
    tk.Message(root, text=description, width=600).grid(row=0, column=0,
        columnspan = 2 * nCol, sticky=tk.N+tk.E+tk.S+tk.W, pady=20)
    # find out number of items etc
    colParam = 0
    for opt in options:
        if opt.has_key('label'):
            colParam += 1
        if opt.has_key('chooseFrom'):
            colParam += len( opt['chooseFrom']) -1
        if opt.has_key('chooseOneOf'):
            colParam += len( opt['chooseOneOf']) -1
    if colParam / nCol * nCol == colParam:
        colParam /= nCol
    else:
        colParam = colParam/nCol + 1
    colCount = 0
    colIndex = 0
    # all entries
    for g in range(len(options)):
        opt = options[g]
        if not opt.has_key('label'):
            continue
        # --------- entryWidget ----------------------------------------------
        # use different entry method for different types
        if opt.has_key('chooseOneOf'):    # single choice
            labelWidgets[g] = tk.Label(root, text=opt['label'])
            labelWidgets[g].grid(column=colIndex*2, row=colCount%colParam+1, padx=10,
                rowspan = len(opt['chooseOneOf']), sticky=tk.E)
            entryWidgets[g] = tk.Listbox(root, width=40, selectmode=tk.SINGLE, \
                exportselection=0, height=len(opt['chooseOneOf']))
            entryWidgets[g].grid(column=colIndex*2+1, row=colCount%colParam+1, padx=10,
                rowspan = len(opt['chooseOneOf']))
            colCount += len(opt['chooseOneOf'])
            for entry in opt['chooseOneOf']:
                entryWidgets[g].insert(tk.END, str(entry))
            if values[g] != None:
                entryWidgets[g].select_set( opt['chooseOneOf'].index(values[g]))
        elif opt.has_key('chooseFrom'):    # multiple choice
            labelWidgets[g] = tk.Label(root, text=opt['label'])
            labelWidgets[g].grid(column=colIndex*2, row=colCount%colParam+1, padx=10,
                rowspan = len(opt['chooseFrom']), sticky=tk.E)
            entryWidgets[g] = tk.Listbox(root, width=40, selectmode=tk.EXTENDED, \
                exportselection=0, height=len( opt['chooseFrom']))
            entryWidgets[g].grid(column=colIndex*2+1, row=colCount%colParam+1, padx=10,
                rowspan = len(opt['chooseFrom']))
            colCount += len(opt['chooseFrom'])
            for entry in opt['chooseFrom']:
                entryWidgets[g].insert(tk.END, str(entry))
            if values[g] != None:
                if type(values[g]) in [types.TupleType, types.ListType]:
                    for val in values[g]:
                        entryWidgets[g].select_set( opt['chooseFrom'].index(val))
                else:
                    entryWidgets[g].select_set( opt['chooseFrom'].index( values[g] ))
        elif (opt.has_key('arg') and opt['arg'][-1] != ':') or \
             (opt.has_key('longarg') and opt['longarg'][-1] != '='):  # true or false
            labelWidgets[g] = tk.Label(root, text=opt['label'])
            labelWidgets[g].grid(column=colIndex*2, row=colCount%colParam+1, padx=10,
                rowspan = 1, sticky=tk.E)
            # replace values[g] by a tk IntVar() because tk.Checkbutton has to store
            # its value in such a variable. values[g].get() will be used to return the
            # state of this Checkbutton.
            # c.f. http://infohost.nmt.edu/tcc/help/pubs/tkinter/control-variables.html
            iv = tk.IntVar()
            iv.set(values[g] == True) # values[g] can be None, True or False
            values[g] = iv
            entryWidgets[g] = tk.Checkbutton(root, height=1,
                     text = "Yes / No", variable=values[g])
            entryWidgets[g].grid(column=colIndex*2+1, row=colCount%colParam+1, padx=10,
                rowspan = 1)
            colCount += 1
            entryWidgets[g].deselect()
        else:
            labelWidgets[g] = tk.Label(root, text=opt['label'])
            labelWidgets[g].grid(column=colIndex*2, row=colCount%colParam+1, padx=10, sticky=tk.E)
            entryWidgets[g] = tk.Entry(root, width=40)
            entryWidgets[g].grid(column=colIndex*2+1, row=colCount%colParam+1, padx=10)
            colCount += 1
             # put default value into the entryWidget
            if values[g] != None:
                # len()>0 to avoid emtpy string with emtpy list
                if type(values[g]) in [types.ListType, types.TupleType] and len(values[g])>1:
                    entryWidgets[g].insert(0, ', '.join(map(str, values[g])))
                else:
                    entryWidgets[g].insert(0,str(values[g]))
        colIndex = colCount /colParam
        entryWidgets[g].bind("<Return>", doGetText)
        entryWidgets[g].bind("<Escape>", doCancel)
    # help button
    helpButton = tk.Button(root, takefocus=1, text="Help")
    helpButton.bind("<Return>"    , doHelp)
    helpButton.bind("<Button-1>", doHelp)
    helpButton.grid( column=0, columnspan=nCol, row = colParam+1, pady=20)
    # ok button
    okButton = tk.Button(root, takefocus=1, text="Run!")
    okButton.bind("<Return>"    , doGetText)
    okButton.bind("<Button-1>", doGetText)
    okButton.grid( column=nCol, columnspan=nCol, row = colParam+1, pady=20)
    # cancel button
    cancelButton = tk.Button(root, takefocus=1, text="Cancel")
    cancelButton.bind("<Return>"    , doCancel)
    cancelButton.bind("<Button-1>", doCancel)
    cancelButton.grid( column=0, columnspan=2*nCol, row = colParam+1, pady=20)
    # ------------------- time for action! -----------------
    # first un-none
    for g in range(len(options)):
        if entryWidgets[g] != None:
            entryWidgets[g].focus_force()
            break
    root.mainloop()    # run it!
    # -------- after the run has completed ----------------------------------
    root.destroy()    # button_click didn't destroy root, so we do it now
    if len(values) == len(options):
        return values
    else:
        return []

#
# This function is adapted from easyGUI.
# It is easy dirty as it can be since I need a
# function to handle all the stuff.
# Maybe a class approach will be used later.
#
def _wxGetParam(options, title = '', description='', details='', checkUnprocessedArgs=True, nCol=1):
    ''' get options from a given options structure '''
    import wx
    if len(options) == 0:
        raise exceptions.ValueError("Empty field names...")    # some behaviors
    # values, not the final result
    # first set them with command line options etc
    values = []
    processedArgs = []
    # process all options
    goto = 0
    for opt in range(0, len(options)):
        if opt < goto:
            values.append(None)
            continue
        p = options[opt]
        # validate p
        for k in p.keys():
            if not k in allowed_keys:
                raise exceptions.ValueError("Unrecognized option entry " + k )
        if p.has_key('separator'):
            val = p['separator']
        else:
            val = _getParamShortArg(p, processedArgs)
            if val == None:
                val = _getParamLongArg(p, processedArgs)
            if val == None:
                val = _getParamConfigFile(p, processedArgs)
            if val == None:
                if p.has_key('default'):
                    val = p['default']
            # ignore jump options
        values.append(val)
    # look if any argument was not processed
    if checkUnprocessedArgs:
        for i in range(1, len(sys.argv)):
            if (not sys.argv[i] in allowed_commandline_options) and (not i in processedArgs):
                raise exceptions.ValueError("Unprocessed command line argument: " + sys.argv[i])
    # parameter 0 prevents wxPython from open a separate window for stdout and stderr
    app = wx.App(0)
    dlg = wx.Dialog(parent=None, id=-1, title=title)
    entryWidgets = [None]*len(options)
    labelWidgets = [None]*len(options)
    # cancel
    def onCancel(event):
        " when ESC is pressed cancel "
        # scripple values by changing its length
        values.append(None)
        dlg.EndModal(wx.ID_CANCEL)
    # done
    def onOK(event):
        ''' get result and convert values '''
        for g in range(len(entryWidgets)):
            if entryWidgets[g] == None:
                continue
            try:
                # get text from different type of entries
                try:    # an entry box or check box
                    val = _getParamValue( options[g], entryWidgets[g].GetValue())
                except:
                    try:    # a list box?
                        val = _getParamValue( options[g],
                            options[g]['chooseOneOf'][int(entryWidgets[g].GetSelection())])
                    except: # a checklist box?
                        items = []
                        for s in range(len(options[g]['chooseFrom'])):
                            if entryWidgets[g].IsChecked(s):
                                items.append( options[g]['chooseFrom'][s])
                        val = _getParamValue( options[g], items)
            except exceptions.Exception, e:
                # incorrect value
                # set to red
                # clear other red colors
                for lab in labelWidgets:
                    if lab != None:
                        lab.SetForegroundColour('black')
                # set this one to red
                labelWidgets[g].SetForegroundColour('red')
                entryWidgets[g].SetFocus()
                return
            else:
                # convert to values
                values[g] = val
        # get all results and return
        dlg.EndModal(wx.ID_OK)
    # help
    def onHelp(event):
        # open another dialog
        dlg1 = wx.Dialog(parent=dlg, id=-1, title='Help for ' + title)
        #global root1
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add( wx.TextCtrl( parent=dlg1, id=-1, size=[600,400],
            style=wx.TE_MULTILINE | wx.TE_READONLY,
            value=usage(options, details)), 0, wx.ALL, 20)
        okButton = wx.Button( parent=dlg1, id=wx.ID_OK, label='OK')
        box.Add( okButton, 0, wx.ALIGN_CENTER | wx.ALL, 20)
        app.Bind( wx.EVT_BUTTON, lambda event:dlg1.EndModal(wx.ID_OK),
            okButton)
        dlg1.SetSizerAndFit(box)
        dlg1.Layout()
        dlg1.ShowModal()
        dlg1.Destroy()
    # format a description to tooltip
    def formatDesc(text):
        # linux can auto wrap, windows can not but sometime wrap
        # at unexpected places... It is safer to wrap at original
        # place.
        return '\n'.join( [x.strip() for x in text.splitlines()] )
    #
    # the main window
    #
    # get font height for dialog arrangement purpose
    # top message
    box = wx.BoxSizer(wx.VERTICAL)
    # do not use a very long description please
    topLabel = wx.StaticText(parent=dlg, id=-1, label='\n'+description)
    box.Add( topLabel, 0, wx.EXPAND | wx.LEFT | wx.ALIGN_CENTER , 50)
    #topLabel.SetForegroundColour('Black')
    topLabel.SetFont(wx.Font(12, wx.SWISS, wx.NORMAL, wx.NORMAL))
    # add a box for all ...
    paraBox = wx.FlexGridSizer(cols=nCol)
    for b in range(nCol):
        paraBox.AddGrowableCol(b)
    #
    # add several FlexGridSizer
    gridBox = []
    for b in range(nCol):
        gridBox.append( wx.FlexGridSizer(cols=2, vgap=5, hgap=20))
        gridBox[-1].AddGrowableCol(0)
        gridBox[-1].AddGrowableCol(1)
        paraBox.Add( gridBox[-1], 1, wx.EXPAND | wx.ALL, 10)
    box.Add(paraBox, 1, wx.EXPAND | wx.ALL, 10)
    # count numbers of valid parameters..
    # chooseFrom count as three
    colParam = 0
    for opt in options:
        if opt.has_key('label') or opt.has_key('separator'):
            colParam += 1
        if opt.has_key('chooseFrom'):
            colParam += len( opt['chooseFrom']) -2
    if colParam / nCol * nCol == colParam:
        colParam /= nCol
    else:
        colParam = colParam/nCol + 1
    colCount = 0
    colIndex = 0
    #print colParam
    # all entries
    for g in range(len(options)):
        opt = options[g]
        if not (opt.has_key('label') or opt.has_key('separator')) :
            continue
        colIndex = colCount / colParam
        if opt.has_key('separator'):
            labelWidgets[g] = wx.StaticText(parent=dlg, id=-1, label=opt['separator'])
            #labelWidgets[g].SetForegroundColour('Blue')
            labelWidgets[g].SetFont(wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
            gridBox[colIndex].Add(labelWidgets[g], 0, wx.ALIGN_LEFT )
            entryWidgets[g] = None
            gridBox[colIndex].Add(wx.StaticText(parent=dlg, id=-1, label=''), 1, wx.ALIGN_LEFT )
            colCount += 1
            continue
        else: # label
            labelWidgets[g] = wx.StaticText(parent=dlg, id=-1, label=opt['label'])
            gridBox[colIndex].Add(labelWidgets[g], 0, wx.ALIGN_LEFT )
        # --------- entryWidget ----------------------------------------------
        # use different entry method for different types
        if opt.has_key('chooseOneOf'):    # single choice
            entryWidgets[g] = wx.Choice(parent=dlg, id=g, choices = opt['chooseOneOf'])
            if opt.has_key('description'):
                entryWidgets[g].SetToolTipString(formatDesc(opt['description']))
            gridBox[colIndex].Add(entryWidgets[g], 1, wx.EXPAND )
            if values[g] != None:
                try:
                    entryWidgets[g].SetSelection(opt['chooseOneOf'].index(values[g]))
                except:
                    raise ValueError('Value: %s is not one of %s.' % (str(values[g]), str(opt['chooseOneOf'])))
            colCount += 1
        elif opt.has_key('chooseFrom'):    # multiple choice
            w,h = labelWidgets[g].GetTextExtent('a')
            entryWidgets[g] = wx.CheckListBox(parent=dlg, id=g, size=(0, (h+4)*len(opt['chooseFrom'])),
                choices = opt['chooseFrom'])
            if opt.has_key('description'):
                entryWidgets[g].SetToolTipString(formatDesc(opt['description']))
            if values[g] != None:
                if type(values[g]) in [types.ListType, types.TupleType]:
                    for val in values[g]:
                        entryWidgets[g].Check( opt['chooseFrom'].index(val))
                else:
                    entryWidgets[g].Check( opt['chooseFrom'].index(values[g]))
            gridBox[colIndex].Add(entryWidgets[g], 1, wx.EXPAND)
            colCount += len(opt['chooseFrom']) -1
        elif (opt.has_key('arg') and opt['arg'][-1] != ':') or \
             (opt.has_key('longarg') and opt['longarg'][-1] != '='):  # true or false
            w,h = labelWidgets[g].GetTextExtent('a')
            entryWidgets[g] = wx.CheckBox(parent=dlg, id=g, label = 'Yes / No')
            if opt.has_key('description'):
                entryWidgets[g].SetToolTipString(formatDesc(opt['description']))
            if values[g] != None:
                entryWidgets[g].SetValue(values[g])
            gridBox[colIndex].Add(entryWidgets[g], 1, wx.EXPAND)
            colCount += 1
        else: # a edit box
            # put default value into the entryWidget
            txt = ''
            if values[g] != None:
             if type(values[g]) in [types.ListType, types.TupleType] and len(values[g])>1:
                 txt =    ', '.join(map(str, values[g]))
             else:
                 txt = str(values[g])
            entryWidgets[g] = wx.TextCtrl(parent=dlg, id=g, value=txt)
            if opt.has_key('description'):
                entryWidgets[g].SetToolTipString(formatDesc(opt['description']))
            gridBox[colIndex].Add(entryWidgets[g], 1, wx.EXPAND )
            colCount += 1
    # help button
    buttonBox = wx.GridSizer(cols=3)
    helpButton = wx.Button(dlg, wx.ID_HELP, 'Help')
    buttonBox.Add( helpButton, 0, wx.ALIGN_CENTER )
    dlg.Bind( wx.EVT_BUTTON, onHelp, helpButton)
    cancelButton = wx.Button(dlg, wx.ID_CANCEL, 'Cancel')
    buttonBox.Add( cancelButton, 0, wx.ALIGN_CENTER )
    dlg.Bind( wx.EVT_BUTTON, onCancel, cancelButton)
    okButton = wx.Button(dlg, wx.ID_OK, 'OK')
    buttonBox.Add( okButton, 0, wx.ALIGN_CENTER)
    dlg.Bind( wx.EVT_BUTTON, onOK, okButton)
    #
    box.Add( buttonBox, 0, wx.ALL | wx.EXPAND, 20)
    # ------------------- time for action! -----------------
    dlg.SetSizerAndFit(box)
    dlg.Layout()
    # first un-none
    for g in range(len(options)):
        if entryWidgets[g] != None:
            entryWidgets[g].SetFocus()
            break
    dlg.ShowModal()
    dlg.Destroy()
    # we do not actually need a main loop
    #app.MainLoop()    # run it!
    # -------- after the run has completed ----------------------------------
    if len(values) == len(options):
        # remove values inserted by separators
        ret = []
        for p in range(len(options)):
            if not options[p].has_key('separator'):
                ret.append( values[p] )
        return ret
    else:
        return []


# get parameter
def getParam(options=[], doc="", details="", noDialog=False, checkUnprocessedArgs=True, verbose=False, nCol=1):
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

        checkUnprocessedArgs: check args, avoid misspelling of arg name

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
                
    if noDialog or '--noDialog' in sys.argv[1:] or '-h' in sys.argv[1:] or '--help' in sys.argv[1:] \
        or True not in map(lambda x:x.has_key('label'), options):
        return _termGetParam(options, doc, verbose)
    else:
        title = os.path.split(sys.argv[0])[-1]
        if useTkinter:
            return _tkGetParam(options, title, doc, details,
                checkUnprocessedArgs, nCol)
        elif useWxPython:
            return _wxGetParam(options, title, doc, details,
                checkUnprocessedArgs, nCol)
        else:
            return _termGetParam(options, doc, verbose)


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
    message += '        --mpi :\n                      Use MPI module\n'
    message += '        --chromMap : \n                Chromosome map for MPI module\n'
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
        if p.has_key('default') and p['default'] != None:
            message +=    '                Default to '
            if type(p['default']) in [types.ListType, types.TupleType]:
                message += ', '.join(map(str, p['default']))
            else:
                message += str(p['default'])
            message += '\n'
    return message


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
                    print >> f, "# shortarg:\t-%s = %s" % (options[p]['arg'][0:-1], str( param[p] ))
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
                    arg = options[p]['longarg'][0:-1]
                else:
                    arg = options[p]['longarg']
                # write out option value, try to make it python readable
                if type(param[p]) == type(''):
                    if "'" in param[p]:
                        print >> f, '%s = "%s"' % (options[p]['longarg'][0:-1], param[p])
                    else:
                        print >> f, "%s = '%s'" % (options[p]['longarg'][0:-1], param[p])
                else:
                    print >> f, "%s = %s" % (options[p]['longarg'][0:-1], str(param[p]))
    f.write("\n\n#The same options can be given by command line (subject to minor changes)\n#    --noDialog ")
    for p in range(0, len(options)):
        if options[p].has_key('label') and options[p].has_key('longarg'):
            if options[p]['longarg'][-1] == '=':
                if str(param[p]).find(",") >= 0:    # has single quote
                    f.write( " --" + options[p]['longarg'][0:-1]
                        + '="' + str( param[p] ) + '"')
                else:
                    f.write( " --" + options[p]['longarg'][0:-1]
                        + "='" + str( param[p] ) + "'")
            else:
                f.write( " --" + options[p]['longarg'] )
    f.write("\n")
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
env_mpi = os.getenv('SIMUMPI')
env_chrom_map = os.getenv('SIMUCHROMMAP')
env_longAllele = os.getenv('SIMUALLELETYPE')
env_debug = os.getenv('SIMUDEBUG')

[par_optimized] = _termGetParam([{'longarg':'optimized', \
    'default':''}], False, False, True)
[par_mpi] = _termGetParam([{'longarg':'mpi', \
    'default':''}], False, False, True)
[par_chrom_map] = _termGetParam([{'longarg':'chromMap', \
    'default':[]}], False, False, True)
[par_quiet] = _termGetParam([{'arg':'q','longarg':'quiet', \
    'default':False}], False, False, True)
[par_useTkinter] = _termGetParam([{'longarg':'useTkinter', \
    'default':False }], False, False, True) 

# remove these parameters from sys.argv
for arg in ['--optimized', '--mpi', '--chromMap', '--quiet', '-q', '--useTkinter']:
    try:
        sys.argv.remove(arg)
    except:
        pass

if par_optimized != '':
    _optimized = par_optimized
elif env_optimized != None:
    _optimized = True
else:     # default to false
    _optimized = False

if par_mpi != '':
    _mpi = par_mpi
elif env_mpi != None:
    _mpi = True
else:     # default to false
    _mpi = False

if par_chrom_map != []:
    _chrom_map = par_chrom_map
elif env_chrom_map != None:
    _chrom_map = eval(env_chrom_map)
else:     # default to []
    _chrom_map = []

if env_longAllele in ['standard', 'short', 'long', 'binary']:
    _longAllele = env_longAllele
else:
    _longAllele = 'standard'

simuOptions = {'Optimized':_optimized, 'MPI':_mpi, 'ChromMap':_chrom_map, 
    'AlleleType':_longAllele, 'Debug':[], 'Quiet':par_quiet}

if env_debug != None:
    simuOptions['Debug'].extend( env_debug.split(',') )


def setOptions(optimized=None, mpi=None, chromMap=[], alleleType=None, quiet=None, debug=[]):
    '''set options before simuPOP is loaded to control which simuPOP module to load,
    and how the module should be loaded.

    optimized: whether or not load optimized version of a module. If not set,
        environmental variable SIMUOPTIMIZED, and commandline option --optimized
        will be used if available. If nothing is defined, standard version will
        be used.

    mpi: currently unused

    chromMap: currently unused

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
    if mpi in [True, False]:
        simuOptions['MPI'] = mpi
    if chromMap != []:
        simuOptions['ChromMap'] = chromMap
    if alleleType in ['standard', 'long', 'binary', 'short']:
        simuOptions['AlleleType'] = alleleType
    if quiet in [True, False]:
        simuOptions['Quiet'] = quiet
    if len(debug) > 0:
        aimuOptions['Debug'].extend(debug)

# short = standard
if simuOptions['AlleleType'] == 'standard':
    simuOptions['AlleleType'] = 'short'
if simuOptions['Optimized'] not in [True, False]:
    simuOptions['Optimized'] = False
if simuOptions['MPI'] not in [True, False]:
    simuOptions['MPI'] = False

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
    except:
        useWxPython = False
    else:
        useWxPython = True


if par_useTkinter or not useWxPython:
    # Tkinter should almost always exists, but ...
    try:
        imp.find_module('Tkinter')
    except:
        print "Tkinter can not be loaded. Please check your Python installation."
        useTkinter = False
    else:
        useTkinter = True
        # this is not possible now because of the use of find_module
        #if TkVersion < 8.0 :
        #    useTkinter = False
