#!/usr/bin/env python

#
# $File: simuUtil.py $
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


"""
simuPOP utilities.

This module provides some commonly used operators
and format conversion utilities.

"""

import exceptions, operator, types, os, sys, re

from simuPOP import *
from simuOpt import simuOptions

def ViewVars(var, gui=None):
    '''
    list a variable in tree format, either in text format or in a
        wxPython window.

    var
        any variable to be viewed. Can be a dw object returned
        by dvars() function

    level
        level of display.

    name
        only view certain variable

    subPop
        whether or not display info in subPop

    '''
    if gui is None:
        gui = simuOptions['GUI']
    #
    if gui == False or gui == 'Tkinter':
        import pprint
        pprint.pprint(var)
        return

    try:
        import wx, wx.py.filling as fill
    except ImportError:
        pprint.pprint(var)
        return

    app = wx.PySimpleApp()
    wx.InitAllImageHandlers()
    if var==None:
        fillFrame = fill.FillingFrame()
    else:
        if type(var) == type( dw({}) ):
            fillFrame = fill.FillingFrame(rootObject=var.__dict__,
                rootLabel='var')
        else:
            fillFrame = fill.FillingFrame(rootObject=var,
                rootLabel='var')
    fillFrame.Show(True)
    app.MainLoop()


# migration rate matrix generators
def MigrIslandRates(r, n):
    '''migration rate matrix

::

         x m/(n-1) m/(n-1) ....
         m/(n-1) x ............
         .....
         .... m/(n-1) m/(n-1) x

    where x = 1-m
    '''
    # n==1?
    if n == 1:
        return [[1]]
    #
    m = []
    for i in range(0,n):
        m.append([r/(n-1.)]*n)
        m[-1][i] = 1-r
    return m


def MigrHierarchicalIslandRates(r1, r2, n):
    '''
    Return the migration rate matrix for a hierarchical island model
    where there are different migration rate within and across groups
    of islands.

    r1
        Within group migration rates. It can be a number or a list of numbers
        for each group of the islands.

    r2
        Across group migration rates which is the probability that someone will
        migrate to a subpopulation outside of his group. A list of r2 could be
        specified for each group of the islands.

    n
        Number of islands in each group. E.g. n=[5, 4] specifies two groups of
        islands with 5 and 4 islands each.

    For individuals in an island, the probability that it remains in the same
    island is 1-r1-r2 (r1, r2 might vary by island groups), that it migrates
    to another island in the same group is r1 and migrates to another island
    outside of the group is r2. Migrate rate to a specific island depends on
    the size of group.
    '''
    if type(n) not in [type(()), type([])]:
        raise ValueError('A list of size of island groups is expected for parameter n')
    nIslands = sum(n)
    if type(r1) in [type(0), type(1.)]:
        r1 = [r1] * len(n)
    elif len(r1) != len(n):
        raise ValueError('If multiple r1 is given, it should be given to all island groups.')
    #
    if type(r2) in [type(0), type(1.)]:
        r2 = [r2] * len(n)
    elif len(r2) != len(n):
        raise ValueError('If multiple r2 is given, it should be given to all island groups.')
    #
    m = []
    for groupIdx, groupSize in enumerate(n):
        nOther = nIslands - groupSize
        groupStart = sum(n[:groupIdx])
        groupEnd = groupStart + groupSize
        for island in range(groupStart, groupEnd):
            m.append([])
            for i in range(groupStart):
                m[-1].append(r2[groupIdx] * 1.0 / nOther)
            for i in range(groupStart, groupEnd):
                if i == island:
                    m[-1].append(1 - r1[groupIdx] - r2[groupIdx])
                else:
                    m[-1].append(r1[groupIdx] * 1.0 / groupSize)
            for i in range(groupEnd, nIslands):
                m[-1].append(r2[groupIdx] * 1.0 / nOther)
    return m


def MigrSteppingStoneRates(r, n, circular=False):
    '''migration rate matrix, circular stepping stone model (X=1-m)

::

           X   m/2               m/2
           m/2 X   m/2           0
           0   m/2 x   m/2 ......0
           ...
           m/2 0 ....       m/2  X

or non-circular

::

           X   m/2               m/2
           m/2 X   m/2           0
           0   m/2 X   m/2 ......0
           ...
           ...              m   X
    '''
    if n < 2:
        raise exceptions.ValueError("Can not define step stone model for n < 2")
    elif n == 2:
        return [[1-r,r],[r,1-r]]
    # the normal case (n>2)
    m = []
    for i in range(0, n):
        m.append([0]*n)
        m[i][i] = 1-r
        m[i][(i+1)%n] = r/2.
        m[i][(i+n-1)%n] = r/2.
    if not circular:
        m[0][1] = r
        m[0][-1] = 0
        m[n-1][0] = 0
        m[n-1][n-2] = r
    return m


def SaveCSV(pop, filename='', fields=[], loci=[], header=True,
    shift=1,
        sexCode={Male: '1', Female: '2'}, affectionCode={True: '1', False: '2'},
        **kwargs):
    '''Save a simuPOP population ``pop`` in csv format.

    pop
        A simuPOP population object. If a string is given, it will be loaded.

    filename
        Output filename.

    fileds
        Information fields to be outputted.

    loci
        If a list of loci is given, only genotype at these loci will be
        written.

    header
        Whether or not a header should be written. These headers will include
        information fields, sex (if ``sexCode`` is not ``None``), affection
        status (if ``affectionCode`` is not ``None``) and loci names. If
        genotype at a locus needs more than one column, ``_1``, ``_2`` etc will
        be appended to locus names.

    genotype
        list of loci to output, default to all.

    combine
        how to combine the markers. Default to None.
        A function can be specified, that takes the form::

             def func(markers):
                 return markers[0]+markers[1]

    shift
        since alleles in simuPOP is 0-based, shift=1 is usually needed to
        output alleles starting from allele 1. This parameter is ignored if
        combine is used.

    '''
    if output != '':
        file = output
    elif outputExpr != '':
        file = eval(outputExpr, globals(), pop.vars() )
    else:
        raise exceptions.ValueError, "Please specify output or outputExpr"
    if loci == []:
        loci = range(0, pop.totNumLoci())
    try:
        out = open( file, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file +" to write."
    # keep the content of pieces in strings first
    content = [''] * pop.numChrom()
    # write out header
    print >> out, 'id, ', ', '.join(fields), ', ',
    if combine is None:
        print >> out, ', '.join(['%s_1, %s_2' % (pop.locusName(loc), pop.locusName(loc)) for loc in loci])
    else:
        print >> out, ', '.join(['%s' % pop.locusName(loc) for loc in loci])
    # write out
    id = 1
    pldy = pop.ploidy()
    for ind in pop.individuals():
        print >> out, id,
        for f in fields:
            if f == 'sex':
                print >> out, ',', sexCode[ind.sex()],
            elif f == 'affection':
                print >> out, ',', affectionCode[ind.affected()],
            else:
                print >> out, ',', ind.info(f),
        for marker in loci:
            if combine is None:
                for p in range(pldy):
                    out.write(", %d" % (ind.allele(marker, p) + shift))
            else:
                out.write(", %d" % combine([ind.allele(marker, p) for p in range(pldy)]))
        print >> out
        id += 1
    out.close()


class _baseProgressBar:
    def __init__(self, message, totalCount):
        '''
        message
            Title of the progress bar

        totalCount
            Total expected steps.

        done
            Message displayed when the job is finished.
        '''
        self.message = message
        self.totalCount = totalCount
        self.percent = 0
        self.completed = False

    def update(self, count):
        '''
        Update the progreebar.
        '''
        count = min(count, self.totalCount)
        self.progress = int(round(100*count/self.totalCount))
        if self.progress <= self.percent:
            return False
        else:
            return True

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if self.completed:
            return False
        else:
            self.completed = True
            return True

class _textProgressBar(_baseProgressBar):
    def __init__(self, message, totalCount, progressChar='.', block=2, done=' Done.\n'):
        '''
        message
            Title of the progress bar

        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        _baseProgressBar.__init__(self, message, totalCount)
        self.percent = 0
        self.progressChar = progressChar
        self.block = block
        self.doneMsg = done
        sys.stdout.write(message)
        sys.stdout.flush()

    def update(self, count):
        '''
        Update the text progreebar.
        '''
        if not _baseProgressBar.update(self, count):
            return
        for p in range(self.percent + 1, self.progress + 1):
            if p == 100:
                self.done()
            elif p % 10 == 0:
                sys.stdout.write(str(p/10))
            elif p % self.block == 0:
                sys.stdout.write(self.progressChar)
        sys.stdout.flush()
        self.percent = self.progress
        if self.percent == 100:
            self.done()

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if not _baseProgressBar.done(self):
            return
        sys.stdout.write(self.doneMsg)
        sys.stdout.flush()


class _tkProgressBar(_baseProgressBar):
    def __init__(self, message, totalCount):
        '''
        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        _baseProgressBar.__init__(self, message, totalCount)
        import Tkinter as tk
        self.width = 300
        self.height = 30
        self.max = 100
        self.fillColor = 'blue'
        self.labelColor = 'black'
        self.label = 'Progress'
        #
        self.app = tk.Tk()
        self.app.title(self.label)
        self.frame = tk.Frame(self.app, bd=0)
        self.canvas = tk.Canvas(self.frame, bd=0, width=self.width+40,
            height = self.height + 70, highlightthickness=0)
        self.label = self.canvas.create_text(20, 20, 
            text='', anchor="w", fill=self.labelColor, font=('Verdana', 10))
        self.scale = self.canvas.create_rectangle(
            20, 50, self.width + 20, 50 + self.height, fill=self.fillColor)
        self.rect = self.canvas.create_rectangle(
            20, 50, self.width + 20, 50 + self.height)
        self.canvas.pack(side='top', fill='x', expand='yes', padx=0)
        self.update(0)
        self.frame.pack(padx=0, pady=0)

    def update(self, count):
        '''
        Update the progreebar.
        '''
        if not _baseProgressBar.update(self, count):
            return
        #
        self.canvas.coords(self.scale, 20, 50,
            20 + self.progress * 1.0 / self.max * self.width, 50 + self.height)
        # Now update the colors
        self.canvas.itemconfig(self.scale, fill=self.fillColor)
        self.canvas.itemconfig(self.label, fill=self.labelColor)
        # And update the label
        if self.progress > 0:
            self.canvas.itemconfig(self.label, text=self.message + "\n%d%% completed." % self.progress)
        else:
            self.canvas.itemconfig(self.label, text=self.message)
        self.canvas.update_idletasks()
        self.app.update()
        #
        self.percent = self.progress
        if self.percent == 100:
            self.done()

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if not _baseProgressBar.done(self):
            return
        self.app.destroy()
        del self.app


class _wxProgressBar(_baseProgressBar):
    def __init__(self, message, totalCount):
        '''
        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        _baseProgressBar.__init__(self, message, totalCount)
        import wx
        self.app = wx.PySimpleApp(0)
        self.dialog = wx.ProgressDialog(
            'Progress', self.message + '\n', self.totalCount,
            style = \
                # wx.PD_CAN_ABORT | \
                # wx.PD_CAN_SKIP | \
                wx.PD_ELAPSED_TIME | \
                # wx.PD_ESTIMATED_TIME | \
                wx.PD_AUTO_HIDE | \
                wx.PD_REMAINING_TIME
            )
        self.dialog.Update(0)

    def update(self, count):
        '''
        Update the progreebar.
        '''
        if not _baseProgressBar.update(self, count):
            return
        self.dialog.Update(count, self.message + "\n%d%% completed." % self.progress)
        self.percent = self.progress
        if self.percent == 100:
            self.done()

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        if not _baseProgressBar.done(self):
            return
        self.dialog.Destroy()
        del self.app



class simuProgress:
    '''
    This class defines a very simple text based progress bar. It will display a
    character (default to "``.``") for each change of progress (default to 2%),
    and a number (1, 2, ..., 9) for each 10% of progress, and print a message
    (default to "``Done.\\n``") when the job is finished.

    This class is used as follows::

        progress = simuProgress("Start simulation", 500)
        for i in range(500):
            progress.update(i+1)
        # if you would like to make sure the done message is displayed.
        progress.done()
    '''
    def __init__(self, message, totalCount, progressChar='.', block=2, done=' Done.\n', gui=None):
        '''
        message
            Title of the progress bar.

        totalCount
            Total expected steps.

        progressChar
            Character to be displayed for each progress.

        block
            display progress at which interval (in terms of percentage)?

        done
            Message displayed when the job is finished.
        '''
        if gui is None:
            self.gui = simuOptions['GUI']
        else:
            self.gui = gui
        if self.gui in ['wxPython', True]:
            try:
                import wx
                self.gui = 'wxPython'
            except ImportError:
                self.gui = 'Tkinter'
        if self.gui == 'Tkinter':
            try:
                import Tkinter
            except ImportError:
                self.gui = False
        if self.gui == 'wxPython':
            self.progressBar = _wxProgressBar(message, totalCount)
        elif self.gui == 'Tkinter':
            self.progressBar = _tkProgressBar(message, totalCount)
        else:
            self.progressBar = _textProgressBar(message, totalCount, progressChar, block, done)

    def update(self, count):
        '''
        Update the progreebar.
        '''
        self.progressBar.update(count)

    def done(self):
        '''
        Finish progressbar, print 'done' message.
        '''
        self.progressBar.done()


class pySubset(pyOperator):
    '''
    This operator rearranges and removes individuals according to their values at
    an information field. Individuals with positive values at this information
    field are moved to the subpopulation specified by the integer value of this
    value. Individuals with negative values are removed. There is no function
    form of this operator because this operator is essentially a wrapper around
    function ``population::setSubPopByIndInfo(field)``.
    '''
    def __init__(self, field, *args, **kwargs):
        '''
        Create a pySubset operator that rearranges and removes individuals
        according to their values at an information field *field*.
        '''
        self.field = field
        pyOperator.__init__(self, func=self.apply,
            param=self.field, *args, **kwargs)

    def apply(self, pop, field):
        pop.setSubPopByIndInfo(field)
        return True

    def __repr_():
        return "<simuPOP::pySubset>"

    def clone():
        return pySubset(self.field)

  
if __name__ == "__main__":
    pass
