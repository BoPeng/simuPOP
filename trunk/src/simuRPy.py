#!/usr/bin/env python

#
# $File: simuRPy.py $
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
This module defines an operator ``varPlotter`` that makes use of the Python
rpy module (http://rpy.sourceforge.net) to plot expressions in R
(http://www.r-project.org), a popular statistical analysis language. Note that
rpy2, the successor of rpy, is currently not supported.
'''

from exceptions import RuntimeError
from math import ceil, sqrt

try:
    import rpy_options
    rpy_options.set_options(VERBOSE = False)
    from rpy import *
except ImportError, e:
    print 'Rpy can not be loaded. Please verify your rpy installation.'
    print 'Note that rpy > 0.99 is needed and rpy2 is not supported'
    raise e

import os

# if under windows, fix a bug with rpy which uses blocking i/o so
# R figure will not be refreshed in time.
if os.name == 'nt':
    r.options(windowsBuffered=False)

from simuPOP import pyOperator

class varPlotter(pyOperator):
    '''
    This class defines a Python operator that uses R to plot the values of a
    Python expression. When this operator is applied to a population at
    different generations, this expression is evaluated with its return values
    expression is evaluated in its local namespace. Its value 
    How this expression is plotted depends on the dimension of the
    return value (if a sequence is returned), number of replicates, whether or
    not historical values (collected over several generations) are plotted,
    and plot type (lines or images).

    The default behavior of this operator is to plot the history of an
    expression. For example, when operator

        varPlotter(var='expr')

    is used in simulator::evolve, the value of ``expr`` will be recorded each
    time when this operator is applied. A line will be draw in a figure with
    x-axis being the generation number. Parameters ``ylim`` can be used to
    specify the range of y-axis.

    If the return value of expression ``expr`` is a sequence (tuple or list),
    parameter ``varDim`` has to be used to indicate the dimension of this
    expression. For example,

        varPlotter(var='expr', varDim=3)

    will plot three lines, corresponding to the histories of each item in the
    array.

    If the expression returns a number and there are several replicates,
    parameter ``numRep``` should be used. In this case, each line will
    correspond to a replicate.

    If the expression returns a vector and there are several replicates, several
    subplots will be used. Parameter ``byRep`` or ``byVal`` should be used
    to tell ``varPlotter`` whether the subplots should be divided by replicate
    or by variable. For example,

        varPlotter(var='expr', varDim=8, numRep=5, byRep=1)

    will use an appropriate layout for your subplots, which is, in this case,
    2x3 for 5 replicates. Each subplot will have 8 lines. If byVal is ``True``,
    there will be 3x3 subplots for 8 items in an array, and each subplot will
    have 5 lines. Note that ``byRep`` or ``byVal`` can also be used when there
    is only one replicate or if the dimension of the expression is one.

    When ``history=False``, histories of each variable will be discarded so
    the figure will always plot the current value of the expression.

    '''
    def __init__(self, expr, win=0, update=1, byRep=False, byVal=False,
        saveAs="", leaveOpen=True, stage=PostMating, begin=0, end=-1, step=1,
        at=[], rep=[], repArgs=['lty', 'col', 'xlab', 'ylab', 'main'],
        valArgs=[], **kwargs):
        '''
        expr
            expression that will be evaluated at each replicate's local
            namespace when the operator is applied. Its value can be a number
            or a list (or tuple) but the type and length of the return value
            should be consistent for all replicates and at all generations.

        win
            Window of generations. If given, only values from generation -win
            to -1 will be plotted.

        update
            Update the figure after specified generations. For example, you can
            evalulate an expression and save its values at every 10 generations
            (parameter ``step=10``) but only draw a figure after every 50
            generations (parameter ``update=50``.

        byRep
            Separate values at different replicates to different subplots.
            
        byVal
            Separate items from sequence results of ``expr`` to different
            subplots. If both ``byRep`` and ``byVal`` are ``True``, the
            subplots will be arranged by variable and then replicates.

        saveAs
            Save figures in files saveAs.ext. If ext (such as pdf, jpeg) is
            given, a corresponding device will be used. Otherwise, a postscript
            driver will be used. Generation number will be inserted before file
            extension so 'figure.eps' will produce files such as 'figure10.eps',
            'figure50.eps'.

        leaveOpen
            Whether or not leave the plot open when plotting is done. Default
            to ``True``. This allows further manipulate of the figures using R
            functions.

        kwargs
            All additional keyword parameters will be passed directly to the
            plot function ``r.plot``. Such parameters can be xlab, ylab,
            main, xlim, ylim, col, lty.... If the name of such an arg appear in
            ``repArgs`` and/or ``valArgs`` and a list of parameters are given,
            this operator will try to use the specified values for each
            replicate and/or dimension of the returned values of ``expr``.
            For example, ``lty=range(1,5)`` will use different line types for
            different lines.

        repArgs
            The name of keyword parameters whose list form will be handled
            specially for each replicate.

        valArgs
            The name of keyword parameters whose list form will be handled
            specially for each dimension of the returned values of ``expr``.
        '''
        # parameters
        self.expr = expr
        self.win = win
        self.update = update
        self.byRep = byRep
        self.byVal = byVal
        self.saveAs = saveAs
        self.leaveOpen = leaveOpen
        self.repArgs = repArgs
        self.valArgs = valArgs
        # internal flags
        self.numRep = 0
        self.numVal = 0
        self.lastPlot = 0
        # data
        self.gen = []
        self.data = []
        # when apply is called, self.plot is called, additional keyword
        # parameters are passed by kwargs.
        pyOperator.__init__(self, func=self.plot, param=kwargs,
            begin=begin, end=end, step=step, at=at,
            rep=rep, subPops=[], infoFields=[])

    def __del__(self) :
        # Close the device if needed.
        if not self.leaveOpen and hasattr(self, 'device'):
            r.dev_off()

    def getDev(self):
        # open a new window
        try:
            # 46754 is the revision number for R 2.8.0
            if int(r.R_Version()['svn rev']) < 46754:
                # For R < 2.8.0, getOption('device') returns a string (such as 'X11')
                r(r.getOption('device') + '()')
            else:
                # For R >= 2.8.0, getOption('device') returns a function
                r('getOption("device")()')
        except:
            raise RuntimeError("Failed to get R version to start a graphical device");
        #
        # get device number
        self.device = r.dev_cur()
        if self.device == 0:
            raise RuntimeError('Can not open new device')

    def pushData(self, gen, rep, data):
        'Push history data to self.data for later retrieval'
        # allocate a list for each replicate.
        while rep + 1 > len(self.data):
            self.data.append([])
        # append gen.
        if len(self.gen) == 0 or self.gen[-1] != gen:
            self.gen.append(gen)
        # check if data type and length are consistent, set self.numVal
        if type(data) in [type(()), type([])]:
            if self.numVal == 0:
                self.numVal = len(data)
            if self.numVal != len(data):
                raise RuntimeError('Data dimension is inconsistent.')
        elif self.numVal > 1:
            raise RuntimeError('Data dimension is inconsistent.')
        # append data
        self.data[rep].append(data)
        # check number of saved generations (self.win)
        if self.win > 0 and len(self.gen) > 1 and rep + 1 == len(self.data) \
            and self.gen[0] + self.win < gen:
            self.gen.pop(0)
            for d in self.data:
                d.pop(0)
        # set self.numRep
        if self.numRep == 0 and len(self.gen) > 1:
            self.numRep = len(self.data)

    def getData(self, rep, dim = 0):
        "Get the dim'th element of the data of replicate rep"
        if type(self.data[rep][0]) in [type(()), type([])]:
            return [x[dim] for x in self.data[rep]]
        else:
            return self.data[rep]

    def getArgs(self, rep, dim, kwargs):
        "Get the single format parameters for keyword parameters"
        ret = {}
        for key,value in kwargs.iteritems():
            if (type(value) not in [type(()), type([])]) or \
                (key not in self.repArgs and key not in self.repArgs):
                ret[key] = value
                continue
            idx = 0
            if key in self.repArgs:
                idx = rep
            if key in self.valArgs:
                idx = idx * self.numVal + dim
            if len(value) >= idx + 1:
                ret[key] = value[idx]
            else:
                ret[key] = value[0]
        return ret

    def plot(self, pop, kwargs):
        "Evaluate expression in pop and save result. Plot all data if needed"
        gen = pop.dvars().gen
        rep = pop.dvars().rep
        # push data 
        self.pushData(gen, rep, pop.evaluate(self.expr))
        # Draw a plot only when
        # 1. There are at least two obervations.
        # 2. rep is the last recorded replicate.
        # 3. we are self.update away from last plot.
        if len(self.gen) <= 1 or rep + 1 != len(self.data) or \
            (self.update >= 1 and gen < self.lastPlot + self.update):
            # do not plot
            return True
        else:
            self.lastPlot = gen
        # create a new graphical device if needed
        if not hasattr(self, 'device'):
            self.getDev()
        # figure out the dimension of data
        if self.numRep == 1:
            self.byRep = False
        if self.numVal == 1:
            self.byVal = False
        # needs subplots?
        nPlots = 1
        if self.byVal:
            nPlots *= self.numVal
        if self.byRep:
            nPlots *= self.numRep
        # call r.par to create subplots
        if nPlots > 1:
            nrow = int(ceil(sqrt(nPlots)))
            ncol = int(ceil(nPlots/float(nrow)))
            if nrow > ncol:
                nrow, ncol = ncol, nrow
            # call r.par to allocate subplots
            r.par(mfrow=[nrow, ncol])
        # now plot.
        if self.byRep:
            # handle each replicate separately
            for i in range(self.numRep):
                if self.byVal:
                    # separate plot for each dim
                    for j in range(self.numVal):
                        r.plot(self.gen, self.getData(i, j), **self.getArgs(i, j, kwargs))
                else:
                    # all var in one subplot
                    r.plot(self.gen, self.getData(i, 0), **self.getArgs(i, 0, kwargs))
                    for j in range(1, self.numVal):
                        r.lines(self.gen, self.getData(i, j), **self.getArgs(i, j, kwargs))
        else:
            # all replicate in one figure
            if self.byVal:
                for i in range(self.numVal):
                    r.plot(self.gen, self.getData(0, i), **self.getArgs(0, i, kwargs))
                    for j in range(1, self.numRep):
                        r.lines(self.gen, self.getData(j, i), **self.getArgs(j, i, kwargs))
            else:
                r.plot(self.gen, self.getData(0, 0), **self.getArgs(0, 0, kwargs))
                for i in range(self.numRep):
                    for j in range(self.numVal):
                        r.lines(self.gen, self.getData(i, j), **self.getArgs(i, j, kwargs))
        self.save(gen)
        return True

    def save(self, gen):
        "Save plots using a device specified by file extension."
        if self.saveAs == "":
            return
        name = self.saveAs
        ext = os.path.splitext(name)[-1]
        try:
            # I need to use this more lengthy form because some
            # functions are not available in, for example, R 2.6.2
            if ext.lower() == '.pdf':
                device = r.pdf
            elif ext.lower() == '.png':
                device = r.png
            elif ext.lower() == '.bmp':
                device = r.bmp
            elif ext.lower() == '.jpg' or ext.lower() == '.jpeg':
                device = r.jpeg
            elif ext.lower() == '.tif' or ext.lower() == '.tiff':
                device = r.tiff
            elif ext.lower() == '.eps':
                device = r.postscript
            else:
                device = r.postscript
        except Exception, e:
            print e
            print 'Can not determine which device to use to save file %s. A postscript driver is used.' % name
            device = r.postscript
        #
        if gen < 0:
            if ext == '':
                name += '.eps'
        else:
            if ext == '':
                name += str(gen) + '.eps'
            else:
                name = name.replace(ext, str(gen) + ext)
        r.dev_print(file=name, device = device)

