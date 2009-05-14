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
    import rpy
except ImportError, e:
    print 'Rpy can not be loaded. Please verify your rpy installation.'
    print 'Note that rpy > 0.99 is needed and rpy2 is not supported'
    raise e

import os

# if under windows, fix a bug with rpy which uses blocking i/o so
# R figure will not be refreshed in time.
if os.name == 'nt':
    rpy.r.options(windowsBuffered=False)

from simuPOP import pyOperator, PostMating

class varPlotter(pyOperator):
    '''
    This class defines a Python operator that uses R to plot the current and
    historical values of a Python expression (``expr``). When this operator is
    applied to populations during evolution, an expression is evaluated at each
    population's local namespace. All results are saved unless parameter
    ``win`` is used to specify a window of generations to keep. The return
    value of the expression can be a number or a sequence, but should have the
    same type and length across all replicates and generations.
    
    A figure will be draw at the end of the last replicate (except for the
    first generation where no line could be drawn) unless the current
    generation is less than ``update`` generations away from the last
    generation at which a figure has been drawn. Historical values for all
    replicates will be displayed as lines with the x xais being generation
    number. Multiple lines for each replicate will be drawn if the results of
    the expression are sequences. These lines could be plotted in the same
    figure, or seperated to subplots by replicates (``byRep``), by each
    dimention of the results (``byDim``), or both. These figure could be saved
    to files in various formats if parameter ``saveAs`` is specified. File
    format is determined by file extension.

    Besides parameters mentioned above, arbitrary keyword parameters could be
    specified and be passed to the underlying R functions ``plot`` and
    ``line``. These parameters could be used to specify line type (``lty``),
    color (``col``), title (``main``), limit of x and y axes (``xlim`` and
    ``ylim``) and many other options (see R manual for details). As a special
    case, multiple values can be passed to each replicate or dimension if the
    name of a parameter ends with ``_rep`` or ``_dim``. For example,
    ``lty_rep=range(1, 5)`` will pass parameters ``lty=1``, ..., ``lty=4`` to
    four replicates.
    '''
    def __init__(self, expr, win=0, update=1, byRep=False, byDim=False,
        saveAs="", leaveOpen=True, legend=[], stage=PostMating, begin=0,
        end=-1, step=1, at=[], rep=[], **kwargs):
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
            
        byDim
            Separate items from sequence results of ``expr`` to different
            subplots. If both ``byRep`` and ``byDim`` are ``True``, the
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

        legend
            labels of the lines. This operator will look for keyword parameters
            such as ``col``, ``lty``, ``lwd``, and ``pch`` and call the
            ``legend`` function to draw a legend. If figure has multiple lines
            for both replicates and dimensions, legends should be given to each
            dimension, and then each replicate.

        **kwargs
            All additional keyword parameters will be passed directly to the
            plot function ``r.plot`` and ``r.line``. Such parameters includes
            but not limited to ``xlab``, ``ylab``, ``main``, ``xlim``,
            ``ylim``, ``col``, ``lty``. Multiple values could be passed to
            different replicates or dimensions of results if suffix ``_rep``
            or ``_dim`` is appended to parameter name. For example,
            ``col_rep=['red', 'blue']`` uses two colors for values from
            different replicates. ``_rep_dim`` and ``_dim_rep`` could also be
            used. If the list has insufficient number of items, existing items
            will be reused. 
        '''
        # parameters
        self.expr = expr
        self.win = win
        self.update = update
        self.byRep = byRep
        self.byDim = byDim
        self.saveAs = saveAs
        self.leaveOpen = leaveOpen
        self.legend = legend
        # internal flags
        self.nRep = 0
        self.reps = []   # allows specification of selected replicates
        self.nDim = 0
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
            rpy.r.dev_off()

    def getDev(self):
        # open a new window
        try:
            # 46754 is the revision number for R 2.8.0
            if int(rpy.r.R_Version()['svn rev']) < 46754:
                # For R < 2.8.0, getOption('device') returns a string (such as 'X11')
                rpy.r(rpy.r.getOption('device') + '()')
            else:
                # For R >= 2.8.0, getOption('device') returns a function
                rpy.r('getOption("device")()')
        except:
            raise RuntimeError("Failed to get R version to start a graphical device");
        #
        # get device number
        self.device = rpy.r.dev_cur()
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
        # check if data type and length are consistent, set self.nDim
        if type(data) in [type(()), type([])]:
            if self.nDim == 0:
                self.nDim = len(data)
            if self.nDim != len(data):
                raise RuntimeError('Data dimension is inconsistent.')
        elif self.nDim > 1:
            raise RuntimeError('Data dimension is inconsistent.')
        # append data
        self.data[rep].append(data)
        # check number of saved generations (self.win)
        if self.win > 0 and len(self.gen) > 1 and rep + 1 == len(self.data) \
            and self.gen[0] + self.win < gen:
            self.gen.pop(0)
            for d in self.data:
                if len(d) > 0:
                    d.pop(0)
        # set self.nRep
        if self.nRep == 0 and len(self.gen) > 1:
            self.reps = [x for x in range(len(self.data)) if len(self.data[x]) > 0]
            self.nRep = len(self.reps)

    def getData(self, rep, dim = 0):
        "Get the dim'th element of the data of replicate rep"
        if type(self.data[rep][0]) in [type(()), type([])]:
            return [x[dim] for x in self.data[rep]]
        else:
            return self.data[rep]

    def getArgs(self, rep, dim, kwargs, **default):
        "Get the single format parameters for keyword parameters."
        ret = {}
        for key,value in kwargs.iteritems():
            if key.endswith('_rep_dim'):
                par = key.rstrip('_rep_dim')
                idx = self.reps.index(rep) * self.nDim + dim
            elif key.endswith('_dim_rep'):
                par = key.rstrip('_dim_rep')
                idx = self.reps.index(rep) + self.nRep * dim
            elif key.endswith('_rep'):
                par = key.rstrip('_rep')
                idx = self.reps.index(rep)
            elif key.endswith('_dim'):
                par = key.rstrip('_dim')
                idx = dim
            else:
                ret[key] = value
                continue
            if (type(value) not in [type(()), type([])]):
                ret[par] = value
            else:
                ret[par] = value[idx % len(value)]
        for key,value in default.iteritems():
            if not ret.has_key(key):
                ret[key] = value
        return ret

    def getLegendArgs(self, legendType, kwargs, **default):
        ret = {}
        for var in ['lty', 'col', 'lwd', 'pch', 'bty']:
            ret[var] = []
            if legendType == '_rep':
                for i in range(self.nRep):
                    arg = self.getArgs(i, 0, kwargs, **default)
                    if arg.has_key(var):
                        ret[var].append(arg[var])
            elif legendType == '_dim':
                for i in range(self.nDim):
                    arg = self.getArgs(0, i, kwargs, **default)
                    if arg.has_key(var):
                        ret[var].append(arg[var])
            elif legendType == '_rep_dim':
                for i in range(self.nRep):
                    for j in range(self.nDim):
                        arg = self.getArgs(i, j, kwargs, **default)
                        if arg.has_key(var):
                            ret[var].append(arg[var])
            if len(ret[var]) == 0:
                ret.pop(var)
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
        if self.nRep == 1:
            self.byRep = False
        if self.nDim == 1:
            self.byDim = False
        # needs subplots?
        nPlots = 1
        if self.byDim:
            nPlots *= self.nDim
        if self.byRep:
            nPlots *= self.nRep
        # call r.par to create subplots
        if nPlots > 1:
            nrow = int(ceil(sqrt(nPlots)))
            ncol = int(ceil(nPlots/float(nrow)))
            if nrow > ncol:
                nrow, ncol = ncol, nrow
            # call r.par to allocate subplots
            rpy.r.par(mfrow=[nrow, ncol])
        # now plot.
        if self.byRep:
            # handle each replicate separately
            for i in self.reps:
                if self.byDim:
                    # separate plot for each dim
                    for j in range(self.nDim):
                        rpy.r.plot(self.gen, self.getData(i, j), **self.getArgs(i, j, kwargs, type='l'))
                else:
                    # all var in one subplot
                    rpy.r.plot(self.gen, self.getData(i, 0), **self.getArgs(i, 0, kwargs, type='l'))
                    for j in range(1, self.nDim):
                        rpy.r.lines(self.gen, self.getData(i, j), **self.getArgs(i, j, kwargs))
                    if len(self.legend) > 0:
                        rpy.r.legend('topright', legend=self.legend,
                            **self.getLegendArgs('_dim', kwargs, bty='n', lty=1))
        else:
            # all replicate in one figure
            if self.byDim:
                for i in range(self.nDim):
                    rpy.r.plot(self.gen, self.getData(self.reps[0], i), **self.getArgs(self.reps[0], i, kwargs, type='l'))
                    for j in self.reps[1:]:
                        rpy.r.lines(self.gen, self.getData(j, i), **self.getArgs(j, i, kwargs))
                    if len(self.legend) > 0:
                        rpy.r.legend('topright', legend=self.legend,
                            **self.getLegendArgs('_rep', kwargs, bty='n', lty=1))
            else:
                rpy.r.plot(self.gen, self.getData(0, 0), **self.getArgs(0, 0, kwargs, type='l'))
                for i in self.reps:
                    for j in range(self.nDim):
                        rpy.r.lines(self.gen, self.getData(i, j), **self.getArgs(i, j, kwargs))
                if len(self.legend) > 0:
                    rpy.r.legend('topright', legend=self.legend,
                        **self.getLegendArgs('_rep_dim', kwargs, bty='n', lty=1))
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
                device = rpy.r.pdf
            elif ext.lower() == '.png':
                device = rpy.r.png
            elif ext.lower() == '.bmp':
                device = rpy.r.bmp
            elif ext.lower() == '.jpg' or ext.lower() == '.jpeg':
                device = rpy.r.jpeg
            elif ext.lower() == '.tif' or ext.lower() == '.tiff':
                device = rpy.r.tiff
            elif ext.lower() == '.eps':
                device = rpy.r.postscript
            else:
                device = rpy.r.postscript
        except Exception, e:
            print e
            print 'Can not determine which device to use to save file %s. A postscript driver is used.' % name
            device = rpy.r.postscript
        #
        if gen < 0:
            if ext == '':
                name += '.eps'
        else:
            if ext == '':
                name += str(gen) + '.eps'
            else:
                name = name.replace(ext, str(gen) + ext)
        rpy.r.dev_print(file=name, device = device)

