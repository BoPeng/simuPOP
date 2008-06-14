#!/usr/bin/env python

############################################################################
#    Copyright (C) 2004 by Bo Peng
#    bpeng@rice.edu
#
#    $LastChangedDate$
#    $Rev$
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
############################################################################



"""
This module helps the use of  rpy  package with simuPOP. It defines an
operator  varPlotter  that can be used to plot population expressions
when  rpy  is installed.
"""

from exceptions import *
from math import ceil, sqrt

try:
    import rpy_options
    rpy_options.set_options(VERBOSE = False)
    from rpy import *
except:
    print "Rpy can not be loaded. Please verify your rpy installation. "
    print "If you are using a rpy version that reuires Python Numeric package. "
    print "Please install Numeric first. "
    raise ImportError("Failed to import rpy package. (You will need rpy>=0.99)")


import os, sys

# if under windows, fix a bug with rpy which uses blocking i/o so
# R figure will not be refreshed in time
if os.name == 'nt':
    r.options(windowsBuffered=False)

from simuPOP import *
from simuUtil import dataAggregator

def rmatrix(mat):
    ''' Convert a Python 2d list to r matrix format
        that can be passed to functions like image
        directly.
    '''
    #return with_mode(NO_CONVERSION, r.do_call)('rbind', mat)
    return r.do_call('rbind', mat)


class _VarPlotter_Base:
    def __init__(self, nplot, update, title, xlab, ylab, axes, lty, col,
        mfrow, plotType, saveAs, leaveOpen, dev='', width=0, height=0):
        """
        initialization function of base properties (layout etc)
        of all plotters
        """
        # save parameters
        self.nplot = nplot
        self.axes = axes
        if lty==[]:
            self.lty = range(1, nplot+1)
        else:
            self.lty = lty
        if col==[]:
            self.col = [0]*nplot
        else:
            self.col = col
        self.mfrow = [int(x) for x in mfrow]
        self.update = update
        self.title = [title]*self.nplot
        self.xlab = xlab
        self.ylab = ylab
        self.plotType = plotType
        self.saveAs = saveAs
        self.leaveOpen = leaveOpen
        self.width = width
        self.height = height
        # layout related variables
        self.first = True
        # R plot device number
        self.dev = dev
        self.device = 1
        # update related variable
        self.skip = -1

    def CanUpdate(self, inc=True):
        if inc:
            self.skip += 1
        if self.update < 1:    # at every gen
            return True
        # if update > 1
        if self.update != 1 and self.skip != self.update:
            return False
        else:
            self.skip = 0
            return True

    def setDev(self):
        # if dev is provided
        if self.first:
            if self.dev != '':
                if self.width != 0 and self.height != 0:
                    r.postscript(file=self.dev, width=self.width, height=self.height)
                else:
                    r.postscript(file=self.dev)
            else:
                # open a new window
                r('get(getOption("device"))()')
                # get device number
                self.device = r.dev_cur()
                assert self.device > 1, 'Can not open new device'
            self.first = False
        else:
            if self.dev == '':
                # use the device of this plotter
                r.dev_set( self.device )

    def __del__(self) :
        # If postscript file, or not leaveOpen
        if not self.leaveOpen or self.dev != '':
            r.dev_off()

    def layout(self):
        # calculate layout
        if self.mfrow != [1,1] and self.mfrow[0]*self.mfrow[1] < self.nplot:
            raise ValueError("mfrow is not enough to hold " + str(self.nplot) + " figures")

        if self.mfrow == [1,1] and self.nplot > 1:
            self.mfrow[0] = int(ceil(sqrt( self.nplot)))
            self.mfrow[1] = int(ceil(self.nplot/float(self.mfrow[0])))
            if self.mfrow[0] > self.mfrow[1]:
                     self.mfrow.reverse()

        # use par in simple case
        if self.plotType    == "plot":
            r.par( mfrow=self.mfrow)
            return

        # else, use layout.
        # 2 3 4 1
        # 5 6 7 1
        # 8 9 0 1
        #
        # where 1 is color bar, 2-9 are plots.
        lt = [0] * ( self.mfrow[0] * (self.mfrow[1]+1))
        for i in range(0, self.mfrow[0]):
            for j in range(0, self.mfrow[1]+1):
                if i*self.mfrow[1]+2+j <= self.nplot + 1:
                    lt[i* (self.mfrow[1]+1)+j] = i*self.mfrow[1]+2+j
                else:
                    lt[i* (self.mfrow[1]+1)+j] = 0
            lt[(i+1)*(self.mfrow[1]+1)-1] = 1

        w = 7 * with_mode(BASIC_CONVERSION, r.par)("csi") * 2.54
        r.assign('lt', lt)
        r('''layout(matrix(lt, byrow=TRUE, nc = %d+1),
            widths = c(rep("1", %d), lcm(%f)))''' % (self.mfrow[1], self.mfrow[1], w))

    def colorBar(self, level, lim):
        " draw a color bar to the right of the plot"
        mar = [0]*4
        mar[0] = [5, 1, 4, 2]

        self.color = r.rainbow(level, start=.7, end=.1)
        r.par(mar=mar)
        levels = with_mode(BASIC_CONVERSION, r.seq)(lim[0], lim[1], length=level)
        r.plot_new()
        r.plot_window(xlim = [0, 1], ylim = lim,
            xaxs = "i", yaxs = "i")
        r.rect(0, levels[:(len(levels)-1)], 1, levels[1:],
            col = self.color)
        r.axis(4)
        r.box()
        mar[3] = 1
        r.par(mar = mar)

    def save(self, gen):
        " save plots "
        if self.saveAs == "":
            return
        if gen < 0:
            r.dev_print(file=self.saveAs + '.eps')
        else:
            r.dev_print(file=self.saveAs + str(gen) + '.eps')

class _VarPlotter_NoHis(_VarPlotter_Base):
    """
    a class to plot variable, potentially from different replicate
    DO NOT PLOT HISTORY.
    The idea is:

     when intiailize set some local variables
     each time during apply, call self.plot(gen, rep, var)
     var will be considered variable fromr eplicate rep.

     If the dimension of var is bigger than one, use par to create more than
     one features. (has to be specified by dim.

     parameters:
         win:         generation window. If 0, plot all generation
         ylim:        limit of y axis. Will be the range of the first replicate if not given.
         update:    update figure after 'update' generations.
         title:     plot title. Will spread to other subplots if a vector is given.
         xlab:    title of x axis. Will spread to other subplots if a vector is given.
         ylab:    title of y axis. Will speread to other subplots if a vector is given.


         byRep:     Whether or not use different subplot for different replicate. Default
                            to false. If set to be true and there are more than one replicate,
                            you should set numRep and maybe mfrow.
         numRep:    number of replicate. If byRep is true and there are more than
                            one replicates, set numRep here and varPlotter will try to set
                            reasonaly subplot dimention, if dim is not given.
         mfrow:     if there are more than one subplots, mdrow=[x,y] specifies rows and
                            columns of the subplots

    """

    def __init__(self, numRep=1, win=0, ylim=[0,0], update=1, plotType="plot", lty=[], col=[],
                             title="variable", xlab="generation", ylab="", axes=True, mfrow=[1,1],
                             byRep=0, saveAs="", leaveOpen=True, level = 20, dev='', width=0, height=0):
        """

        """
        self.byRep = byRep

        if byRep:
            _VarPlotter_Base.__init__(self, numRep, update, title, xlab,
                ylab, axes, lty, col, mfrow, plotType, saveAs, leaveOpen, dev, width, height)
        else:
            _VarPlotter_Base.__init__(self, 1, update, title, xlab,
                ylab, axes, lty, col, mfrow, plotType, saveAs, leaveOpen, dev, width, height)

        self.numRep = numRep

        if self.byRep and self.numRep > 1:
            for d in range(0, self.numRep):
                self.title[d] = self.title[d] + ", rep " + str(d)

        self.win = win
        self.ylim = ylim
        self.xlim = [0,0]
        if ylim == [0,0]:
            self.dynamicLim = 1
        else:
            self.dynamicLim = 0
        self.level = level    # for image

    def plot(self, pop, expr):
        gen = pop.gen()
        rep = pop.rep()
        data = pop.evaluate(expr)
        #
        if type(data) == type(0) or type(data) == type(0.):
            data = [data]

        if rep >= self.numRep:
            raise ValueError("Replicate out of range. (forget to set numRep?)")

        # plot data
        if rep==-1 or rep == 0:
            # inc ...
            if not self.CanUpdate(True):
                return True
        else:    # regular replicate
            if not self.CanUpdate(False):
                return True
        self.setDev()
        if rep == 0:
            self.layout()

        if self.dynamicLim and self.plotType != 'image':
            self.xlim[0] = 0
            self.xlim[1] = len(data)-1

            self.ylim[0] = data[0]
            self.ylim[1] = data[0]
            for i in range(1, len(data)):
                if data[i] < self.ylim[0]:
                    self.ylim[0] = data[i]
                if data[i] > self.ylim[1]:
                    self.ylim[1] = data[i]
        elif self.plotType == 'image':     # if plot type is image?
            self.zlim=[data[0][0], data[0][0]]
            for i in range(1, len(data)):
                for j in range(1, len(data[i])):
                    if data[i][j] < self.zlim[0]:
                        self.zlim[0] = data[i][j]
                    if data[i][j] > self.zlim[1]:
                        self.zlim[1] = data[i][j]

        if self.plotType == 'plot':
            if self.byRep or rep <= 0 :
                r.plot( range(0, len(data)), data,
                    xlab=self.xlab, ylab=self.ylab,axes=self.axes,
                    main=self.title[rep], type='l', xlim=self.xlim,
                    ylim=self.ylim, lty=1)
            else:
                r.lines( range(0, len(data)), data, type='l', lty=rep+1)
        elif self.plotType == 'image':
            if rep == 0:
                self.setDev()
                self.colorBar(self.level, self.zlim)
            d=[]
            for i in range(0,len(data)):
                d.extend(data[i])
            r.image( z= with_mode(NO_CONVERSION, r.t)(
                with_mode(NO_CONVERSION, r.matrix)(d, ncol=len(data))),
                xlab=self.xlab, axes=self.axes,
                ylab=self.ylab, main=self.title[rep], col=self.color)
        else:
            raise ValueError("Only 'plot' or 'image' is supported as plotType")

        if self.saveAs != "":
            self.save(self.saveAs, gen)
        return True

# with history
class _VarPlotter_His(_VarPlotter_Base):
    """
    a class to plot variable, potentially from different replicate
    The idea is:

     when intiailize set some local variables
     each time during apply, call self.plot(gen, rep, var)
     var will be considered variable fromr eplicate rep.

     If the dimension of var is bigger than one, use par to create more than
     one features. (has to be specified by dim.

     parameters:
         win:         generation window. If 0, plot all generation
         ylim:        limit of y axis. Will be the range of the first replicate if not given.
         update:    update figure after 'update' generations.
         title:     plot title. Will spread to other subplots if a vector is given.
         xlab:    title of x axis. Will spread to other subplots if a vector is given.
         ylab:    title of y axis. Will speread to other subplots if a vector is given.


         byRep:     Whether or not use different subplot for different replicate. Default
                            to false. If set to be true and there are more than one replicate,
                            you should set numRep and maybe mfrow.
         numRep:    number of replicate. If byRep is true and there are more than
                            one replicates, set numRep here and varPlotter will try to set
                            reasonaly subplot dimention, if dim is not given.
         byVal:     Whether or not use different subplot for each dimension of val.
                            Default to false. If set to be true and the dimension of data
                            is more than one, you should set varDim and maybe mfrgw.
         varDim:    dimension of variable. Need to be set if byVal is true.
         mfrow:     if there are more than one subplots, mdrow=[x,y] specifies rows and
                            columns of the subplots

         separate: if true, use different panel (within a figure) to display each
                            item of an array.

         plotType: usual plot or image
         level:        level of color, default to 20

    """

    def __init__(self, varDim=1, numRep=1, win=0, ylim=[0,0], update=1, lty=[], col=[],
                             title="variable", xlab="generation", ylab="", axes=True, mfrow=[1,1],
                             byRep=0, byVal=0, separate=False, plotType="plot", level=20,
                             saveAs="", leaveOpen=True, dev='', width=0, height=0):
        """
        initialize _withHis plotter
        """
        if (not byRep) and numRep > 1 and (not byVal) and varDim > 1:
            raise ValueError("Can not plot multi-dimension data for multiple replicate in one figure.\n" + \
                "Use byRep or byVal variable.")

        if byVal and byRep:
            raise ValueError("Can not separate variables both by replicate and by var dimension.")

        self.byRep = byRep
        self.byVal = byVal

        self.data = []
        nplot = 1
        if byRep:
            self.numRep = numRep
            nplot = numRep
            for rep in range(0, numRep):
                self.data.append(dataAggregator(maxRecord=win, recordSize=varDim))
        elif byVal:
            nplot = varDim
            for d in range(0, varDim):
                self.data.append(dataAggregator(maxRecord=win,recordSize=numRep))
        else:
            # otherwise, all data in one figure
            if numRep > 1:
                self.data.append(dataAggregator(win, numRep))
                self.byVal = 1
            else:
                self.data.append(dataAggregator(win, varDim))
                self.byRep = 1

        if lty==[]:
            lty = [1]*numRep
        if col==[]:
            col = ['black']*numRep
        # initialize base plotter
        _VarPlotter_Base.__init__(self, nplot, update, title,
            xlab, ylab, axes, lty, col, mfrow, plotType, saveAs, leaveOpen, dev, width, height)
        #
        self.varDim = varDim
        self.numRep = numRep
        self.separate = separate
        #
        if len(self.data) > 1 :
            for d in range(0, len(self.data)):
                if byRep:
                    self.title[d] = self.title[d] + ", rep " + str(d)
                else:
                    self.title[d] = self.title[d] + ", " + str(d)
        #
        self.win = win
        self.ylim = ylim
        self.xlim = [0,0]
        self.level = level
        if ylim == [0,0]:
            self.dynamicLim = 1
        else:
            self.dynamicLim = 0

    def clear(self):
        self.data=[]

    def plot(self, pop, expr):
        gen = pop.gen()
        rep = pop.rep()
        data = pop.evaluate(expr)
        _data = data
        # now start!
        if type(_data) == type(0) or type(_data) == type(0.):
            _data = [_data]

        if type(data) == type({}):
            raise ValueError("Variable is a dictionary. A vector is expected.")

        if len(_data) != self.varDim:
            print "Variable dimension is wrong (forget to set varDim?) varDim: ", \
                self.varDim, " data: ", str(data)
            if len(_data) > self.varDim:
                _data = data[:self.varDim]
            else:
                for i in range(0, self.varDim - len(_data)):
                    _data.append(0)
            print "Change it to ", _data

        if rep >= self.numRep:
            raise ValueError("Replicate out of range. (forget to set numRep?)")

        if self.byRep: # push all data in to an replicate-specific aggregator
            if rep >= self.numRep:
                raise ValueError("Replicate number out of range. (Use numRep parameter.")
            self.data[rep].push( gen, _data)
        else: # separate data
            for d in range(0, self.varDim):
                self.data[d].push( gen, _data[d], rep)

        # plot data, only plot at the last replicate
        if rep==-1 or rep == self.numRep-1:
            # last replicate, but can not update
            if not self.CanUpdate() or len( self.data[0].gen) == 1:
                return True
        # only update at the last replicate, so...
        else:
            return True

        self.setDev()
        self.layout()

        self.xlim[0] = self.data[0].gen[0]
        self.xlim[1] = self.data[0].gen[-1]

        if self.win > 0 and self.xlim[1] - self.xlim[0] < self.win:
            self.xlim[1] = self.xlim[0] + self.win

        if self.dynamicLim:
            self.ylim = self.data[0].dataRange()
            for i in range(1, len(self.data)):
                yl = self.data[i].dataRange()
                if yl[0] < self.ylim[0]:
                    self.ylim[0] = yl[0]
                if yl[1] > self.ylim[1]:
                    self.ylim[1] = yl[1]

        if self.byRep:
            if self.plotType == "image":                # use image
                self.colorBar(self.level, self.ylim)
                for rep in range(0, self.numRep):
                    r.image(x=self.data[rep].gen, y=range(0, self.data[rep].recordSize),
                        xlim=self.xlim,
                        z=with_mode(NO_CONVERSION, r.t)(
                            with_mode(NO_CONVERSION, r.matrix)(
                                self.data[rep].flatData(), byrow=True, ncol=len(self.data[rep].gen)
                                )
                            ), xlab=self.xlab,
                        ylab=self.ylab, main=self.title[rep], col=self.color)
            elif self.separate == False or self.data[rep].recordSize == 1:
                for rep in range(0, self.numRep):
                    r.plot( self.data[rep].gen, self.data[rep].data[0],
                        xlab=self.xlab, ylab=self.ylab,axes=self.axes,
                        main=self.title[rep], type='l', xlim=self.xlim,
                        ylim=self.ylim, lty=self.lty[0], col=self.col[0])
                    for i in range( 1, self.data[rep].recordSize ):
                        r.lines( self.data[rep].gen, self.data[rep].data[i], type='l', lty=self.lty[i],
                            col=self.col[i])
            else:    # use panel for each item of the plot
                height = self.ylim[1] - self.ylim[0]
                newlim = [self.ylim[0], self.ylim[0] + height* self.data[0].recordSize]
                for rep in range(0, self.numRep):
                    r.plot( self.data[rep].gen, self.data[rep].data[0],
                        xlab=self.xlab, ylab=self.ylab,axes=self.axes,
                        main=self.title[rep], type='l', xlim=self.xlim,
                        ylim=newlim, lty=self.lty[0], col=self.col[0])
                    r.abline(h=[newlim[0] + height*x for x in range(0, self.data[rep].recordSize)])
                    for i in range( 1, self.data[rep].recordSize ):
                        r.lines( self.data[rep].gen,
                            [height*i+x for x in self.data[rep].data[i]], type='l', lty=self.lty[i],
                            col=self.col[i])

        if self.byVal:
            if self.plotType == "image":
                self.colorBar(self.level,self.ylim)
                for rep in range(0, self.numRep):
                    r.image(x=self.data[rep].gen, y=range(0, self.data[rep].recordSize),
                        xlim=self.xlim, z= r.t(r.matrix(self.data[rep].flatData(),
                        byrow=True, ncol=len(self.data[rep].data[0]))), xlab=self.xlab,
                        axes=self.axes,
                        ylab=self.ylab, main=self.title[rep], col=self.color)
            elif self.separate == False or self.data[rep].recordSize == 1:
                for v in range(0, self.varDim):
                    r.plot( self.data[v].gen, self.data[v].data[0],
                        xlab=self.xlab, ylab=self.ylab,axes=self.axes,
                        main=self.title[v], type='l', xlim=self.xlim,
                        ylim=self.ylim, lty=self.lty[0], col=self.col[0])
                    for i in range( 1, self.data[v].recordSize ):
                        r.lines( self.data[v].gen, self.data[v].data[i], type='l',
                            lty=self.lty[i], col=self.col[i])
            else:
                height = self.ylim[1] - self.ylim[0]
                newlim = [self.ylim[0], self.ylim[0] + height* self.data[0].recordSize]
                for v in range(0, self.varDim):
                    r.plot( self.data[v].gen, self.data[v].data[0],
                        xlab=self.xlab, ylab=self.ylab,axes=self.axes,
                        main=self.title[v], type='l', xlim=self.xlim,
                        ylim=newlim, lty=self.lty[0], col=self.col[0])
                    r.abline(h=[newlim[0] + height*x for x in range(0, self.data[rep].recordSize)])
                    for i in range( 1, self.data[v].recordSize ):
                        r.lines( self.data[v].gen,
                            [height*i + x for x in self.data[v].data[i]], type='l',
                            lty=self.lty[i], col=self.col[i])
        self.save(gen)
        return True

# currently, type can be plot, image
#
class varPlotter(pyOperator):
    '''
    Plotting with history

    plot a number in the form of a variable or expression, use
        >>> varPlotter(var='expr')

    plot a vector in the same window and there is only one replicate in
    the simulator, use

        >>> varPlotter(var='expr', varDim=len)

    where len is the dimension of your variable or expression. Each line
    in the figure represents the history of an item in the array.

    plot a vector in the same window and there are several replicates, use

        >>> varPlotter(var='expr', varDim=len, numRep=nr, byRep=1)

    varPlotter will try to use an appropriate layout for your subplots
    (for example, use 3x4 if  numRep=10  ). You can also specify parameter
    mfrow to change the layout.

    if you would like to plot each item of your array variables in a subplot,
    use

        >>> varPlotter(var='expr', varDim=len, byVal=1)

    or in case of a single replicate

        >>> varPlotter(var='expr', varDim=len, byVal=1, numRep=nr)


    There will be numRep lines in each subplot.

    Use option   history=False  to plot with history. Parameters   byVal  ,  varDim  etc. will be ignored.

    Other options are

    title, xtitle, ytitle: title of your figure(s). title is default to your expression, xtitle is defaulted to generation.

    win: window of generations. I.e., how many generations to keep in a figure. This is useful when you want to keep track of only recent changes.

    update: update figure after update generations. This is used when you do not want to update the figure at every generation.

    saveAs: save figures in files saveAs#gen.eps. For example, if saveAs='demo', you will get files demo1.eps, demo2.eps etc.

    separate: plot data lines in separate panels.

    image: use R image function to plot image, instead of lines.

    level: level of image colors (default to 20).

    leaveOpen: whether or not leave the plot open when plotting is done. Default to True.

    '''
    def __init__(self, expr, history=True, varDim=1, numRep=1, win=0, ylim=[0,0],
        update=1, title="", xlab="generation", ylab="",axes=True, lty=[], col=[],
        mfrow=[1,1], separate=False, byRep=False, byVal=False, plotType="plot",
        level=20, saveAs="", leaveOpen=True, dev='', width=0, height=0,
        *args, **kwargs):
        self.expr = expr
        if history == True:
            self.plotter = _VarPlotter_His(numRep=numRep, varDim=varDim, win=win, ylim=ylim,
                update=update, title=title, xlab=xlab, ylab=ylab, axes=axes, mfrow=mfrow,
                byRep=byRep, byVal=byVal, separate=separate, plotType=plotType,
                level=level, saveAs=saveAs, lty=lty, col=col,
                leaveOpen=leaveOpen, dev=dev, width=width, height=height)
        else:
            self.plotter = _VarPlotter_NoHis(numRep=numRep, win=win, ylim=ylim,
                update=update, title=title, xlab=xlab, ylab=ylab, axes=axes, mfrow=mfrow,
                byRep=byRep, plotType=plotType, level=level, saveAs=saveAs, lty=lty, col=col,
                leaveOpen=leaveOpen, dev=dev, width=width, height=height)
        # when apply is called, self.plotter.plot is called
        # *args, **kwargs may be things like rep=, gro=, at=, begin=
        pyOperator.__init__(self, func=self.plotter.plot,
            param=self.expr, *args, **kwargs)


