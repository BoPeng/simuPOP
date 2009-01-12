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
This module helps the use of ``rpy`` package with simuPOP. It defines an
operator ``varPlotter`` that can be used to plot population expressions
when ``rpy`` is installed.
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

class _dataAggregator:
    """
    collect variables so that plotters can plot them all at once

    Usage

            a = _dataAggregator( maxRecord=0, recordSize=0)

                maxRecord
                    if more data is pushed, the old ones are discarded

                recordSize
                    size of record

            a.push(gen, data, idx=-1)

                gen
                    generation number

                data
                    one record (will set recordSize if the first time), or

                idx
                    if idx!=-1, set data at idx.

            a.clear()
            a.range()    # return min, max of all data
            a.data[i]    # column i of the data
            a.gen        #
            a.ready()    # if all column has the same length, so data is ready


    Internal data storage::
            self.gen   
            self.data   column1
            column2

    each record is pushed at the end of

    """
    def __init__(self, maxRecord=0, recordSize=0):
        """maxRecord: maxRecorddow size. I.e., maximum generations of data to keep"""
        self.gen = []
        self.data = []
        self.maxRecord = maxRecord
        self.recordSize = recordSize

    def __repr__(self):
        s = str(self.gen) + "\n"
        for i in range(0, len(self.data)):
            s += str(self.data[i]) + "\n"
        return s

    def clear(self):
        self.gen = []
        self.data = []

    def ready(self):
        return self.recordSize>0 and len(gen)>0 and len( data[0] ) == len( data[-1] )

    def flatData(self):
        res = []
        for d in self.data:
            res.extend( d )
        return res

    def dataRange(self):
        if len(self.gen) == 0:
            return [0,0]

        y0 = min( [ min(x) for x in self.data] )
        y1 = max( [ max(x) for x in self.data] )
        return [y0,y1]

    def push(self, _gen, _data, _idx=-1 ):
        # first add data to allData
        if len(self.gen) == 0:     # the first time
            self.gen = [ _gen ]
            if _idx == -1:        # given a full array of data
                if self.recordSize == 0:
                    self.recordSize = len(_data)
                elif self.recordSize != len(_data):
                    raise exceptions.ValueError("Data length does not equal specfied record size")
                for i in range(self.recordSize):
                    self.data.append( [_data[i]] )
                return
            elif _idx == 0:     # the only allowed case
                if type(_data) in [type(()), type([])]:
                    raise exceptions.ValueError("If idx is specified, _data should not be a list.")
                self.data = [ [_data] ]
                return
            else:                                                # data out of range
                raise exceptions.ValueError("Appending data with wrong idx")
        elif len(self.gen) == 1:             # still the first generation
            if self.gen[-1] == _gen:        # still working on this generation
                if _idx == -1:    # give a full array?
                    raise exceptions.ValueError("Can not reassign data from this generation")
                elif self.recordSize != 0 and    _idx >= self.recordSize:
                    raise exceptions.ValueError("Data exceeding specified record size")
                elif _idx == len(self.data):    # append
                    if type(_data) in [type(()), type([])]:
                        raise exceptions.ValueError("If idx is specified, _data should not be a list.")
                    self.data.append( [_data] )
                elif _idx < len(self.data):    # change exsiting one?
                    raise exceptions.ValueError("You can not change exisiting data")
                else:                                                # data out of range
                    raise exceptions.ValueError("Appending data with wrong idx")
            else:                                                    # go to the next one!
                if self.recordSize == 0:         # not specified
                    self.recordSize = len(self.data)
                elif self.recordSize != len(self.data):
                    raise exceptions.ValueError("The first row is imcomplete")
                self.gen.append( _gen )
                if _idx == -1:        # given a full array of data
                    if self.recordSize != len(_data):
                        raise exceptions.ValueError("Data length does not equal specfied record size")
                    for i in range(self.recordSize):
                        self.data[i].append( _data[i] )
                    return
                elif _idx == 0:     # the only allowed case
                    if type(_data) in [type(()), type([])]:
                        raise exceptions.ValueError("If idx is specified, _data should not be a list.")
                    self.data[0].append(_data)
                    return
                else:                                                # data out of range
                    raise exceptions.ValueError("Appending data with wrong idx")
        else:     # already more than one record
            # trim data if necessary
            if self.maxRecord > 0 :
                if _gen - self.gen[0] >= self.maxRecord:
                    self.gen = self.gen[1:]
                    for i in range(0, self.recordSize):
                        self.data[i] = self.data[i][1:]
            if self.gen[-1] == _gen:     # still this generation
                if _idx == -1:    # give a full array?
                    raise exceptions.ValueError("Can not reassign data from this generation")
                elif _idx >= self.recordSize:
                    raise exceptions.ValueError("Data exceeding specified record size")
                elif _idx < len(self.data):    # change exsiting one?
                    if type(_data) in [type(()), type([])]:
                        raise exceptions.ValueError("If idx is specified, _data should not be a list.")
                    self.data[_idx].append( _data )
                else:                                                # data out of range
                    raise exceptions.ValueError("Appending data with wrong idx")
            else:                                                    # go to the next one!
                self.gen.append( _gen )
                if _idx == -1:        # given a full array of data
                    if self.recordSize != len(_data):
                        raise exceptions.ValueError("Data length does not equal specfied record size")
                    for i in range(self.recordSize):
                        self.data[i].append( _data[i] )
                    return
                elif _idx == 0:     # the only allowed case
                    if type(_data) in [type(()), type([])]:
                        raise exceptions.ValueError("If idx is specified, _data should not be a list.")
                    self.data[0].append(_data)
                    return
                else:                                                # data out of range
                    raise exceptions.ValueError("Appending data with wrong idx")



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
                try:
                    # 46754 is the revision number for R 2.8.0
                    if int(r.R_Version()['svn rev']) < 46754:
                        # For R < 2.8.0, getOption('device') returns a string (such as 'X11')
                        r(r.getOption('device') + '()')
                    else:
                        # For R >= 2.8.0, getOption('device') returns a function
                        r('getOption("device")()')
                except:
                    raise SystemError("Failed to get R version to start a graphical device");
                #
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
        mar = [5, 1, 4, 2]

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
                self.data.append(_dataAggregator(maxRecord=win, recordSize=varDim))
        elif byVal:
            nplot = varDim
            for d in range(0, varDim):
                self.data.append(_dataAggregator(maxRecord=win,recordSize=numRep))
        else:
            # otherwise, all data in one figure
            if numRep > 1:
                self.data.append(_dataAggregator(win, numRep))
                self.byVal = 1
            else:
                self.data.append(_dataAggregator(win, varDim))
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


# currently, type can be plot, image
# without history image
class _VarPlotter_NoHis_Image(_VarPlotter_NoHis):
    """
    a class to plot image, potentially from different replicate
    DO NOT PLOT HISTORY.

    """

    def plot(self, pop, expr):
        gen = pop.dvars().gen
        rep = pop.dvars().rep
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

        self.zlim=[data[0][0], data[0][0]]
        for i in range(1, len(data)):
            for j in range(1, len(data[i])):
                if data[i][j] < self.zlim[0]:
                    self.zlim[0] = data[i][j]
                if data[i][j] > self.zlim[1]:
                    self.zlim[1] = data[i][j]
       
       
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

        if self.saveAs != "":
            self.save(self.saveAs, gen)
        return True


# without history plot
class _VarPlotter_NoHis_Plot(_VarPlotter_NoHis):
    """
    a class to plot variable, potentially from different replicate
    DO NOT PLOT HISTORY.

    """

    def plot(self, pop, expr):
        gen = pop.dvars().gen
        rep = pop.dvars().rep
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

        if self.dynamicLim:
            self.xlim[0] = 0
            self.xlim[1] = len(data)-1

            self.ylim[0] = data[0]
            self.ylim[1] = data[0]
            for i in range(1, len(data)):
                if data[i] < self.ylim[0]:
                    self.ylim[0] = data[i]
                if data[i] > self.ylim[1]:
                    self.ylim[1] = data[i]
                    
        if self.byRep or rep <= 0 :
            r.plot( range(0, len(data)), data,
                xlab=self.xlab, ylab=self.ylab,axes=self.axes,
                main=self.title[rep], type='l', xlim=self.xlim,
                ylim=self.ylim, lty=1)
        else:
            r.lines( range(0, len(data)), data, type='l', lty=rep+1)

        if self.saveAs != "":
            self.save(self.saveAs, gen)
        return True

     
# with history Image 
class _VarPlotter_His_Image(_VarPlotter_His):
    """
    a class to plot image, potentially from different replicate
    
    """

    def plot(self, pop, expr):
        gen = pop.dvars().gen
        rep = pop.dvars().rep
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
          
        if self.byVal: 
			self.colorBar(self.level,self.ylim)
			for rep in range(0, self.numRep):
				r.image(x=self.data[rep].gen, y=range(0, self.data[rep].recordSize),
					xlim=self.xlim, z= r.t(r.matrix(self.data[rep].flatData(),
					byrow=True, ncol=len(self.data[rep].data[0]))), xlab=self.xlab,
					axes=self.axes,
					ylab=self.ylab, main=self.title[rep], col=self.color)
          
        self.save(gen)
        return True


# with history plot 
class _VarPlotter_His_Plot(_VarPlotter_His):
    """
    a class to plot variable, potentially from different replicate

    """

    def plot(self, pop, expr):
        gen = pop.dvars().gen
        rep = pop.dvars().rep
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
            if self.separate == False or self.data[rep].recordSize == 1:
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
            if self.separate == False or self.data[rep].recordSize == 1:
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


#
class varPlotter(pyOperator):
    '''
    This class defines a Python operator that uses R to plot a simuPOP express.
    During the evolution, this express is evaluated in each replicate's local
    namespace. How this expression is plotted depends on the dimension of the
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
    subplots will be used. Parameter ``byRep`` or ``byVar`` should be used
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
    def __init__(self, expr, history=True, varDim=1, numRep=1, win=0, ylim=[0,0],
        update=1, title="", xlab="generation", ylab="",axes=True, lty=[], col=[],
        mfrow=[1,1], separate=False, byRep=False, byVal=False, plotType="plot",
        level=20, saveAs="", leaveOpen=True, dev='', width=0, height=0,
        *args, **kwargs):
        '''
        expr
            expression that will be evaluate at each replicate's local namespace
            when the operator is applied.

        history
            whether or not record and plot the history of an expression. Default
            to True.

        varDim
            If the return value of ``expr`` is a sequence, ``varDim`` should be
            set to the length of this sequence. Default to 1.

        numRep
            Number of replicates of the simulator. Default to 1.

        win
            Window of generations. I.e., how many generations to keep in a figure.
            This is useful when you want to keep track of only recent changes of
            an expression. The default value is 0, which will keep all histories.

        ylim
            The range of y-axis.

        update
            Update figure after update generations. This is used when you do not
            want to update the figure every time when this operator is applied.

        title, xlab, ylab
            Title, label at x and y axes of your figure(s). xtitle is defaulted
            to 'generation'.

        axes
            Whether or not plot axes. Default to ``True``.

        lty
            A list of line type for each line in the figure.

        col
            A list of colors for each line in the figure.

        level
            level of image colors (default to 20).

        saveAs
            save figures in files saveAs#gen.eps. For example, if saveAs='demo',
            you will get files demo1.eps, demo2.eps etc.

        separate
            plot data lines in separate panels.

        image
            use R image function to plot image, instead of lines.

        leaveOpen
            whether or not leave the plot open when plotting is done. Default to True.

        '''
        self.expr = expr
        if history == True:
			if plotType == "image":
				self.plotter = _VarPlotter_His_Image(numRep=numRep, varDim=varDim, win=win, ylim=ylim,
					update=update, title=title, xlab=xlab, ylab=ylab, axes=axes, mfrow=mfrow,
					byRep=byRep, byVal=byVal, separate=separate, plotType=plotType,
					level=level, saveAs=saveAs, lty=lty, col=col,
					leaveOpen=leaveOpen, dev=dev, width=width, height=height)
			else:
				self.plotter = _VarPlotter_His_Plot(numRep=numRep, varDim=varDim, win=win, ylim=ylim,
					update=update, title=title, xlab=xlab, ylab=ylab, axes=axes, mfrow=mfrow,
					byRep=byRep, byVal=byVal, separate=separate, plotType=plotType,
					level=level, saveAs=saveAs, lty=lty, col=col,
					leaveOpen=leaveOpen, dev=dev, width=width, height=height)
        else:
            if plotType == "image":
				self.plotter = _VarPlotter_NoHis_Image(numRep=numRep, win=win, ylim=ylim,
					update=update, title=title, xlab=xlab, ylab=ylab, axes=axes, mfrow=mfrow,
					byRep=byRep, plotType=plotType, level=level, saveAs=saveAs, lty=lty, col=col,
					leaveOpen=leaveOpen, dev=dev, width=width, height=height)
            else:
				self.plotter = _VarPlotter_NoHis_Plot(numRep=numRep, win=win, ylim=ylim,
					update=update, title=title, xlab=xlab, ylab=ylab, axes=axes, mfrow=mfrow,
					byRep=byRep, plotType=plotType, level=level, saveAs=saveAs, lty=lty, col=col,
					leaveOpen=leaveOpen, dev=dev, width=width, height=height)
        #
        # when apply is called, self.plotter.plot is called
        # *args, **kwargs may be things like rep=, gro=, at=, begin=
        pyOperator.__init__(self, func=self.plotter.plot,
            param=self.expr, *args, **kwargs)


