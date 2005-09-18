"""
simuPOP utilities.
This module provides some commonly used operators
and some utility class for plotting.
"""

from exceptions import *

try:
  from rpy import *
except:
  print "Rpy can not be loaded. Please verify your rpy installation. "
  print "If you are using a rpy version that reuires Python Numeric package. "
  print "Please install Numeric first. "
  raise exceptions.ImportError("Failed to import rpy package.")
  

import sys, math, random

from simuPOP import *
from simuUtil import Aggregator


def rmatrix(mat):
  ''' VERY IMPORTANT 
    convert a Python 2d list to r matrix format
    that can be passed to functions like image
    directly.
  '''
  return with_mode(NO_CONVERSION, r.do_call)('rbind', mat)
 

class VarPlotter_Base:
  def __init__(self, nplot, update, title, xlab, ylab, axes, 
    mfrow, plotType, saveAs, leaveOpen, dev='', width=0, height=0):
    """
    initialization function of base properties (layout etc)
    of all plotters
    """
    # save parameters
    self.nplot = nplot
    self.axes = axes
    self.mfrow = mfrow
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
        r('get(getOption("device"))()')
        self.device = r.dev_cur()
      self.first = False
    else:
      if self.dev == '':
        r.dev_set( self.device)

     
  def __del__(self) :
    if ( not self.leaveOpen) or self.dev != '':
      r.dev_off()

  def layout(self):
    # calculate layout
    if self.mfrow != [1,1] and self.mfrow[0]*self.mfrow[1] < self.nplot:
      raise ValueError("mfrow is not enough to hold " + str(self.nplot) + " figures")

    if self.mfrow == [1,1] and self.nplot > 1:
      self.mfrow[0] = int(math.ceil(math.sqrt( self.nplot )))
      self.mfrow[1] = int(math.ceil(self.nplot/float(self.mfrow[0]) ))
      if self.mfrow[0] > self.mfrow[1]:
           self.mfrow.reverse()
           
    # use par in simple case
    if self.plotType  == "plot":
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
        if  i*self.mfrow[1]+2+j <= self.nplot + 1:
          lt[i* (self.mfrow[1]+1)+j] = i*self.mfrow[1]+2+j
        else:
          lt[i* (self.mfrow[1]+1)+j] = 0  
      lt[(i+1)*(self.mfrow[1]+1)-1] = 1

    w = 7 * r.par("csi") * 2.54
    r.layout(r.matrix(lt, byrow=True, nc = self.mfrow[1]+1),
      widths = ["1"]*self.mfrow[1] + [ r.lcm(w)])
    
  def colorBar(self, level, lim):
    " draw a color bar to the right of the plot"
    mar = [0]*4
    mar[0] = [5, 1, 4, 2]
    
    self.color = r.rainbow(level, start=.7, end=.1)
    r.par(mar=mar)
    levels = r.seq(lim[0], lim[1], length=level)
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

 
class VarPlotter_NoHis(VarPlotter_Base):
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
     win:     generation window. If 0, plot all generation
     ylim:    limit of y axis. Will be the range of the first replicate if not given.
     update:  update figure after 'update' generations.
     title:   plot title. Will spread to other subplots if a vector is given.
     xlab:  title of x axis. Will spread to other subplots if a vector is given.
     ylab:  title of y axis. Will speread to other subplots if a vector is given.
     
              
     byRep:   Whether or not use different subplot for different replicate. Default 
              to false. If set to be true and there are more than one replicate, 
              you should set numRep and maybe mfrow.            
     numRep:  number of replicate. If byRep is true and there are more than 
              one replicates, set numRep here and varPlotter will try to set
              reasonaly subplot dimention, if dim is not given.
     mfrow:   if there are more than one subplots, mdrow=[x,y] specifies rows and
              columns of the subplots
  
  """
  
  def __init__(self, numRep=1, win=0, ylim=[0,0], update=1, plotType="plot",
               title="variable", xlab="generation", ylab="", axes=True, mfrow=[1,1],
               byRep=0, saveAs="", leaveOpen=True, level = 20, dev='', width=0, height=0):
    """
    
    """
    self.byRep = byRep
    
    if byRep == 1:
      VarPlotter_Base.__init__(self, numRep, update, title, xlab,
        ylab, axes, mfrow, plotType, saveAs, leaveOpen, dev, width, height)
    else:
      VarPlotter_Base.__init__(self, 1, update, title, xlab,
        ylab, axes, mfrow, plotType, saveAs, leaveOpen, dev, width, height)
       
    self.numRep = numRep
       
    if self.byRep == 1 and self.numRep > 1:
      for d in range(0, self.numRep):
        self.title[d] = self.title[d] + ", rep " + str(d)
             
    self.win = win
    self.ylim = ylim
    self.xlim = [0,0]
    if ylim == [0,0]:
      self.dynamicLim = 1
    else:
      self.dynamicLim = 0
    self.level = level  # for image
    
  def plot(self, gen, rep, data ):
    if type(data) == type(0) or type(data) == type(0.):
      data = [data]
  
    if rep >= self.numRep:
      raise ValueError("Replicate out of range. (forget to set numRep?)")
    
    # plot data
    if rep==-1 or rep == self.numRep - 1:
      self.CanUpdate()  # effectively increase skip
    if self.CanUpdate(False):  # check status, does not inc counter
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
      else:
        self.zlim=[data[0][0],data[0][0]]
        for i in range(1, len(data)):
          for j in range(1, len(data[i])):
            if data[i][j] < self.zlim[0]:
              self.zlim[0] = data[i][j]
            if data[i][j] > self.zlim[1]:
              self.zlim[1] = data[i][j]
  
      if self.plotType == 'plot':
        if self.byRep == 1 or rep <= 0 :
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
          d.extend( data[i])
        r.image( z= r.t(r.matrix( d, ncol=len(data))), 
          byrow=True, xlab=self.xlab, axes=self.axes,
          ylab=self.ylab, main=self.title[rep], col=self.color)
      else: 
        raise ValueError("Only 'plot' or 'image' is supported as plotType") 
           
    if self.saveAs != "":
      self.save(self.saveAs, gen)

# with history
class VarPlotter_His(VarPlotter_Base):
  """
  a class to plot variable, potentially from different replicate
  The idea is:
  
   when intiailize set some local variables
   each time during apply, call self.plot(gen, rep, var)
   var will be considered variable fromr eplicate rep.
 
   If the dimension of var is bigger than one, use par to create more than 
   one features. (has to be specified by dim.
  
   parameters:
     win:     generation window. If 0, plot all generation
     ylim:    limit of y axis. Will be the range of the first replicate if not given.
     update:  update figure after 'update' generations.
     title:   plot title. Will spread to other subplots if a vector is given.
     xlab:  title of x axis. Will spread to other subplots if a vector is given.
     ylab:  title of y axis. Will speread to other subplots if a vector is given.
     
              
     byRep:   Whether or not use different subplot for different replicate. Default 
              to false. If set to be true and there are more than one replicate, 
              you should set numRep and maybe mfrow.            
     numRep:  number of replicate. If byRep is true and there are more than 
              one replicates, set numRep here and varPlotter will try to set
              reasonaly subplot dimention, if dim is not given.
     byVal:   Whether or not use different subplot for each dimension of val. 
              Default to false. If set to be true and the dimension of data 
              is more than one, you should set varDim and maybe mfrgw.
     varDim:  dimension of variable. Need to be set if byVal is true.
     mfrow:   if there are more than one subplots, mdrow=[x,y] specifies rows and
              columns of the subplots

     separate: if true, use different panel (within a figure) to display each
              item of an array.

     plotType: usual plot or image
     level:    level of color, default to 20
  
  """
  
  def __init__(self, varDim=1, numRep=1, win=0, ylim=[0,0], update=1, 
               title="variable", xlab="generation", ylab="", axes=True, mfrow=[1,1],
               byRep=0, byVal=0, separate=False, plotType="plot", level=20,
               saveAs="", leaveOpen=True, dev='', width=0, height=0):
    """
    initialize _withHis plotter
    """    
    if byRep == 0 and numRep > 1 and byVal == 0 and varDim > 1:
      raise ValueError("Can not plot multi-dimension data for multiple replicate in one figure.\n" +
        "Use byRep or byVal variable.")
      
    if byVal == 1 and byRep == 1:
      raise ValueError("Can not separate variables both by replicate and by var dimension.")
      
    self.byRep = byRep
    self.byVal = byVal

    self.data = []
    nplot = 1
    if byRep == 1:
      self.numRep = numRep
      nplot = numRep
      for rep in range(0, numRep):
        self.data.append(Aggregator(win=win, width=varDim))
    elif byVal == 1:
      nplot = varDim
      for d in range(0, varDim):
        self.data.append(Aggregator(win=win,width=numRep))
    else:
      # otherwise, all data in one figure
      if numRep > 1:
        self.data.append(Aggregator(win, numRep))
        self.byVal = 1
      else:
        self.data.append(Aggregator(win, varDim))
        self.byRep = 1
  
    VarPlotter_Base.__init__(self, nplot, update, title,
      xlab, ylab, axes, mfrow, plotType, saveAs, leaveOpen, dev, width, height)

    self.varDim = varDim
    self.numRep = numRep
    self.separate = separate
    
    if len(self.data) > 1 :
      for d in range(0, len(self.data)):
        if byRep == 1:
          self.title[d] = self.title[d] + ", rep " + str(d)
        else:
          self.title[d] = self.title[d] + ", " + str(d)
    
    self.win = win
    self.ylim = ylim
    self.xlim = [0,0]
    self.level = level
    if ylim == [0,0]:
      self.dynamicLim = 1
    else:
      self.dynamicLim = 0
    self.first = True
    
  def clear(self):
    self.data=[]
    r.dev_off()

    
  def plot(self, gen, rep, data ):
    _data = data    # do not change data itself
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
      
    if self.byRep == 1: # push all data in to an replicate-specific aggregator
      if rep >= self.numRep:
        raise ValueError("Replicate number out of range. (Use numRep parameter.")
      self.data[rep].push( gen, _data)
    else: # separate data
      for d in range(0, self.varDim):
        self.data[d].push( gen, _data[d], rep)
    
    # plot data, only plot at the last replicate
    if rep==-1 or rep == self.numRep-1:
      if not self.CanUpdate() or len( self.data[0].gen) == 1:
        return
    else: 
      return

    self.setDev()
    self.layout()
    
    self.xlim[0] = self.data[0].gen[0]
    self.xlim[1] = self.data[0].gen[-1]

    if self.win > 0 and self.xlim[1] - self.xlim[0] < self.win:
      self.xlim[1] = self.xlim[0] + self.win

    if self.dynamicLim:
      self.ylim = self.data[0].ylim()
      for i in range(1, len(self.data)):
        yl = self.data[i].ylim()
        if yl[0] < self.ylim[0]:
          self.ylim[0] = yl[0]
        if yl[1] > self.ylim[1]:
          self.ylim[1] = yl[1]

    if self.byRep == 1:
      if self.plotType == "image":   # use image
        self.colorBar(self.level, self.ylim)
        for rep in range(0, self.numRep):
          r.image(x=self.data[rep].gen, y=range(0, self.data[rep].width),
            xlim=self.xlim, z= r.t(r.matrix( self.data[rep].flatData(),axes=self.axes,
            byrow=True, ncol=len(self.data[rep].gen))),xlab=self.xlab,
            ylab=self.ylab, main=self.title[rep], col=self.color)
      elif self.separate == False or self.data[rep].width == 1:
        for rep in range(0, self.numRep):
          r.plot( self.data[rep].gen, self.data[rep].data[0],
            xlab=self.xlab, ylab=self.ylab,axes=self.axes,
            main=self.title[rep], type='l', xlim=self.xlim,
            ylim=self.ylim, lty=1)
          for i in range( 1, self.data[rep].width ):        
            r.lines( self.data[rep].gen, self.data[rep].data[i], type='l', lty=i+1)
      else:  # use panel for each item of the plot
        height = self.ylim[1] - self.ylim[0]
        newlim = [self.ylim[0], self.ylim[0] + height* self.data[0].width]
        for rep in range(0, self.numRep):
          r.plot( self.data[rep].gen, self.data[rep].data[0],
            xlab=self.xlab, ylab=self.ylab,axes=self.axes,
            main=self.title[rep], type='l', xlim=self.xlim,
            ylim=newlim, lty=1)
          r.abline(h=[newlim[0] + height*x for x in range(0, self.data[rep].width)])
          for i in range( 1, self.data[rep].width ):        
            r.lines( self.data[rep].gen, 
              [height*i+x for x in self.data[rep].data[i]], type='l', lty=i+1)
         
          
    if self.byVal == 1:
      if self.plotType == "image":
        self.colorBar(self.level,self.ylim)
        for rep in range(0, self.numRep):
          r.image(x=self.data[rep].gen, y=range(0, self.data[rep].width), 
            xlim=self.xlim, z= r.t(r.matrix( self.data[rep].flatData(),axes=self.axes,
            byrow=True, ncol=len(self.data[rep].data[0]))),xlab=self.xlab,
            ylab=self.ylab, main=self.title[rep], col=self.color)
      elif self.separate == False or self.data[rep].width == 1:
        for v in range(0, self.varDim):
          r.plot( self.data[v].gen, self.data[v].data[0],
            xlab=self.xlab, ylab=self.ylab,axes=self.axes,
            main=self.title[v], type='l', xlim=self.xlim,
            ylim=self.ylim, lty=1)
          for i in range( 1, self.data[v].width ):        
            r.lines( self.data[v].gen, self.data[v].data[i], type='l', lty=i+1)
      else:
        height = self.ylim[1] - self.ylim[0]
        newlim = [self.ylim[0], self.ylim[0] + height* self.data[0].width]
        for v in range(0, self.varDim):
          r.plot( self.data[v].gen, self.data[v].data[0],
            xlab=self.xlab, ylab=self.ylab,axes=self.axes,
            main=self.title[v], type='l', xlim=self.xlim,
            ylim=newlim, lty=1)
          r.abline(h=[newlim[0] + height*x for x in range(0, self.data[rep].width)])
          for i in range( 1, self.data[v].width ):        
            r.lines( self.data[v].gen, 
              [height*i + x for x in self.data[v].data[i]], type='l', lty=i+1)
        
    self.save(gen)
          
# currently, type can be plot, image
#
def varPlotter(expr, history=1, varDim=1, numRep=1, win=0, ylim=[0,0],
  update=1, title="", xlab="generation", ylab="",axes=True, 
  mfrow=[1,1], separate=False, byRep=0, byVal=0, plotType="plot",
  level=20, saveAs="", leaveOpen=True, dev='', width=0, height=0, *args, **kwargs):
  # deal with additional arguments
  parm = ''
  for (k,v) in kwargs.items():
    parm += str(k) + '=' + str(v) + ', '
  #
  # pyExec( preStmts = 'rpv=VarPlotter(options....)',
  #   stmts = 'rpv.plot(gen,rep,expr', ... other options... )
  if title == "":
    title = expr
  # give a random name to the plotter
  plotterName = 'plt' + str( random.randint(0,1000000) )
  if history == 1: # with histroy, use VarPlotter_His
    cmd = '''pyEval(exposePop=1, name="""varPlotter %s""", %s \
      preStmts= \"\"\"%s=VarPlotter_His(varDim=%d, numRep=%d, win=%d, ylim=[%f,%f],  
      update=%d, title=\'\'\'%s\'\'\', xlab=\'\'\'%s\'\'\',  
      ylab=\'\'\'%s\'\'\', axes=%d,mfrow=[%d,%d], byRep=%d, 
      byVal=%d, separate=%d, plotType="%s", level=%d, saveAs=\'\'\'%s\'\'\', leaveOpen=%d, 
      dev="%s", width=%d, height=%d, 
      )\"\"\",
      postStmts=\"\"\"del %s\"\"\",
      stmts=r\"\"\"%s.plot(gen, rep, eval(\'\'\'%s\'\'\'))\"\"\")''' % \
    ( expr, parm, plotterName, varDim, numRep, win, ylim[0], ylim[1], update, title, \
    xlab, ylab, axes, mfrow[0], mfrow[1], byRep, byVal, separate, plotType, level, \
    saveAs, leaveOpen, dev, width, height, plotterName, plotterName, expr) 
  else: # without history, use VarPlotter_NoHis
    cmd = '''pyEval(exposePop=1, name="""varPlotter %s""", %s \
      preStmts= \"\"\"%s=VarPlotter_NoHis(numRep=%d, win=%d, ylim=[%f,%f], 
      update=%d, title=\'\'\'%s\'\'\', xlab=\'\'\'%s\'\'\',  
      ylab=\'\'\'%s\'\'\', axes=%d,mfrow=[%d,%d], byRep=%d, plotType="%s", level=%d,
      saveAs=\'\'\'%s\'\'\', leaveOpen=%d, dev="%s", width=%d, height=%d )\"\"\",
      postStmts= \"\"\"del %s\"\"\",
      stmts=r\"\"\"%s.plot(gen, rep, eval(\'\'\'%s\'\'\'))\"\"\")''' % \
    ( expr, parm, plotterName, numRep, win, ylim[0], ylim[1], update, title, \
    xlab, ylab, axes, mfrow[0], mfrow[1], byRep, plotType, level,\
    saveAs, leaveOpen, dev, width, height, plotterName, plotterName, expr) 
  #print cmd
  return eval(cmd)

# VarPlotter wrapper   
def VarPlot(value, history=0, varDim=1, numRep=1, win=0, ylim=[0,0],
  update=1, title="", xlab="generation", ylab="", axes=True,
  mfrow=[1,1], separate=False, byRep=0, byVal=0, plotType="plot",
  level=20, saveAs=""):
  # no history.
  VarPlotter_NoHis( numRep=1, win=0, ylim=[0,0],
    update=1, title=title, xlab=xlab, ylab=ylab, plotType=plotType, axes=axes,
    mfrow=[1,1], byRep=0, saveAs=saveAs, level = level).plot(0,0,value)


# 
# genomic contro, use whole subpopulation as case/control
# case: subpop index of the cases
# control: subpop indx of the controls.
# loci considered.
#
# To use this function, genoNum for 1-1, 1-2, 2-2 should
# exist. I.e., genoFreq[ [0 1 1 1 2 2 2] ] should
# be specified.
#
# the result is put in stat(rep).GC as a dictinary
# its keys are explained in R/gap/gcontrol .
# 
def GControl(pop, case=0, control=1, loci=[]):
  # first collect data, in the form of
  #
  # r0 r1 r2 s0 s1 s2
  #
  # for each loci.
  # here, we strictly use alleles  1 2
  # and treat 11 as r0, 12/21 as r1 and 22 as r2.
  #
  # so r0 is basically stat(rep,case).genoNum('1-1') and so on
  #
  #
  R = pop.dvars(case)
  S = pop.dvars(control)
  A = pop.dvars()
  def listGet(lst, idx, default={}):
    try:
      return lst[idx]
    except:
      return default
      
  if loci==[]:
    loci = range(0, pop.totNumLoci())
  data = []
  for loc in loci:
    data.extend( [R.genoNum[loc][1].get(1,0),
      R.genoNum[loc][1].get(2,0),
      listGet(R.genoNum[loc],2,{}).get(2,0),
      S.genoNum[loc][1].get(1,0),
      S.genoNum[loc][1].get(2,0),
      listGet(S.genoNum[loc],2,{}).get(2,0)])
  # now put data in R
  r.assign("data", data)
  r("""
    # import genetic analysis package library
    library(gap)
    data = matrix(data, ncol=6, byrow=TRUE)
    res = gcontrol(data)
    lambda = (median( res$x2 ) / 0.675)
    """)
  # get result
  A.GC = r('res')
  A.GC['lambda'] = r.get("lambda")
  return A.GC

# a wrapper of GControl
def gControl(case=0,control=1,loci=[],**kwargs):
  # deal with additional arguments
  parm = ''
  for (k,v) in kwargs.items():
    parm += str(k) + '=' + str(v) + ', '
  opt = ("""pyEval(%s stmts=r\'\'\'GControl(loci=""" + \
    str(loci) + """, case=%d, control=%d, rep=rep )\'\'\')""")\
    % ( parm, case, control) 
  #print opt
  return eval(opt)  
               
